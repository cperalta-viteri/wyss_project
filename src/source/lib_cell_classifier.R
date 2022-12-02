

source("/home/dnf/Desktop/final_model_julien_torch/allen20/src/julien/draft/torch/lib_torch.R")


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Define the model
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
cell_classifier <- nn_module(
  "cell_classifier",
  initialize = function(feature_names,class_names) {
    self$feature_names <- feature_names
    self$class_names <- class_names

    #nn_batch_norm1d(length(feature_names))

    self$layer1_dropout <- nn_dropout(0.5)
    self$layer1_linear <- nn_linear(in_features = length(feature_names), out_features = 1024L,bias = FALSE)
    self$prune_mask <- nn_buffer(torch_ones_like(self$layer1_linear$weight,torch_bool()))
    
    self$layer2 <- nn_sequential(
      nn_dropout(0.2), nn_linear(in_features = 1024L, out_features = 512L), nn_relu(),
      nn_dropout(0.2), nn_linear(in_features = 512L, out_features = 256L), nn_relu(),
      nn_dropout(0.2), nn_linear(in_features = 256L, out_features = length(class_names))
    )
    
    
  },
  
  compute_prune_mask = function(n=50L) {
    lo <- self$layer1_linear$weight$topk(n,largest = FALSE)[[1]][,n,drop=FALSE]
    hi <- self$layer1_linear$weight$topk(n,largest = TRUE)[[1]][,n,drop=FALSE]
    self$prune_mask$set_((self$layer1_linear$weight<=lo) | (self$layer1_linear$weight>=hi))
  },
  
  forward = function(rpm) {
    with_no_grad({
      self$layer1_linear$weight$mul_(self$prune_mask)
      self$layer1_linear$weight$div_(self$layer1_linear$weight$norm(dim=2L,keepdim=TRUE))
    })
    
    # preprocessing
    #rpm <- trn_ds[1:10]$x
    log_rpm <- rpm$log1p()

    sf <- torch_matmul(log_rpm,self$layer1_linear$weight$abs()$transpose(1,2))
    sf <- sf$clamp(0.01)
    
    log_rpm %>% 
      self$layer1_dropout() %>%
      self$layer1_linear() %>%
      torch_div(sf) %>%
      nnf_hardtanh() %>%
      self$layer2()
  },
  
  # calculate the loss
  criterion = function(y_pred,y_true) {
    loss <- nnf_multi_margin_loss(y_pred, y_true) # calculate loss
    loss
  },
  
  # Fit the model on given training data
  fit = function(trn_dl,val_dl=NULL,num_epoch=30L,device='cpu') {
    self$to(device=device)
    
    log <- epoch_logger$new()  
    optimizer <- optim_adam(self$parameters,lr=0.0001,weight_decay = 0)
    for (epoch in seq(num_epoch)) {
      #-#-#-#-#-#-#-#-#-#-#
      # Epoch training
      #-#-#-#-#-#-#-#-#-#-#
      self$train() # activate train mode
      trn_log <- batch_logger$new(sprintf("Epoch %d:trn",log$epoch))
      coro::loop(for (b in trn_dl) {
        optimizer$zero_grad()
        y_pred <- self$forward(b$x$to(device=device))
        loss <- self$criterion(y_pred, b$y$to(device=device))
        loss$backward()   # have gradients get calculated
        optimizer$step()  # have gradients get applied
        trn_log$add(loss$item(),as_array(y_pred$argmax(2L)$to(device="cpu")),as_array(b$y$to(device="cpu")))
      })
      log$add("loss",mean(trn_log$losses),"trn")
      log$add("acc",trn_log$accuracy(),"trn")
      cat("\n")
      gc()
      
      #-#-#-#-#-#-#-#-#-#-#
      # Epoch validation
      #-#-#-#-#-#-#-#-#-#-#
      self$eval() # activate eval mode
      val_log <- batch_logger$new(sprintf("Epoch %d:val",log$epoch))
      if (!is.null(val_dl)) {
        coro::loop(for (b in val_dl) {
          y_pred <- self$forward(b$x$to(device=device))
          loss <- self$criterion(y_pred, b$y$to(device=device))
          val_log$add(loss$item(),as_array(y_pred$argmax(2L)$to(device="cpu")),as_array(b$y$to(device="cpu")))
        })
        log$add("loss",mean(val_log$losses),"val")
        log$add("acc",val_log$accuracy(),"val")
        cat("\n")
      }
      
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      # Epoch end
      #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
      #torch_save(self,"out/last_model.torch")
      #saveRDS(log,file="out/last_log.rds")
      #plot(log)
      log$end_of_epoch()
    }
    list(log=log,trn_log=trn_log,val_log=val_log)
  }
)


