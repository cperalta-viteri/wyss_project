



library(torch)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# A linear layer with rescaled output and pruning options
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
nn_linear_rescale <- nn_module(
  "linear_rescale",
  
  initialize = function(in_features,out_features,prune_mask=NULL) {
    self$linear <- nn_linear(in_features = in_features, out_features = out_features,bias = FALSE)
    if (is.null(prune_mask)) {
      prune_mask <- torch_ones_like(self$linear$weight,torch_bool())
    }
    stopifnot(identical(dim(self$linear$weight),dim(prune_mask)))
    self$prune_mask <- nn_buffer(prune_mask)
  },
  
  forward = function(x) {
    # apply pruning mask and normalize weights to 1
    with_no_grad({
      self$linear$weight$mul_(self$prune_mask)
      d <- self$linear$weight$norm(dim=2L,keepdim=TRUE)
      d$clamp_min_(1e-6)
      self$linear$weight$div_(d)
    })
    
    # compute normalization factors for each node of layer1 and each instance
    sf <- torch_matmul(x,self$linear$weight$abs()$transpose(-2L,-1L))
    sf <- sf$clamp_min_(0.01)
    
    # Compute layer output
    self$linear(x)/sf
  },
  
  compute_prune_mask = function(n=50L) {
    lo <- self$linear$weight$topk(n,largest = FALSE)[[1]][,n,drop=FALSE]
    hi <- self$linear$weight$topk(n,largest = TRUE)[[1]][,n,drop=FALSE]
    self$prune_mask$set_((self$linear$weight<=lo) | (self$linear$weight>=hi))
    self$prune_mask
  }
)



#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Define an advanced model with layer1 scaling and pruning
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
advanced_cell_classifier <- nn_module(
  "advanced_cell_classifier",
  
  initialize = function(feature_names,class_names,...) {
    self$feature_names <- feature_names
    self$class_names <- class_names
    self$input_dropout <- nn_dropout(0.5)
    self$layer1 <- nn_linear_rescale(in_features = length(self$feature_names),out_features = 4096L,...)
    self$layers <- nn_sequential(
      nn_tanh(),
      nn_dropout(0.3), nn_linear(in_features = 4096L, out_features = 1024L), nn_relu(),
      nn_dropout(0.3), nn_linear(in_features = 1024L, out_features = 512L), nn_relu(),
      nn_dropout(0.1), nn_linear(in_features = 512L, out_features = 256L), nn_relu(),
      nn_dropout(0.1), nn_linear(in_features = 256L, out_features = length(self$class_names))
    )
  },
  forward = function(x) {x |> self$input_dropout() |> self$layer1() |> self$layers()},
  loss = function(input,target) {nnf_multi_margin_loss(input,target)},
  set_optimizers = function(lr1=1e-5,lr2=lr1) {list(
  	  	l1 = optim_adam(self$layer1$parameters, lr = lr1,weight_decay=0.1),
  		  l2 = optim_adam(self$layers$parameters, lr = lr2)
  )},
  compute_prune_mask = function(...) {self$layer1$compute_prune_mask()}
)



#' Predict cell type using given torch model
#' @param N expression matrix, ideally log1p(RPM/100)/log(2) counts, but raw counts seems to work too
#' @param model_file path to the file containng the model
#' @param ignore.case if TRUE ignore case when aligning feature names (=gene names). It is useful for easy Human-Mouse gene-name matching.
#' @return a prediction matrix giving each input cell a score for each cell type
torch_predict_cell_type <- function(N,model_file,ignore.case=FALSE) {
  # Load the model
  model2 <- torch_load(model_file,device = "cpu")
  
  # Align the input matrix to desired matrix format, adding 0 when genes are missing
  if (isTRUE(ignore.case)) {
    i <- match(toupper(model2$feature_names),toupper(rownames(N)),nomatch = 0)  
  } else {
    i <- match(model2$feature_names,rownames(N),nomatch = 0)  
  }
  stopifnot(mean(i>0)>=0.2) # check that 20% of model genes are in the input matrix
  m <- rbind(N[1,,drop=FALSE]*0,N)[i+1,]
  
  # Compute prediction by batch of 128 cells
  f <- factor((seq(ncol(m))-1L) %/% 128L)
  P <- lapply(levels(f),function(i) {
    model2(torch_tensor(as.matrix(t(m[,f==i,drop=FALSE]))))
  })
  
  # Concatenate prediction and return
  P <- torch_cat(P)
  P <- as_array(P)
  colnames(P) <- model2$class_names
  P  
}



library(DelayedArray)
sc_dataset <- dataset(
  name = "DelayedArray_dataset",
  initialize = function(x,y) {
      self$x <- x
      self$y <- y
      stopifnot(identical(length(y),nrow(x)))
  },
  .length = function() {length(self$y)},
  .getbatch = function(index) {list(
        x = torch_tensor(as.matrix(self$x[index,,drop=FALSE]),torch_float()),
        y = torch_tensor(self$y[index],torch_long())
  )},
  .getitem = function(index) {
      if (is.list(index)) {index <- unlist(index)}
      self$.getbatch(index)
  }
)




