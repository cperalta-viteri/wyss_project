



trn_ds <- tensor_dataset(
  x = torch_tensor(data.matrix(iris[1:4])),
  y = torch_tensor(as.integer(iris$Species))
)



#' Ordinal regression loss
#' @param C a square 2D float tensor representing the cost matrix
#' @param margin a scalar positive value
nnf_ordinal_regression_loss <- function(pred,target,C,margin=1.0) {
  pred <- pred$squeeze(-1L)
  target <- target$squeeze(-1L)
  stopifnot(pred$ndim == 1L)
  stopifnot(target$ndim == 1L)
  stopifnot(ncol(C) == nrow(C))

  ij <- torch_combinations(torch_arange(1L,length(pred)+1L, dtype = torch_long()))
  y <- target[ij]
  l <- nnf_margin_ranking_loss(
    input1 = pred[ij[,1L]], 
    input2 = pred[ij[,2L]], 
    target = sign(y[,1L] - y[,2L]),
    reduction = "none",margin = margin
  )
  w <- torch_gather(C[y[,1L],],2L,y[,2L,drop=FALSE])$squeeze(-1L)
  sum(l*w) / sum(w!=0)
}


nn_ordinal_regression_loss <- nn_module(
  "nn_ordinal_regression_loss",
  initialize = function(C) {
    self$C <- C
  },
  forward = function(input, target){
    nnf_ordinal_regression_loss(input, target, C = self$C)
  }
)


ordinal_regressor <- nn_module(
  "ordinal_regressor",
  initialize = function() {
    self$layers <- nn_sequential(
      nn_dropout(),
      nn_linear(in_features = 4L, out_features = 50L), nn_tanh(),
      nn_linear(in_features = 50L, out_features = 50L), nn_relu(),
      nn_linear(in_features = 50L, out_features = 50L), nn_relu(),
      nn_linear(in_features = 50L, out_features = 1L)
    ) 
  },
  forward = function(x) {
    x |> self$layers()
  }
)

C <- 1 - torch_eye(nlevels(iris$Species))
fm <- ordinal_regressor |>
  setup(loss = nn_ordinal_regression_loss(C=C),optimizer = optim_adam) |>
  set_opt_hparams(lr=1e-3,weight_decay=0) |>
  fit(trn_ds,
      epochs = 200L,
      dataloader_options = list(batch_size=32L,drop_last=FALSE),
      accelerator = accelerator(cpu = TRUE)
  )

boxplot(as_array(predict(fm,trn_ds)$to(device="cpu"))~iris$Species)
plot(fm)




