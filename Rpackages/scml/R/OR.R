

#' Ordinal regression loss
#' @param C a square 2D float tensor representing the cost matrix
#' @param margin a scalar positive value
#' @return loss value
#' @export
#' @import torch
nnf_ordinal_regression_loss <- function(pred,target,C,margin=1.0) {
  pred <- pred$squeeze(-1L)
  target <- target$squeeze(-1L)
  stopifnot(pred$ndim == 1L)
  stopifnot(target$ndim == 1L)
  stopifnot(identical(ncol(C),nrow(C)))

  ij <- torch_combinations(torch_arange(1L,length(pred)+1L, dtype = torch_long()))
  y <- target[ij]
  l <- nnf_margin_ranking_loss(
    input1 = pred[ij[,1L]],
    input2 = pred[ij[,2L]],
    target = sign(y[,1L] - y[,2L]),
    reduction = "none",margin = margin
  )
  w <- torch_gather(C[y[,1L],],2L,y[,2L,drop=FALSE])$squeeze(-1L)
  l$dot(w) / torch_maximum(sum(w!=0),1L)
}

#' @export
#' @import torch
nn_ordinal_regression_loss <- torch::nn_module(
  "nn_ordinal_regression_loss",
  initialize = function(C) {
    self$C <- C
  },
  forward = function(input, target){
    nnf_ordinal_regression_loss(input, target, C = self$C)
  }
)
