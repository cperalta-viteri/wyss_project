
#' @title nnf_ordinal_regression_loss
#' @description  Ordinal regression loss function.
#' @param pred Predicted values.
#' @param target Target values.
#' @param C a square 2D float tensor representing the cost matrix.
#' @param margin a scalar positive value.
#' @return loss value
#' @export
#' @import torch
nnf_ordinal_regression_loss <- function(pred,target,C,margin=1.0) {
  pred <- pred$flatten()
  target <- target$flatten()
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

#' @title nn_ordinal_regression_loss
#' @description  Ordinal regression loss module
#' @param C a square 2D float tensor representing the cost matrix.
#' @export
#' @import torch
#' @examples
#' nn_ordinal_regression_loss(1 - torch::torch_eye(3L))
nn_ordinal_regression_loss <- torch::nn_module(
  "nn_ordinal_regression_loss",
  initialize = function(C) {
    self$C <- C
  },
  forward = function(input, target){
    nnf_ordinal_regression_loss(input, target, C = self$C)
  }
)
