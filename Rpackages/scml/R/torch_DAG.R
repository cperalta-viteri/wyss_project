
#' @title nnf_dag_loss
#' @description Direct Acyclic Graph loss function.
#' @param pred Predicted values.
#' @param target Target values.
#' @param e a 2 column integer matrix (edge list of the DAG), with parent index
#'          in the first column and child index in the second
#' @param w weight vector for each node of the dag (how much to penalize a node)
#' @param margin a scalar positive value.
#' @return loss value
#' @export
#' @import torch
nnf_dag_loss <- function(pred,target,e,w,margin=1.0) {
  e <- torch_tensor(e)
  w <- torch_tensor(w)
  stopifnot(ncol(e)==2L)
  y2 <- target$index_select(-1L,e[,2L])  # y_true of the child
  f2 <- pred$index_select(-1L,e[,2L])  # child value
  ub <- torch_full_like(pred,+Inf)$index_reduce(-1L,e[,1L],f2$masked_fill(!y2,+Inf),"amin",include_self = FALSE)
  lb <- torch_full_like(pred,-Inf)$index_reduce(-1L,e[,1L],f2$masked_fill(y2,-Inf),"amax",include_self = FALSE)
  l <- (lb + margin - ub)$maximum(0)

  #nub <- torch_zeros_like(pred)$index_add(-1L,e[,1L],y2*w[e[,2L]])
  #nlb <- torch_zeros_like(pred)$index_add(-1L,e[,1L],(!y2)*w[e[,2L]])

  l$mv(w)$mean()
}


#' @title nn_dag_loss
#' @description Direct Acyclic Graph loss module
#' @param e a 2 column integer matrix (edge list of the DAG), with parent index
#'          in the first column and child index in the second
#' @param w weight vector for each node of the dag (how much to penalize a node)
#' @export
#' @import torch
nn_dag_loss <- torch::nn_module(
  "nn_ordinal_regression_loss",
  initialize = function(e,w) {
    self$e <- e
    self$w <- w
  },
  forward = function(input, target){
    nnf_dag_loss(input, target, e = self$e, w = self$w)
  }
)
