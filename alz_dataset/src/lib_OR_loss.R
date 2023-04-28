
pacman::p_load(torch, luz)

# nnf_ordinal_regression_loss <- function(pred,target,C=NULL){
#   # combinations
#   ij <- torch_combinations(torch_arange(1,length(pred)+1, dtype = torch_long()))
#   # margin ranking loss
#   l <- nnf_margin_ranking_loss(input1 = pred[ij[,1]]$squeeze(-1L),
#                                input2= pred[ij[,2]]$squeeze(-1L),
#                                target = sign(target[ij[,1]] - target[ij[,2]])$squeeze(-1L),
#                                reduction = "none",margin = 1)
#
#   if (!is.null(C)){
#     # gathering the penalty
#     #browser()
#     p <- torch_diag(c[target[ij[,1]]$squeeze(-1L), target[ij[,2]]$squeeze(-1L)])
#     # weighted loss
#     loss <- l$dot(p$squeeze(-1L))
#   } else {loss <- l$sum()}
#
# }

#' Ordinal regression loss
#' @param C a square 2D float tensor representing the cost matrix
#' @param margin a scalar positive value
nnf_ordinal_regression_loss <- function(pred,target,C,margin=1.0) {
  pred <- pred$squeeze(-1L)
  target <- target$squeeze(-1L)
  stopifnot(target$ndim == 1L)
  stopifnot(pred$ndim == 1L)
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
  #browser()
  l$dot(w) / torch_maximum(sum(w!=0),torch_ones(1)$to(device=l$device))
}

nn_ordinal_regression_loss <- nn_module(
  "nn_ordinal_regression_loss",
  initialize = function(C=NULL){
    self$C <- C
  },
  forward = function(input, target){
    nnf_ordinal_regression_loss(input, target, C = self$C)
  }

)

