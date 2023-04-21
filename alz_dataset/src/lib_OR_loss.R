
pacman::p_load(torch, luz)

nnf_ordinal_regression_loss <- function(pred,target,c=NULL){
  # combinations
  ij <- torch_combinations(torch_arange(1,length(pred)+1, dtype = torch_long()))
  # margin ranking loss
  l <- nnf_margin_ranking_loss(input1 = pred[ij[,1]]$squeeze(-1L),
                               input2= pred[ij[,2]]$squeeze(-1L),
                               target = sign(target[ij[,1]] - target[ij[,2]])$squeeze(-1L),
                               reduction = "none",margin = 1)

  if (!is.null(c)){
    # gathering the penalty
    p <- torch_diag(c[target[ij[,1]]$squeeze(-1L), target[ij[,2]]$squeeze(-1L)])
    # weighted loss
    loss <- l$dot(p$squeeze(-1L))
  } else {loss <- l$sum()}

}

nn_ordinal_regression_loss <- nn_module(
  "nn_ordinal_regression_loss",
  initialize = function(c=NULL){
    self$c <- c
  },
  forward = function(input, target){
    nnf_ordinal_regression_loss(input, target, c = self$c)
  }

)

