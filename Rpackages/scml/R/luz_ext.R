

#' @title luz_callback_prune
#' @param param Model parameters
#' @param n Nodes
#' @param at_epoch Epoch at which compute pruning mask.
#' @import torch
#' @import luz
#' @export
luz_callback_prune <- luz::luz_callback(
  initialize = function(param,n,at_epoch) {
    self$param <- param
    self$n <- n
    self$at_epoch <- at_epoch
  },
  on_train_batch_after_step = function() {
    if (ctx$epoch >= self$at_epoch) {
      w <- ctx$model$parameters[[self$param]]
      if (is.null(self$prune_mask)) {
        lo <- w$topk(self$n,largest = FALSE)[[1]][,self$n,drop=FALSE]
        hi <- w$topk(self$n,largest = TRUE)[[1]][,self$n,drop=FALSE]
        self$prune_mask <- (w>lo) & (w<hi)
      }
      w[self$prune_mask] <- 0
    }
  }
)
