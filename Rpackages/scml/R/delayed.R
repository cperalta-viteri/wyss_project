


#' A torch dataset which accepts DelayedArrays
#' @export
#' @import torch
#' @import DelayedArray
delayed_dataset <- torch::dataset(
  name = "delayed_dataset",
  initialize = function(...,dtypes=list(torch_float())) {
      self$arrays <- list(...)
      self$dtypes <- dtypes
  },
  .length = function() {nrow(self$arrays[[1]])},
  .getbatch = function(index) {
      mapply(function(A,dtype) {
        torch_tensor(as.matrix(A[index,,drop=FALSE]),dtype)
      },self$arrays,self$dtypes)
  },
  .getitem = function(index) {
      if (is.list(index)) {index <- unlist(index)}
      self$.getbatch(index)
  }
)



#' Predict cell type using given torch model
#' @param N expression matrix, ideally log1p(RPM/100)/log(2) counts, but raw counts seems to work too
#' @param model_file path to the file containng the model
#' @param ignore.case if TRUE ignore case when aligning feature names (=gene names). It is useful for easy Human-Mouse gene-name matching.
#' @return a prediction matrix giving each input cell a score for each cell type
#' @export
#' @import torch
#' @import luz
predict_delayed <- function(N,luz_model,ignore.case=FALSE) {
  # Align the input matrix to desired matrix format, adding 0 when genes are missing
  if (isTRUE(ignore.case)) {
    i <- match(toupper(luz_model$model$feature_names),toupper(rownames(N)),nomatch = 0)
  } else {
    i <- match(luz_model$model$feature_names,rownames(N),nomatch = 0)
  }
  stopifnot(mean(i>0)>=0.2) # check that 20% of model genes are in the input matrix
  m <- rbind(N[1,,drop=FALSE]*0,N)[i+1,]

  # Compute prediction by batch
  pred <- predict(
    luz_model,
    delayed_dataset(t(m),rep(1L,ncol(m))),
    dataloader_options = list(batch_size=4096L)
  )
  pred
}



