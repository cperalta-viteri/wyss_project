

#' @title delayed_dataset
#' @description A torch dataset which accepts DelayedArrays.
#' @param x Count matrix.
#' @param y Cell type array.
#' @param dtypes torch types.
#' @param ...  ###łæßðđðsdfgw
#' @export
#' @import torch
#' @import DelayedArray
delayed_xy_dataset <- torch::dataset(
  name = "delayed_dataset",
  initialize = function(x,y=rep(0L,nrow(m)),x_dtype=torch_float(),y_dtype=torch_long()) {
    stopifnot(identical(length(y), nrow(x)))
    stopifnot(identical(length(dim(x)),2L))
    self$x <- x; self$y <- y
    self$x_dtype <- x_dtype; self$y_dtype <- y_dtype
  },
  .length = function() {length(self$y)},
  .getbatch = function(index) {
    list(
      x = torch_tensor(as.matrix(self$x[index,,drop=FALSE]), self$x_dtype),
      y = torch_tensor(self$y[index], self$y_dtype)
    )
  },
  .getitem = function(index) {
      if (is.list(index)) {index <- unlist(index)}
      self$.getbatch(index)
  }
)



#' @title predict_delayed
#' @description Predict cell type given a torch model
#' @param N expression matrix, ideally log1p(RPM/100)/log(2) counts, but raw counts seems to work too
#' @param luz_model model to be fitted
#' @param ignore.case if TRUE ignore case when aligning feature names (=gene names). It is useful for easy Human-Mouse gene-name matching.
#' @return a prediction matrix giving each input cell a score for each cell type
#' @export
#' @importFrom stats predict
predict_delayed <- function(N,luz_model,ignore.case=FALSE,batch_size=1024L,...) {
  # Align the input matrix to desired matrix format, adding 0 when genes are missing
  if (isTRUE(ignore.case)) {
    i <- match(toupper(luz_model$model$feature_names),toupper(colnames(N)),nomatch = 0)
  } else {
    i <- match(luz_model$model$feature_names,colnames(N),nomatch = 0)
  }
  stopifnot(mean(i>0)>=0.2) # check that 20% of model genes are in the input matrix
  m <- cbind(N[,1,drop=FALSE]*0,N)[,i+1]

  # Compute prediction by batch
  pred <- predict(
    luz_model,
    delayed_xy_dataset(m),
    dataloader_options = list(batch_size=batch_size),
    ...
  )
  pred
}





