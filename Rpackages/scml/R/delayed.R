

#' @title delayed_dataset
#' @description A torch dataset which accepts DelayedArrays.
#' @param x Count matrix.
#' @param y Cell type array.
#' @param dtypes torch types.
#' @param ...  ###łæßðđðsdfgw
#' @export
#' @import torch
#' @import DelayedArray
x_y_delayed_dataset <- torch::dataset(
  name = "delayed_dataset",
  initialize = function(x,y,dtypes=list(torch_float())) {
      self$x <- x ; self$y <- y
      self$dtypes <- dtypes
      stopifnot(identical(length(y), nrow(x)))
  },
  .length = function() {nrow(self$arrays[[1]])},
  .getbatch = function(index) {
      # mapply(function(A,dtype) {
      #   torch_tensor(as.matrix(A[index,,drop=TRUE]),dtype)
      # },self$arrays,self$dtypes)
    list(x=torch_tensor(as.matrix(self$x[index,]), self$dtypes[[1]]),
         y=torch_tensor(as.matrix(self$y[index]), self$dtypes[[2]])$flatten()) ##loss function does not accept 2D target array
  },
  .getitem = function(index) {
      if (is.list(index)) {index <- unlist(index)}
      self$.getbatch(index)
  }
)



# NN, optimizer=optim_adam, metrics=list(), loss, hparams=list(), epoch,
#                      train_data, test_data, dataloader_options=list(), callbacks=list()

#' @title train_torch_NN
#' @description Allows user-friendly training set up of a torch model
#' @param NN Torch module.
#' @param optimizer Torch optimizer, by default set to `optim_adam`.
#' @param metrics Luz metrics to evaluate model performance. See Luz documentation for reference.
#' @param loss Torch loss.
#' @param hparams Module's parameters, passed as list: `hparams=list(feature_names=rownames(train), class_names=train$subclass_label)`.
#' @param weight_decay Weight decay parameter for optimizer. Default set to 0.
#' @param lr Learning rate for optimizer. Default set to 100.
#' @param epoch Number of epochs. Default set to 100.
#' @param train_data Dataset used for fitting the model.
#' @param test_data Dataset used model performance assessment.
#' @param dataloader_options Options to be passed to `torch::dataloader()` as list. By default it is setting the batch size to 100:
#' `dataloader_options=list(batch_size=100)`.
#' @param callbacks Argument list containing Luz callbacks.
#' @export

train_torch_NN <- function(NN, optimizer=optim_adam, metrics=NULL, loss=NULL, hparams=NULL,
                           weight_decay=0, lr=0.001, epoch=100L,
                           train_data, test_data, dataloader_options=list(batch_size=100), callbacks=list()){
  NN |>
    setup(optimizer=optim_adam,
          metrics = metrics,
          loss=loss
    ) |>
    set_hparams(unlist(hparams)) |>
    set_opt_hparams(weight_decay=weight_decay, lr=lr) |>
    fit(
      epoch = epoch,
      data = train_data,
      valid_data = test_data,
      dataloader_options=dataloader_options,
      accelerator = accelerator(cpu=TRUE),
      callbacks = callbacks
    )
}



#' @title predict_delayed
#' @description Predict cell type given a torch model
#' @param N expression matrix, ideally log1p(RPM/100)/log(2) counts, but raw counts seems to work too
# @param model_file path to the file containing the model
#' @param luz_model model to be fitted
#' @param ignore.case if TRUE ignore case when aligning feature names (=gene names). It is useful for easy Human-Mouse gene-name matching.
#' @return a prediction matrix giving each input cell a score for each cell type
#' @export
#' @import torch
#' @import luz
#' @import stats
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
  pred <- stats::predict(
    luz_model,
    delayed_dataset(t(m),rep(1L,ncol(m))),
    dataloader_options = list(batch_size=4096L)
  )
  pred
}



