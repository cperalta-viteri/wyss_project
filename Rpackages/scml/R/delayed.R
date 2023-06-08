

#' @title delayed_dataset
#' @description A torch dataset which accepts DelayedArrays.
#' @param x Count matrix.
#' @param y Cell type array.
#' @param dtypes torch types.
#' @export
#' @import torch
#' @import DelayedArray
delayed_classifier_dataset <- torch::dataset(
  name = "delayed_dataset",
  initialize = function(x, y = rep(0L, nrow(x)), x_dtype = torch_float(), y_dtype = torch_long()) {
    stopifnot(identical(length(y), nrow(x)))
    stopifnot(identical(length(dim(x)), 2L))
    self$x <- x
    self$y <- y
    self$x_dtype <- x_dtype
    self$y_dtype <- y_dtype
  },
  .length = function() {
    length(self$y)
  },
  .getbatch = function(index) {
    list(
      x = torch_tensor(as.matrix(self$x[index, , drop = FALSE]), self$x_dtype),
      y = torch_tensor(self$y[index], self$y_dtype)
    )
  },
  .getitem = function(index) {
    if (is.list(index)) {
      index <- unlist(index)
    }
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
predict_delayed <- function(N, luz_model, ignore.case = FALSE, batch_size = 1024L, ...) {
  # Align the input matrix to desired matrix format, adding 0 when genes are missing
  if (isTRUE(ignore.case)) {
    i <- match(toupper(luz_model$model$feature_names), toupper(colnames(N)), nomatch = 0)
  } else {
    i <- match(luz_model$model$feature_names, colnames(N), nomatch = 0)
  }
  stopifnot(mean(i > 0) >= 0.2) # check that 20% of model genes are in the input matrix
  m <- cbind(N[, 1, drop = FALSE] * 0, N)[, i + 1]

  # Compute prediction by batch
  pred <- predict(
    luz_model,
    delayed_classifier_dataset(m),
    dataloader_options = list(batch_size = batch_size),
    ...
  )
  pred
}





#' @title train_torch_NN
#' @description Allows user-friendly training set up of a torch model.
#' @param x Count matrix.
#' @param y Array containing the cell type class label.
#' @param epoch Total number of epoch to train
#' @param prune_at_epoch The epoch at which we prune the model. Set to NA to train without pruning.
#' @param batch_size Size of batch of data loaded at each training step.
#' @param prune_size number of top/down gene to keep for each node of layer1 when pruning
#' @param lr a numeric vector of decreasing learning rates
#' @param auto_resume_path a character to specify a file to store weights at each epoch,
#'                         so training can be resumed when stoped. Set to NA to disable.
#' @param ... additional parameters are passed to fit() call
#' @export
train_delayed_classifier <- function(x, y, epoch = 100L, prune_at_epoch = epoch %/% 2L, batch_size = 128L, prune_size = pmin(50L, ncol(x)), lr = 10^(-3:-6), auto_resume_path = NA, ...) {
  # stopifnot(epoch>=2L)
  if (!identical(lr, sort(lr, decreasing = TRUE))) stop("lr must be a decreasing numeric vector")
  y <- as.factor(y)


  # Define the function which compute lr from a given epoch number
  lr_fun <- if (!is.na(prune_at_epoch)) {
    stopifnot(prune_at_epoch <= epoch)
    stepfun(
      c(seq(1L, prune_at_epoch, length.out = length(lr) + 1L), seq(prune_at_epoch + 1L, max(epoch, prune_at_epoch + 1L), length.out = length(lr) + 1L)),
      c(max(lr), lr, min(lr), lr, min(lr))
    )
  } else {
    stepfun(
      seq(1L, epoch, length.out = length(lr) + 1L),
      c(max(lr), lr, min(lr))
    )
  }
  # plot(lr_fun(seq(epoch)));abline(v=prune_at_epoch)

  fm <- nn_cell_scorer |>
    setup(
      optimizer = optim_adam,
      metrics = list(luz_metric_accuracy()),
      loss = nn_multi_margin_loss(margin = 1)
    ) |>
    set_hparams(feature_names = colnames(x), class_names = levels(y)) |>
    set_opt_hparams(weight_decay = 0, lr = 1) |>
    fit(
      epoch = epoch,
      data = delayed_classifier_dataset(x, as.integer(y)),
      dataloader_options = list(batch_size = batch_size, drop_last = FALSE, shuffle = TRUE),
      ...,
      callbacks = c(
        list(
          luz_callback_lr_scheduler(lr_lambda, verbose = TRUE, lr_fun)
        ),
        if (is.na(prune_at_epoch)) NULL else luz_callback_prune("layer1.linear.weight", n = prune_size, at_epoch = prune_at_epoch),
        if (is.na(auto_resume_path)) NULL else luz_callback_auto_resume(path = auto_resume_path)
      )
    )
  fm
}
