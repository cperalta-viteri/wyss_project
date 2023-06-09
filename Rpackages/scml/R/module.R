

#' @title nn_linear_rescale
#' @description  A linear layer a rescaled output.
#' @param in_features Input number of features.
#' @param out_features Output number of features
#' @import torch
#' @export
nn_linear_rescale <- nn_module(
  initialize = function(in_features,out_features) {
    self$linear <- nn_linear(in_features = in_features, out_features = out_features,bias = FALSE)
  },
  forward = function(x) {
    # compute normalization factors for each node of layer1 and each instance
    sf <- torch_matmul(x$abs(),self$linear$weight$abs()$transpose(-2L,-1L))
    sf <- sf$clamp_(0.01,100) # maximal and minimal rescaling factor
    self$linear(x)/sf
  }
)


#' @title nn_cell_scorer
#' @description Cell type classifier NN.
#' @param feature_names Names of the features.
#' @param class_names Class of the cells
#' @param ... additional argument are forwarded to nn_linear_rescale()
#' @import torch
#' @export
nn_cell_scorer <- torch::nn_module(
  initialize = function(feature_names,class_names,...) {
    self$feature_names <- feature_names
    self$class_names <- class_names
    self$input_dropout <- nn_dropout(0.5)
    self$layer1 <- nn_linear_rescale(in_features = length(self$feature_names),out_features = 4096L,...)
    self$layers <- nn_sequential(
      nn_tanh(),
      nn_dropout(0.3), nn_linear(in_features = 4096L, out_features = 1024L), nn_relu(),
      nn_dropout(0.3), nn_linear(in_features = 1024L, out_features = 512L), nn_relu(),
      nn_dropout(0.1), nn_linear(in_features = 512L, out_features = 256L), nn_relu(),
      nn_dropout(0.1), nn_linear(in_features = 256L, out_features = length(self$class_names))
    )
  },
  forward = function(x) {x |> self$input_dropout() |> self$layer1() |> self$layers()},
)


#' @title train_torch_NN
#' @description Allows user-friendly training set up of a torch model.
#' @param x Count matrix.
#' @param y Array containing the cell type class label.
#' @param batch_size Size of data loaded at training.
#' @param epoch Number of epochs.
#' @param prune_size
#' @param lr
#' @export
train_cell_scorer <- function(x, y, n_epoch=100L, batch_size=128L, prune_size=50L, lr=c(1e-3,1e-4,1e-5,1e-6)) {
  y <- as.factor(y)

  sch_epoch <- n_epoch %/% 2L
  prune_epoch <- sch_epoch - 1L

  # stop("lr not implemented yet: use approx to implement")

  fm <- nn_cell_scorer |>
    setup(optimizer=optim_adam,
          metrics = list(luz_metric_accuracy()),
          loss = nn_multi_margin_loss(margin = 1)
    ) |>
    set_hparams(feature_names=colnames(x), class_names=levels(y)) |>
    set_opt_hparams(weight_decay=0, lr=0.001) |>
    fit(
      epoch = n_epoch,
      data = delayed_xy_dataset(x,as.integer(y)),
      valid_data = delayed_xy_dataset(x,as.integer(y)),
      dataloader_options = list(batch_size=batch_size,drop_last=FALSE),
      # accelerator = accelerator(cpu=TRUE),
      callbacks = list(
        luz_callback_lr_scheduler(lr_lambda, verbose=TRUE, function(epoch) {0.1^(((epoch * 4L) %/% sch_epoch) %% 4L)}),
        luz_callback_prune("layer1.linear.weight", n=prune_size, at_epoch=prune_epoch)
      )
    )
  fm
}

