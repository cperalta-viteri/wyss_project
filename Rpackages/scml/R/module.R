

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




