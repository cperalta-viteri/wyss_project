# source("R/delayed.R")

#' @title nn_linear_rescale
#' @description  A linear layer returning a re-scaled output.
#' @param in_features Input number of features.
#' @param out_features Output number of features.
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
#' @description Cell type classifier NN-based.
#' @param feature_names Names of the features.
#' @param class_names Class of the cells.
#' @param input_dropout_rate Dropout rate to apply on input values when training
#' @param n A numeric integer vector of length>0 specifying the number of feature in each internal layer.
#' @param dropout_rates A numeric vector of length>0 of dropouts to apply before each internal layer.
#' @import torch
#' @export
nn_cell_scorer <- torch::nn_module(
  inherit = nn_sequential,
  initialize = function(feature_names,class_names,input_dropout_rate=0.5,n=256L,dropout_rates=0.25) {
    super$initialize()

    # Process parameters
    self$feature_names <- as.character(feature_names)
    self$class_names <- as.character(class_names)
    n <- as.integer(n)
    stopifnot("n must be a non-empty integer vector"=length(n)>0L)
    stopifnot("dropout_rates must be a non-empty numeric vector"=length(dropout_rates)>0L)
    dropout_rates <- as.numeric(dropout_rates)
    dropout_rates <- rep_len(dropout_rates,length(n)) # Recycle to length(n)

    # First layer
    self$add_module("input_dropout",nn_dropout(input_dropout_rate))
    self$add_module("input_linear_rescaling",nn_linear_rescale(in_features = length(self$feature_names),out_features = n[1L]))
    self$add_module("input_tanh_indicators",nn_tanh())

    # Internal layers
    for(i in seq_along(n)[-1L]) {
      self$add_module(paste0(i,":dropout"),nn_dropout(dropout_rates[i-1L]))
      self$add_module(paste0(i,":linear"),nn_linear(in_features = n[i-1L], out_features = n[i]))
      self$add_module(paste0(i,":relu"),nn_relu())
    }

    # Last layer
    self$add_module("final_dropout",nn_dropout(dropout_rates[length(dropout_rates)]))
    self$add_module("final_linear",nn_linear(in_features = n[length(n)], out_features = length(self$class_names)))
  }
)





