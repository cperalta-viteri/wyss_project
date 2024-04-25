## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE------------------------------------------------------------
library(scml)
library(torch)
library(luz)

## ----message=FALSE------------------------------------------------------------
fm <- train_delayed_classifier(
  data.matrix(iris[1:4]),
  iris$Species,
  pre_pruning_epoch = 100L,post_pruning_epoch = 100L,
  batch_size = 150L,
  valid_data = 0.1,
  input_dropout_rate = 0.1,
  accelerator = accelerator(cpu=TRUE)
)

## -----------------------------------------------------------------------------
plot(fm)

## ----eval=FALSE, warning=FALSE------------------------------------------------
#  pred <- predict_delayed(data.matrix(iris[1:4]),fm,accelerator = accelerator(cpu=TRUE))
#  table(iris$Species,max.col(pred),useNA="i")

## ----message=FALSE------------------------------------------------------------
fm <- train_delayed_classifier(
  data.matrix(iris[1:4]),
  iris$Species,
  y_type = "ordinal",
  batch_size = 150L,
  valid_data = 0.1,
  input_dropout_rate = 0.25,
  pre_pruning_epoch = 100L,post_pruning_epoch = 100L,
  n = 50,
  weight_decay = 0.9,
  accelerator = accelerator(cpu=TRUE)
)
plot(fm)

## -----------------------------------------------------------------------------
pred <- predict_delayed(data.matrix(iris[1:4]),fm,accelerator = accelerator(cpu=TRUE))
plot(pred~iris$Species)

