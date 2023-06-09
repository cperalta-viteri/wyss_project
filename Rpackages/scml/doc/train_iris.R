## ---- include = FALSE---------------------------------------------------------
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
  batch_size = 150L,
  valid_data = 0.1,
  epoch = 200L,
  prune_at_epoch = NA,
  accelerator = accelerator(cpu=TRUE)
)

## -----------------------------------------------------------------------------
plot(fm)

## ----eval=FALSE, warning=FALSE------------------------------------------------
#  pred <- predict_delayed(data.matrix(iris[1:4]),fm)
#  table(iris$Species,max.col(pred))

