## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = TRUE
)

## ----message=FALSE------------------------------------------------------------
library(scml)
library(torch)
library(luz)
library(ABAData)
library(SingleCellExperiment)
DelayedArray::setAutoBPPARAM(BiocParallel::MulticoreParam(workers = 6L,progressbar = TRUE))

## ----eval=FALSE---------------------------------------------------------------
#  make_celltype_ds <- function() {
#    # Load 10x and smart-seq cells
#    allenSSq <- ABAHumanSSv4("../../data/ABAMouseSSv4/");allenSSq$source <- "ssv4"
#    allen10x <- ABAMouseTENx("../../data/ABAMouseTENx/");allen10x$source <- "10x"
#  
#    # Subset to common genes only
#    allen10x <- allen10x[intersect(rownames(allen10x),rownames(allenSSq))]
#    allenSSq <- allenSSq[rownames(allen10x)]
#  
#    # Bind the datasets
#    mcols(allen10x) <- mcols(allenSSq) <- NULL
#    mcol <- intersect(names(colData(allen10x)),names(colData(allenSSq)))
#    colData(allen10x) <- colData(allen10x)[mcol]
#    colData(allenSSq) <- colData(allenSSq)[mcol]
#    allen <- cbind(allen10x,allenSSq)
#    allen$source <- factor(allen$source)
#  
#    # Split cells into 3 datasets (trn/val/tst)
#    set.seed(12345L)
#    allen$mlset <- sample(
#      factor(c("trn","val","tst"),c("trn","val","tst")),
#      size = ncol(allen),
#      replace = TRUE,
#      prob = c(0.7,0.15,0.15)
#    )
#  
#    # Save allen dataset into an H5 format optimized for learning
#    HDF5Array::saveHDF5SummarizedExperiment(
#      allen,"~/Downloads/ABAMouseLearning",chunkdim = c(nrow(allen),1L),replace = TRUE
#    )
#  }
#  make_celltype_ds()

## ----eval=FALSE---------------------------------------------------------------
#  # Load dataset
#  allen <- HDF5Array::loadHDF5SummarizedExperiment("~/Downloads/ABAMouseLearning")
#  allen <- allen[,!is.na(allen$subclass_label)]
#  
#  # Compute target labels to predict
#  allen$y_true <- paste(allen$class_label,allen$neighborhood_label,allen$subclass_label,sep=" --- ") |>
#    factor()
#  stopifnot(all(!is.na(allen$y_true)))
#  
#  # Normalize counts
#  allen <- scuttle::logNormCounts(allen,center_size_factors=FALSE,size.factors=allen$sum/1e4)
#  
#  # Balance training set
#  trn_ds <- local({
#    trn_ds <- allen[,allen$mlset %in% "trn"]
#    trn_ids <- splitAsList(colnames(trn_ds),paste(trn_ds$y_true,trn_ds$source)) |>
#      lapply(sample,1000,replace=TRUE) |>
#      unlist() |>
#      sample()
#    trn_ds[,trn_ids]
#  })
#  print(table(trn_ds$y_true,trn_ds$source))

## ----eval=FALSE---------------------------------------------------------------
#  model <- train_delayed_classifier(
#    x = t(logcounts(trn_ds)),
#    y = trn_ds$y_true,
#    accelerator = accelerator(cpu=TRUE),
#    valid_data = 0.1,
#    weight_decay = 1e-4
#  )
#  luz_save(model,"~/Downloads/aba_mouse_subclass_classifier.luz")

## -----------------------------------------------------------------------------
model <- luz_load("~/Downloads/aba_mouse_subclass_classifier.luz")
plot(model)

## ----eval=FALSE---------------------------------------------------------------
#  pred <- predict_delayed(t(logcounts(trn_ds)),model,accelerator = accelerator(cpu = TRUE))
#  mean(trn_ds$y_true==colnames(pred)[max.col(pred)])

## ----eval=FALSE---------------------------------------------------------------
#  table(trn_ds$y_true,colnames(pred)[max.col(pred)])

