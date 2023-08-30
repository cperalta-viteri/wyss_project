
library(scml)
library(torch)
library(luz)
library(SingleCellExperiment)
library(HDF5Array)
DelayedArray::setAutoBPPARAM(BiocParallel::MulticoreParam(workers = 6L,progressbar = TRUE))

dir.create("/scratch/ABAMouseSSv4")
file.copy("~/scratch/ABAMouseSSv4/sce","/scratch/ABAMouseSSv4",recursive = TRUE)
dir.create("/scratch/ABAMouseTENx")
file.copy("~/scratch/ABAMouseTENx/sce","/scratch/ABAMouseTENx",recursive = TRUE)


aba_ds <- function() {
  # Load 10x and smart-seq cells
  allenSSq <- loadHDF5SummarizedExperiment("/scratch/ABAMouseSSv4/sce");allenSSq$source <- "ssv4"
  allen10x <- loadHDF5SummarizedExperiment("/scratch/ABAMouseTENx/sce");allen10x$source <- "10x"

  # Subset to common genes only
  allen10x <- allen10x[intersect(rownames(allen10x),rownames(allenSSq))]
  allenSSq <- allenSSq[rownames(allen10x)]

  # Bind the datasets
  mcols(allen10x) <- mcols(allenSSq) <- NULL
  mcol <- intersect(names(colData(allen10x)),names(colData(allenSSq)))
  colData(allen10x) <- colData(allen10x)[mcol]
  colData(allenSSq) <- colData(allenSSq)[mcol]
  allen <- cbind(allen10x,allenSSq)
  allen$source <- factor(allen$source)

  # Split cells into 3 datasets (trn/val/tst)
  set.seed(12345L)
  allen$mlset <- sample(
    factor(c("trn","val","tst"),c("trn","val","tst")),
    size = ncol(allen),
    replace = TRUE,
    prob = c(0.7,0.15,0.15)
  )

  allen
}


# Load ABA dataset
allen <- aba_ds()
allen <- allen[,!is.na(allen$subclass_label)] # Keep only cells with an assigned sublcass
allen <- allen[,!is.na(allen$cluster_label)]  # Keep only cells with an assigned cluster

# Normalize counts
allen <- scuttle::logNormCounts(allen,center_size_factors=FALSE,size.factors=allen$sum/1e4)

# Compute target labels to predict
allen$y_true <- paste(allen$class_label,allen$neighborhood_label,allen$subclass_label,allen$cluster_label,sep=" --- ") |>
  factor()
stopifnot(all(!is.na(allen$y_true)))

# Balance training set
trn_ds <- local({
  trn_ds <- allen[,allen$mlset %in% "trn"]
  trn_ids <- splitAsList(colnames(trn_ds),paste(trn_ds$y_true,trn_ds$source)) |>
    lapply(sample,1000,replace=TRUE) |>
    unlist() |>
    sample()
  trn_ds[,trn_ids]
})
print(table(trn_ds$y_true,trn_ds$source))





#logcounts(trn_ds) <- as.matrix(logcounts(trn_ds))
model <- train_delayed_classifier(
  x = t(logcounts(trn_ds)),
  y = trn_ds$y_true,
  pre_pruning_epoch = 25,post_pruning_epoch = 25,
  #accelerator = accelerator(cpu=TRUE),
  valid_data = 0.1,
  #weight_decay = 1e-4,
  auto_resume_path = file.path(params$model_dir,"aba_mouse_cluster_classifier.tmp"),
  #initial_weights_path = "~/Downloads/w.mw"
)
luz_save(model,"~/scratch/aba_mouse_cluster_classifier.luz")






