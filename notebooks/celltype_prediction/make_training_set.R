#Sys.setenv(PYTORCH_ENABLE_MPS_FALLBACK=1)
#Sys.getenv("PYTORCH_ENABLE_MPS_FALLBACK")

library(torch)
library(luz)
source("notebooks/celltype_prediction/lib_cell_classifier.R")
source("notebooks/01_data_collection/lib_sce.R")




x <- loadHDF5SingleCellExperiment("notebooks/01_data_collection/allen_brain_atlas/out/mouse_ctx_ssv4.h5")
x <- scuttle::logNormCounts(x)

x$y_true <- factor(paste(x$class_label,x$neighborhood_label,x$subclass_label,sep=" --- "))
x$loss_group <- "all"

x$y_true <- factor(x$region_label)
x$loss_group <- x$subclass_label

x$y_true <- factor(x$y_true,names(which(table(x$y_true)>=10)))
x <- x[,!is.na(x$y_true)]
set.seed(1234L)
trn_cells <- splitAsList(colnames(x),x$y_true) |>
  lapply(sample,size=100,replace=TRUE) |>
  unlist() |>
  sample()
x_trn <- x[,trn_cells]
x_trn <- x_trn[rowMeans(assay(x_trn)>0)>0.01]
ds <- list(
  x = as.matrix(t(logcounts(x_trn))),
  y = x_trn$y_true,
  loss_group = x_trn$loss_group
)
#saveRDS(ds,"notebooks/celltype_prediction/subclass_dataset.rds")
saveRDS(ds,"notebooks/celltype_prediction/region_dataset.rds")





