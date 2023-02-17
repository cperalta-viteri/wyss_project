#Sys.setenv(PYTORCH_ENABLE_MPS_FALLBACK=1)
#Sys.getenv("PYTORCH_ENABLE_MPS_FALLBACK")

library(torch)
library(luz)
source("notebooks/celltype_prediction/lib_cell_classifier.R")
source("notebooks/01_data_collection/lib_sce.R")




x <- loadHDF5SingleCellExperiment("notebooks/01_data_collection/allen_brain_atlas/out/mouse_ctx_ssv4.h5")
x <- logNormCounts(x)

X <- x[1:100,1:100]
ds <- list(x=as.matrix(t(logcounts(X))),y=factor(X$cluster_label))
saveRDS(ds,"notebooks/celltype_prediction/dataset.rds")
