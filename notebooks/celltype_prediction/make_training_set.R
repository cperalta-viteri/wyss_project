#Sys.setenv(PYTORCH_ENABLE_MPS_FALLBACK=1)
#Sys.getenv("PYTORCH_ENABLE_MPS_FALLBACK")

#library(torch)
#library(luz)
source("lib_cell_classifier.R")
source("../01_data_collection/lib_sce.R")




x <- loadHDF5SingleCellExperiment("../01_data_collection/allen_brain_atlas/out/mouse_ctx_ssv4.h5")

x <- scuttle::logNormCounts(x)

# 10_Lamp5 | 100_Sst | 91_Sst | 70_Sst | 236_Car3 | 272_NP SUB
#  125    |   97     |  126   |  119   | 191      |  80

v <- c("10_Lamp5", "100_Sst", "91_Sst", "70_Sst", "236_Car3", "272_NP SUB")

#X <- x[1:100,1:100]
X <- x[, x$cluster_label %in% v]
ds <- list(x=as.matrix(t(logcounts(X))),y=factor(X$cluster_label))
saveRDS(ds,"dataset.rds")

