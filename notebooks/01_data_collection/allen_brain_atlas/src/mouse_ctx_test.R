
library(scater)
source("../lib_sce.R")


e10x <- loadHDF5SingleCellExperiment("out/mouse_ctx_10x.h5")
essv4 <- loadHDF5SingleCellExperiment("out/mouse_ctx_ssv4.h5")

ggcells(e10x[,sample(ncol(e10x),50000)]) + 
  geom_point(aes(x=UMAP.1,y=UMAP.2,color=subclass_label),size=1) + 
  coord_equal() + 
  theme_classic() +
  theme(legend.position="none")


