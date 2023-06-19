

install.packages(c('Matrix','curl','XML','rjson','jsonlite','rhdf5','igraph','reshape2','Rcpp','devtools'))
install.packages(c('LowRankQP','bmrm','limma'))
install.packages(c('scales','RColorBrewer','wesanderson'))
install.packages(c('ggplot2','patchwork','ggrepel','ggraph','cowplot','ggnewscale'))
install.packages(c('jpeg','png','ggpubr'))
install.packages(c('irlba', 'umap', 'Rtsne'))
install.packages(c('tidyverse','tidygraph'))
install.packages(c('reticulate','Seurat'))
install.packages(c('bs4Dash', 'reactlog', 'rentrez', 'httr2'))


install.packages(c('torch','luz'))
torch::install_torch()

devtools::install_github('bimberlab/cellhashR', ref = 'master', dependencies = TRUE)
devtools::install_github('cellgeni/sceasy')
devtools::install_github("AllenInstitute/cocoframer")


# bioconductor packages
install.packages("BiocManager", repos="https://cran.rstudio.com")
BiocManager::install(c('ShortReads','Rsamtools','rtracklayer'))
BiocManager::install(c('S4Vectors','IRanges','GenomicRanges','GenomicFeatures','GenomicAlignments','SummarizedExperiment'))
BiocManager::install(c('Rhdf5lib','rhdf5','HDF5Array'))
BiocManager::install(c('SingleCellExperiment','SpatialExperiment','LoomExperiment'))
BiocManager::install(c('DropletUtils','scuttle','scran','scater','uwot','AUCell','SingleR'))
BiocManager::install(c('biomaRt'))
BiocManager::install(c('edgeR','GEOquery','zellkonverter','sceasy'))
BiocManager::install(c('splatter'))
BiocManager::install(c('demuxmix'))









