
library(IRanges)
library(HDF5Array)
library(rhdf5)
library(SingleCellExperiment)
library(zellkonverter)


source("../lib_allen.R")

#reticulate::use_virtualenv("/home/christian97/Documents/Bioinfo_platf/wyss_project/docker_envi/python_3_8_venv")

Sys.setenv(RETICULATE_PYTHON="/home/christian97/Documents/Bioinfo_platf/wyss_project/docker_envi/python_3_8_venv/bin/python")


make_h5ad <- function(dir="data/", name="allen_mouse_TENx.h5ad") {
  # Process dendrogram
  dend <- read_dend_json(file.path(dir,"dend.json.gz"))
  
  # Process cells meta-data
  meta <- read.csv(file.path(dir,"metadata.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names=NULL)
  meta <- meta[sapply(meta,function(x) length(unique(x)))>1]
  stopifnot(all(meta$cell_type_accession_label %in% dend$cell_set_accession))
  
  
  tsne <- read.csv(file.path(dir,"tsne.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="sample_name")
  # Checking that meta and tsne have the same number of cells
  meta <- local({
    #stopifnot(setequal(meta$sample_name,rownames(tsne)))
    # total cells = 1169320 | meta contains 1169213 cell, ~100 less
    if (!all(rownames(tsne) %in% meta$sample_name)) {
      ## this can be used to filter out:
      #tsne <- tsne[!rownames(tsne) %in% rownames(tsne[c(setdiff(rownames(tsne), meta$sample_name)), ]),]
      tsne$sample_name <- rownames(tsne)
      meta <- merge(meta,tsne, by="sample_name", all=TRUE)
    }
    rownames(meta) <- meta$sample_name
    meta
  })
  
  # Process medians and trimmed means
  medians <- data.matrix(read.csv(file.path(dir,"medians.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  trimmed_means <- data.matrix(read.csv(file.path(dir,"trimmed_means.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  stopifnot(identical(dimnames(medians),dimnames(trimmed_means)))
  
  M <- HDF5Array(file.path(dir,"matrix.hdf5"),"/data/counts")
  colnames(M) <- as.character(HDF5Array(file.path(dir,"matrix.hdf5"), name = "data/gene"))
  rownames(M) <- as.character(HDF5Array(file.path(dir, "matrix.hdf5"), name= "/data/samples")) # M does not have identical sample names as present in tsne data
  # run: identical(rownames(M), tsne$sample_name) # FALSE
  
  
  M <- t(M)
  meta <- meta[order(order(match(colnames(M), rownames(meta)))),]
  tsne <- tsne[order(order(match(colnames(M), rownames(tsne)))),]
  stopifnot(all(rownames(M) == rownames(medians)))
  stopifnot(all(rownames(M) == rownames(trimmed_means)))
  stopifnot(all(colnames(M) == rownames(meta)))
  
  #browser()
  sce <- SingleCellExperiment(assays = list(counts=M),
                              rowData = data.frame(medians=medians, trim_means=trimmed_means,check.names = FALSE),
                              colData = meta,
                              metadata = list(taxonomy=dend),
                              reducedDims = list(umap=tsne))
  
  
  writeH5AD(sce = sce, file = file.path(dir, name))
}


make_h5ad()




#test <- sce[1:1000, 1:1000]
#writeH5AD(test, "data/sce.h5ad")
#h5ad <- readH5AD("data/sce.h5ad", use_hdf5 = TRUE)



library(sceasy)
library(reticulate)
loompy <- reticulate::import('loompy')

























