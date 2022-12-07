
library(IRanges)
library(HDF5Array)
source("../lib_allen.R")



make_hdf5 <- function(dir="data") {
  # Process dendrogram
  dend <- read_dend_json(file.path(dir,"dend.json"))
  
  # Process cells meta-data
  meta <- read.csv(file.path(dir,"metadata.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names=NULL)
  meta <- meta[sapply(meta,function(x) length(unique(x)))>1]
  stopifnot(all(meta$cell_type_accession_label %in% c("",dend$cell_set_accession)))
  rownames(meta) <- meta$sample_name
  
  # Add tsne coordinates
  meta <- local({
    tsne <- read.csv(file.path(dir,"tsne.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="sample_name")
    stopifnot(all(rownames(tsne) %in% meta$sample_name))
    cbind(meta,tsne[match(meta$sample_name,rownames(tsne)),])
  })
  
  # Process medians and trimmed means
  medians <- data.matrix(read.csv(file.path(dir,"medians.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  trimmed_means <- data.matrix(read.csv(file.path(dir,"trimmed_means.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  stopifnot(identical(dimnames(medians),dimnames(trimmed_means)))
  
  M <- read_sparseMatrix_from_csv(file.path(dir,"matrix.csv.gz"))
  M <- t(M)
  stopifnot(all(rownames(M) == rownames(medians)))
  stopifnot(all(rownames(M) == rownames(trimmed_means)))
  stopifnot(all(colnames(M) == rownames(meta)))
  
  
  
  #save(dend,medians,trimmed_means,meta,M,file = "out/objs.RData")
  #load("out/objs.RData")
  h5_file <- "out/mat.h5"
  h5write(dend,h5_file,"/dendrogram")
  writeHDF5Array(medians,h5_file,"/medians",with.dimnames = TRUE)
  writeHDF5Array(trimmed_means,h5_file,"/trimmed_means",with.dimnames = TRUE)
  h5write(meta,h5_file,"/cell_meta")
  writeTENxMatrix(M,h5_file,"/matrix")
}



