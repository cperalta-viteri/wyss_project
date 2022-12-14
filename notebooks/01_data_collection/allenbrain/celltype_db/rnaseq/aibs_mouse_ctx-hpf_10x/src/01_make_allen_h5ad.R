
library(IRanges)
library(HDF5Array)
library(rhdf5)
library(SingleCellExperiment)
library(zellkonverter)


source("../lib_allen.R")

make_h5ad <- function(dir="data/", name="new_h5ad") {
  # Process dendrogram
  dend <- read_dend_json(file.path(dir,"dend.json.gz"))
  
  # Process cells meta-data
  meta <- read.csv(file.path(dir,"metadata.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names=NULL)
  meta <- meta[sapply(meta,function(x) length(unique(x)))>1]
  stopifnot(all(meta$cell_type_accession_label %in% dend$cell_set_accession))
  
  # Add tsne coordinates
  meta <- local({
    tsne <- read.csv(file.path(dir,"tsne.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="sample_name")
    #stopifnot(setequal(meta$sample_name,rownames(tsne)))
    
    # total cells = 1169320 | meta contains 1169213 cell, ~100 less
    ## this can be used to filter out:
    #if (setequal(rownames(tsne), meta$sample_name) != TRUE) {
      #tsne <- tsne[!rownames(tsne) %in% rownames(tsne[c(setdiff(rownames(tsne), meta$sample_name)), ]),]
    #}
    tsne$sample_name <- rownames(tsne)
    merge(meta,tsne, by="sample_name", all=TRUE)
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
  stopifnot(all(rownames(M) == rownames(medians)))
  stopifnot(all(rownames(M) == rownames(trimmed_means)))
  # we cannot compare indices with gene names
  #stopifnot(all(colnames(M) == rownames(meta)))
  
  sce <- SingleCellExperiment(assays = list(counts=M), rowData = list(medians=medians), colData = list(cell_data=meta), metadata = list(taxonomy=dend))
  
  
  writeH5AD(sce = sce, file = name)
}


test <- sce[1:1000, 1:1000]
h5ad <- readH5AD(paste(dir, "SCE", sep = ""))









