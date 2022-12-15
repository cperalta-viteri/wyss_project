
library(IRanges)
library(HDF5Array)
library(SingleCellExperiment)
library(zellkonverter)
library(rhdf5)
source("../lib_allen.R")



make_h5ad <- function(dir="data/", name="allen_mouse_smartSeq.h5ad") {
  # Process dendrogram
  dend <- read_dend_json(file.path(dir,"dend.json"))
  
  # Process cells meta-data
  meta <- read.csv(file.path(dir,"metadata.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names=NULL)
  meta <- meta[sapply(meta,function(x) length(unique(x)))>1]
  stopifnot(all(meta$cell_type_accession_label %in% c("",dend$cell_set_accession)))
  # put it below, after tsne merging
  #rownames(meta) <- meta$sample_name
  
  # Add tsne coordinates
  meta <- local({
    tsne <- read.csv(file.path(dir,"tsne.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="sample_name")
    if (!all(rownames(tsne) %in% meta$sample_name)) {
      tsne$sample_name <- rownames(tsne)
      meta <- merge(meta,tsne, by="sample_name", all=TRUE)
    }
    #stopifnot(all(rownames(tsne) %in% meta$sample_name))
    #cbind(meta,tsne[match(meta$sample_name,rownames(tsne)),])
    rownames(meta) <- meta$sample_name
    meta
  })
  
  # Process medians and trimmed means
  medians <- data.matrix(read.csv(file.path(dir,"medians.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  trimmed_means <- data.matrix(read.csv(file.path(dir,"trimmed_means.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  stopifnot(identical(dimnames(medians),dimnames(trimmed_means)))
  
  # empty create hdf5 file
  if (file.exists("data/matrix2.hdf5")){file.remove("data/matrix2.hdf5")}
  h5createFile("data/matrix2.hdf5")
  # creating groups inside hdf5 file
  h5createGroup("data/matrix2.hdf5","/matrix")
  h5createGroup("data/matrix2.hdf5","/matrix/features")
  # writing information in the groups
  ## exon/i == /matrix/indices
  writeHDF5Array(HDF5Array("data/matrix.hdf5","/data/t_exon/i"), "data/matrix2.hdf5","/matrix/indices")
  ## exon/p == /matrix/indptr
  writeHDF5Array(HDF5Array("data/matrix.hdf5","/data/t_exon/p"),"data/matrix2.hdf5","/matrix/indptr")
  ## exon/x == matrix/data
  writeHDF5Array(HDF5Array("data/matrix.hdf5","/data/t_exon/x"),"data/matrix2.hdf5","/matrix/data")
  ## exon/dims == matrix/shape
  writeHDF5Array(HDF5Array("data/matrix.hdf5","/data/t_exon/dims"),"data/matrix2.hdf5","/matrix/shape")
  # /sample_names is /matrix/barcode
  writeHDF5Array(HDF5Array("data/matrix.hdf5","/sample_names"),"data/matrix2.hdf5","/matrix/barcodes")
  # /gene_names is /matrix/features/name
  writeHDF5Array(HDF5Array("data/matrix.hdf5","/gene_names"),"data/matrix2.hdf5","/matrix/features/id")
  
  M <- TENxMatrix(file.path(dir,"matrix2.hdf5"))
  #M <- t(M)
  meta <- meta[order(order(match(colnames(M), rownames(meta)))),]
  
  file.remove("data/matrix2.hdf5")
  
  stopifnot(all(rownames(M) == rownames(medians)))
  stopifnot(all(rownames(M) == rownames(trimmed_means)))
  stopifnot(all(colnames(M) == rownames(meta)))
  
  sce <- SingleCellExperiment(assays = list(counts=M),
                              rowData = data.frame(medians=medians, trim_means=trimmed_means,check.names = FALSE),
                              colData = meta,
                              metadata = list(taxonomy=dend))
  
  
  writeH5AD(sce, file.path(dir, name))
}

make_h5ad()

