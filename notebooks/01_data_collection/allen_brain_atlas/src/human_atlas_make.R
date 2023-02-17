

source("src/lib_allen.R")
source("../lib_sce.R")


read_human_ctx_sseq <- function(dir="data/human_ctx_sseq") {
  # Process dendrogram
  dend <- read_dend_json(file.path(dir,"dend.json"))
  
  # Process cells meta-data
  meta <- read.csv(file.path(dir,"metadata.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names=NULL)
  meta <- meta[sapply(meta,function(x) length(unique(x)))>1]
  stopifnot(all(meta$cell_type_accession_label %in% c("",dend$cell_set_accession)))
  rownames(meta) <- meta$sample_name
  
  # Add tsne coordinates
  tsne <- read.csv(file.path(dir,"tsne.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="sample_name") |>
    data.matrix()
  tsne <- tsne[match(meta$sample_name,rownames(tsne)),]
  rownames(tsne) <- meta$sample_name

  # Process medians and trimmed means
  medians <- data.matrix(read.csv(file.path(dir,"medians.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  trimmed_means <- data.matrix(read.csv(file.path(dir,"trimmed_means.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  stopifnot(identical(dimnames(medians),dimnames(trimmed_means)))
  
  # Read counts
  M <- read_sparseMatrix_from_csv(file.path(dir,"matrix.csv.gz"))
  M <- t(M)


  # check data alignment
  stopifnot(all(rownames(M) == rownames(medians)))
  stopifnot(all(rownames(M) == rownames(trimmed_means)))
  stopifnot(all(colnames(M) == rownames(meta)))

  # Build the object  
  e <- SingleCellExperiment(
    list(counts = M),
    colData = meta,
    rowData = DataFrame(symbol=rownames(M),medians=I(medians),trimmed_means=I(trimmed_means)),
    reducedDims = list(TSNE=tsne),
    metadata = list(taxonomy=dend)
  )
  e
}



read_human_m1_10x <- function(dir="data/human_m1_10x") {
  # Process dendrogram
  dend <- read_dend_json(file.path(dir,"dend.json"))
  
  # Process cells meta-data
  meta <- read.csv(file.path(dir,"metadata.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names=NULL)
  meta <- meta[sapply(meta,function(x) length(unique(x)))>1]
  stopifnot(all(meta$cell_type_accession_label %in% dend$cell_set_accession))
  rownames(meta) <- meta$sample_name
  
  # Add tsne coordinates
  tsne <- read.csv(file.path(dir,"tsne.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="sample_name") |>
    data.matrix()
  tsne <- tsne[match(meta$sample_name,rownames(tsne)),]
  rownames(tsne) <- meta$sample_name
  
  # Process medians and trimmed means
  medians <- data.matrix(read.csv(file.path(dir,"medians.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  trimmed_means <- data.matrix(read.csv(file.path(dir,"trimmed_means.csv.gz"),check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  stopifnot(identical(dimnames(medians),dimnames(trimmed_means)))
  
  # Read counts
  M <- read_sparseMatrix_from_csv(file.path(dir,"matrix.csv.gz"))
  M <- t(M)
  
  
  # check data alignment
  stopifnot(all(rownames(M) == rownames(medians)))
  stopifnot(all(rownames(M) == rownames(trimmed_means)))
  stopifnot(all(colnames(M) == rownames(meta)))
  
  # Build the object  
  e <- SingleCellExperiment(
    list(counts = M),
    colData = meta,
    rowData = DataFrame(symbol=rownames(M),medians=I(medians),trimmed_means=I(trimmed_means)),
    reducedDims = list(TSNE=tsne),
    metadata = list(taxonomy=dend)
  )
  e
}


saveHDF5SingleCellExperiment(read_human_ctx_sseq(),"out/human_ctx_sseq.h5")
saveHDF5SingleCellExperiment(read_human_m1_10x(),"out/human_m1_10x.h5")


