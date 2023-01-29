

source("src/lib_allen.R")
source("../lib_sce.R")


optimize_DF_size <- function(X) {
  stopifnot(is(X,"DataFrame"))
  X[] <- lapply(X,function(v) {
    if (is.vector(v)) {
      if (is.character(v)) {
        f <- factor(v)
        if (object.size(f) < object.size(v)) v <- f
      }
      f <- Rle(v)
      if (object.size(f) < object.size(v)) v <- f
    }
    v
  })
  X
}


read_mouse_ctx_10x <- function(dir="data/GSE185862") {
  files <- list(
    metadata = file.path(dir,"GSE185862_metadata_10x.csv.gz"),
    umap = file.path(dir,"GSE185862_umap.2d.10x.csv.gz"),
    h5 = file.path(dir,"GSE185862_expression_matrix_10x.hdf5")
  )
  
  # Process dendrogram
  dend <- read_dend_json(file.path(dir,"dend.json"))
  
  # Process cells meta-data
  meta <- read.csv(files$metadata,check.names = FALSE,stringsAsFactors = FALSE,row.names=1L)
  stopifnot(all(meta$cell_type_accession_label %in% c("",dend$cell_set_accession)))
  rownames(meta) <- meta$sample_name
  
  # Load umap coordinates
  umap <- read.csv(files$umap,check.names = FALSE,stringsAsFactors = FALSE,row.names = 1) |>
    data.matrix()
  umap <- umap[match(meta$sample_name,rownames(umap)),]
  rownames(umap) <- meta$sample_name
  
  M <- HDF5Array(files$h5,"data/counts")
  colnames(M) <- as.character(h5read(files$h5,"data/gene"))
  rownames(M) <- as.character(h5read(files$h5,"data/samples"))
  M <- t(M)
  
  stopifnot(all(rownames(meta) %in% colnames(M)))
  M <- M[,match(rownames(meta),colnames(M),nomatch = 0)]
  
  e <- SingleCellExperiment(
    assays = list(counts = M),
    rowData = DataFrame(symbol=rownames(M)),
    colData = meta,
    metadata = list(taxonomy=dend),
    reducedDims = list(UMAP=umap)
  )
  e
}

read_mouse_ctx_ssv4 <- function(dir="data/GSE185862") {
  files <-  list(
    metadata = file.path(dir,"GSE185862_metadata_ssv4.csv.gz"),
    umap = file.path(dir,"GSE185862_umap2d_ssv4.csv.gz"),
    h5 = file.path(dir,"GSE185862_expression_matrix_SSv4.hdf5")
  )
  
  # Process dendrogram
  dend <- read_dend_json(file.path(dir,"dend.json"))
  
  # Process cells meta-data
  meta <- read.csv(files$metadata,check.names = FALSE,stringsAsFactors = FALSE,row.names=1L)
  stopifnot(all(meta$cell_type_accession_label %in% c("",dend$cell_set_accession)))
  rownames(meta) <- meta$sample_name
  
  # Load umap coordinates
  umap <- read.csv(files$umap,check.names = FALSE,stringsAsFactors = FALSE,row.names = 1) |>
    data.matrix()
  rownames(umap) <- sub("^Smartseq_cells\\.","",rownames(umap))
  umap <- umap[match(meta$sample_name,rownames(umap)),]
  rownames(umap) <- meta$sample_name
  
  M <- Matrix::sparseMatrix(
    i = as.integer(h5read(files$h5,"/data/t_exon/i")),
    p = as.integer(h5read(files$h5,"/data/t_exon/p")),
    x = as.vector(h5read(files$h5,"/data/t_exon/x")),
    dims = as.vector(h5read(files$h5,"/data/t_exon/dims")),
    index1 = FALSE
  )
  colnames(M) <- as.character(h5read(files$h5,"/sample_names"))
  rownames(M) <- as.character(h5read(files$h5,"/gene_names"))
  
  stopifnot(all(rownames(meta) %in% colnames(M)))
  M <- M[,match(rownames(meta),colnames(M),nomatch = 0)]
  e <- SingleCellExperiment(
    assays = list(counts = M),
    rowData = DataFrame(symbol=rownames(M)),
    colData = meta,
    metadata = list(taxonomy=dend),
    reducedDims = list(UMAP=umap)
  )
  e
}


make_mouse_ctx_10x <- function() {
  e <- read_mouse_ctx_10x()
  saveHDF5SingleCellExperiment(e,"out/mouse_ctx_10x.h5",sparse.assays = character())
}

make_mouse_ctx_ssv4 <- function() {
  e <- read_mouse_ctx_ssv4()
  assay(e,withDimnames = FALSE) <- writeHDF5Array(assay(e),H5type = "H5T_STD_I32LE",as.sparse = FALSE)
  saveHDF5SingleCellExperiment(e,"out/mouse_ctx_ssv4.h5")
}




  
  
