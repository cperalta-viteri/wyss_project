library(HDF5Array)
library(SingleCellExperiment)
library(DelayedArray)
#library(rjson)
#library(BiocParallel)
#library(scuttle)
#library(Matrix)



#' Read json file storing the allen taxonomy and convert into a data.frame
#' @param path to the json file containing the taxonomy in allen20 format
#' @return a data.frame
read_dend_json <- function(json_file) {
  x <- rjson::fromJSON(file=json_file,simplify = FALSE)

  # recursive function converting into dcf format
  to_dcf <- function(n,path="") {
    # Check if the node is either a leaf or an internal node
    stopifnot(xor(setequal(names(n),"leaf_attributes"),setequal(names(n),c("children","node_attributes"))))

    # retrieve attributes
    att <- n[[grep("_attributes",names(n))]][[1]]
    stopifnot(hasName(att,"cell_set_accession"))
    att$path <- paste0(path,"/",att$cell_set_accession) # add path

    l <- mapply(to_dcf,n$children,att$path)
    c(paste(names(att),att,sep=":",collapse="\n"),unlist(l))
  }

  dcf <- to_dcf(x)
  dcf <- read.dcf(textConnection(paste(dcf,collapse="\n\n")))
  dcf <- as.data.frame(dcf)

  # write and reparse the data.frame to convert numeric values
  local({
    con <- textConnection(NULL,"w")
    write.table(dcf,file=con,sep="\t",quote=FALSE,na = "")
    read.table(text=textConnectionValue(con),sep="\t",quote="",comment="",check.names = FALSE,stringsAsFactors = FALSE)
  })
}


.download_url <- function(url,dir) {
  local_path <- file.path(dir,basename(url))
  download.file(url, local_path)
  local_path
}


.add_qc_metrics <- function(e) {
  BPPARAM <- BiocParallel::MulticoreParam(6,progressbar = TRUE)
  e <- scuttle::computeLibraryFactors(e,BPPARAM=BPPARAM)
  e <- scuttle::addPerCellQCMetrics(e,BPPARAM=BPPARAM)
  e <- scuttle::addPerFeatureQCMetrics(e,BPPARAM=BPPARAM)
  e
}

#' Save a SingleCellExperiment object into a folder with generic formats (H5 and TSV only)
#' Note: only counts() assay is saved, raise a warning if more assays available
.save_SCE <- function(e,dir,prefix) {
  saveRDS(rowData(e), file.path(dir,paste0(prefix,"_rowData.rds")))
  saveRDS(colData(e), file.path(dir,paste0(prefix,"_colData.rds")))
  saveRDS(reducedDims(e), file.path(dir,paste0(prefix,"_reducedDims.rds")))
  options(DelayedArray.block.size=1e9) # 1GB block size.
  writeHDF5Array(
    counts(e),
    file = file.path(dir,paste0(prefix,"_counts.h5")),
    name = "counts",
    chunkdim = NULL
  )
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Allen Brain SSv4 Mouse Atlas
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Download the files from GEO as well as the cell-type hierarchy
.download_ABAMouseSSv4 <- function(dir) {
  #dir <- "~/Downloads"
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/dend.json',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_expression_matrix_SSv4.hdf5',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_metadata_ssv4.csv.gz',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_umap2d_ssv4.csv.gz',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_CTX-HPF-73538sample_metadata.xlsx',dir)
}

.load_ABAMouseSSv4 <- function(dir) {
  #dir <- "../../notebooks/01_data_collection/allen_brain_atlas/data/GSE185862"
  files <- list(
    metadata = file.path(dir,"GSE185862_metadata_ssv4.csv.gz"),
    umap = file.path(dir,"GSE185862_umap2d_ssv4.csv.gz"),
    h5 = file.path(dir,"GSE185862_expression_matrix_SSv4.hdf5")
  )

  # Process cell-type hierarchy
  dend <- read_dend_json(file.path(dir,"dend.json"))

  # Load count matrix
  M <- Matrix::sparseMatrix(
    i = as.integer(h5read(files$h5,"/data/t_exon/i")),
    p = as.integer(h5read(files$h5,"/data/t_exon/p")),
    x = as.vector(h5read(files$h5,"/data/t_exon/x")),
    dims = as.vector(h5read(files$h5,"/data/t_exon/dims")),
    index1 = FALSE
  )
  colnames(M) <- as.character(h5read(files$h5,"/sample_names"))
  rownames(M) <- as.character(h5read(files$h5,"/gene_names"))

  # Process cells meta-data
  meta <- read.csv(files$metadata,check.names = FALSE,stringsAsFactors = FALSE,row.names=1L)
  stopifnot(all(meta$cell_type_accession_label %in% c("",dend$cell_set_accession)))
  rownames(meta) <- meta$sample_name
  stopifnot(all(rownames(meta) %in% colnames(M)))
  meta <- meta[colnames(M),]
  rownames(meta) <- colnames(M)

  # Load umap coordinates
  umap <- read.csv(files$umap,check.names = FALSE,stringsAsFactors = FALSE,row.names = 1) |>
    data.matrix()
  umap <- umap[match(colnames(M),rownames(umap)),]
  rownames(umap) <- colnames(M)

  e <- SingleCellExperiment(
    assays = list(counts = M),
    rowData = DataFrame(symbol=rownames(M)),
    colData = meta,
    metadata = list(taxonomy=dend),
    reducedDims = list(UMAP=umap)
  )

  e
}


#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Allen Brain 10x Mouse Atlas
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Download the files from GEO as well as the cell-type hierarchy
.download_ABAMouseTENx <- function(dir) {
  #dir <- "~/Downloads"
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/dend.json',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_expression_matrix_10x.hdf5',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_metadata_10x.csv.gz',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_umap.2d.10x.csv.gz',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_CTX-HPF-73538sample_metadata.xlsx',dir)
}

.load_ABAMouseTENx <- function(dir) {
  #dir <- "../../notebooks/01_data_collection/allen_brain_atlas/data/GSE185862/"
  files <- list(
    metadata = file.path(dir,"GSE185862_metadata_10x.csv.gz"),
    umap = file.path(dir,"GSE185862_umap.2d.10x.csv.gz"),
    h5 = file.path(dir,"GSE185862_expression_matrix_10x.hdf5")
  )

  # Process cell-type hierarchy
  dend <- read_dend_json(file.path(dir,"dend.json"))

  # Load count matrix
  M <- HDF5Array(files$h5,"data/counts")
  colnames(M) <- as.character(h5read(files$h5,"data/gene"))
  rownames(M) <- as.character(h5read(files$h5,"data/samples"))
  M <- t(M)

  # Process cells meta-data
  meta <- read.csv(files$metadata,check.names = FALSE,stringsAsFactors = FALSE,row.names=1L)
  stopifnot(all(meta$cell_type_accession_label %in% c("",dend$cell_set_accession)))
  rownames(meta) <- meta$sample_name
  stopifnot(all(rownames(meta) %in% colnames(M)))
  meta <- meta[colnames(M),]
  rownames(meta) <- colnames(M)

  # Load umap coordinates
  umap <- read.csv(files$umap,check.names = FALSE,stringsAsFactors = FALSE,row.names = 1) |>
    data.matrix()
  umap <- umap[match(colnames(M),rownames(umap)),]
  rownames(umap) <- colnames(M)

  e <- SingleCellExperiment(
    assays = list(counts = M),
    rowData = DataFrame(symbol=rownames(M)),
    colData = meta,
    metadata = list(taxonomy=dend),
    reducedDims = list(UMAP=umap)
  )

  e
}



#' @export
prepare_ABAMouseTENx <- function(outdir=".") {
  tmpdir <- tempdir()
  .download_ABAMouseTENx(tmpdir)
  e <- .add_qc_metrics(.load_ABAMouseTENx(tmpdir))
  .save_SCE(e,outdir,"ABAMouseTENx")
}

#' @export
ABAMouseTENx <- function(dir) {
}


#' @export
prepare_ABAMouseSSv4 <- function(outdir=".") {
  tmpdir <- tempdir()
  .download_ABAMouseSSv4(tmpdir)
  e <- .add_qc_metrics(.load_ABAMouseSSv4(tmpdir))
  .save_SCE(e,outdir,"ABAMouseSSv4")
}

#' @export
ABAMouseSSv4 <- function(dir) {
}






