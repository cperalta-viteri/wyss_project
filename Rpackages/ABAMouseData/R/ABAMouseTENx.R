

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
  #.download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_sample_group_mapping.csv.gz',dir)
  #.download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_CTX-HPF-73538sample_metadata.xlsx',dir)
  #.download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_CTX_HPF_2021_10x_readme.txt',dir)
}

#' @import HDF5Array
#' @import SingleCellExperiment
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
ABAMouseTENx <- function(outdir="ABAMouseTENx") {
  SingleCellExperimentCachedDataset(outdir,.download_ABAMouseTENx,.load_ABAMouseTENx)
}







