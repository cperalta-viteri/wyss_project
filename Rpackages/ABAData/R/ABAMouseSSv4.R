

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Allen Brain SSv4 Mouse Atlas
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

# Download the files from GEO as well as the cell-type hierarchy
.download_ABAMouseSSv4 <- function(dir) {
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_mouse_ctx-hpf_10x/dend.json',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_metadata_ssv4.csv.gz',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_umap2d_ssv4.csv.gz',dir)
  .download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_expression_matrix_SSv4.hdf5',dir)
  #.download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_sample_group_mapping.csv.gz',dir)
  #.download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_CTX-HPF-73538sample_metadata.xlsx',dir)
  #.download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_CTX_HPF_2021_SS_readme.txt',dir)
}

#' @import Matrix
#' @import SingleCellExperiment
.load_ABAMouseSSv4 <- function(dir) {
  #dir <- "ABAMouseSSv4/download/"
  files <- list(
    metadata = file.path(dir,"GSE185862_metadata_ssv4.csv.gz"),
    umap = file.path(dir,"GSE185862_umap2d_ssv4.csv.gz"),
    h5 = file.path(dir,"GSE185862_expression_matrix_SSv4.hdf5")
  )

  # Process cell-type hierarchy
  dend <- aba_read_dend_json(file.path(dir,"dend.json"))

  # Load count matrix
  M <- Matrix::sparseMatrix(
    i = as.integer(rhdf5::h5read(files$h5,"/data/t_exon/i")),
    p = as.integer(rhdf5::h5read(files$h5,"/data/t_exon/p")),
    x = as.vector(rhdf5::h5read(files$h5,"/data/t_exon/x")),
    dims = as.vector(rhdf5::h5read(files$h5,"/data/t_exon/dims")),
    index1 = FALSE
  )
  colnames(M) <- as.character(rhdf5::h5read(files$h5,"/sample_names"))
  rownames(M) <- as.character(rhdf5::h5read(files$h5,"/gene_names"))

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
  umap <- umap[match(colnames(M),sub("^Smartseq_cells.","",rownames(umap))),]
  rownames(umap) <- colnames(M)
  colnames(umap) <- NULL

  e <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = M),
    rowData = S4Vectors::DataFrame(symbol=rownames(M)),
    colData = meta,
    metadata = list(taxonomy=dend),
    reducedDims = list(UMAP=umap)
  )

  e
}

#' Allen Brain SmartSeq V4 Mouse Atlas
#'
#' Return a SingleCellExperiment object with 73k cells annotated
#' @param outdir name of the directory to store the dataset
#' @export
ABAMouseSSv4 <- function(outdir="ABAMouseSSv4") {
  SingleCellExperimentCachedDataset(outdir,.download_ABAMouseSSv4,.load_ABAMouseSSv4)
}




