
# Download the files from GEO as well as the cell-type hierarchy
.download_ABAHumanSSv4 <- function(dir) {
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/intron.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/exon.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/trimmed_means.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/medians.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/dend.json',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/tsne.csv',dir)
  #.download_url('https://transcriptomic-viewer-downloads.s3-us-west-2.amazonaws.com/human/taxonomy.zip',dir)
  #.download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/30/80/30803016-9515-45fd-a6af-8595902c176f/readme_human_ctx.txt',dir)
  #.download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b1/9f/b19f5e10-c6b2-4cf4-80fc-182e69a4511f/sample-exp_component_mapping_human_smart-seq_oct_2019.zip',dir)
}



#' @import Matrix
#' @import SingleCellExperiment
.load_ABAHumanSSv4 <- function(dir) {
  #dir <- "../../notebooks/01_data_collection/allen_brain_atlas/data/GSE185862/"
  files <- list(
    metadata = file.path(dir,"metadata.csv"),
    tsne = file.path(dir,"tsne.csv"),
    h5 = file.path(dir,"matrix.csv"),
    medians = file.path(dir,"medians.csv"),
    trimmed_means = file.path(dir,"trimmed_means.csv")
  )
  # browser()
  # Process cell-type hierarchy
  dend <- .read_dend_json(file.path(dir,"dend.json"))

  # Load count matrix
  M <- .read_sparseMatrix_from_csv(files$h5)
  M <- t(M)

  # Process cells meta-data
  meta <- read.csv(files$metadata,check.names = FALSE,stringsAsFactors = FALSE,row.names=NULL)
  meta <- meta[sapply(meta,function(x) length(unique(x)))>1]
  stopifnot(all(meta$cell_type_accession_label %in% c("",dend$cell_set_accession)))
  rownames(meta) <- meta$sample_name

  # Load tsne coordinates
  tsne <- read.csv(files$tsne,check.names = FALSE,stringsAsFactors = FALSE,row.names = "sample_name") |>
    data.matrix()
  tsne <- tsne[match(colnames(M),rownames(tsne)),]
  rownames(tsne) <- colnames(M)

  # Process medians and trimmed means
  medians <- data.matrix(read.csv(files$medians,check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  trimmed_means <- data.matrix(read.csv(files$trimmed_means,check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  stopifnot(identical(dimnames(medians),dimnames(trimmed_means)))

  # check data alignment
  stopifnot(all(rownames(M) == rownames(medians)))
  stopifnot(all(rownames(M) == rownames(trimmed_means)))
  stopifnot(all(colnames(M) == rownames(meta)))

  # Build the object
  e <- SingleCellExperiment::SingleCellExperiment(
    list(counts = M),
    colData = meta,
    rowData = S4Vectors::DataFrame(symbol=rownames(M),medians=I(medians),trimmed_means=I(trimmed_means)),
    reducedDims = list(TSNE=tsne),
    metadata = list(taxonomy=dend)
  )
  e

}


#' @title ABAHumanTENx
#' @description Imports atlas data into R.
#' @param outdir Path to directory from where to read the data files.
#' @export
ABAHumanSSv4 <- function(outdir="ABAHumanSS") {
  SingleCellExperimentCachedDataset(outdir,.download_ABAHumanSSv4, .load_ABAHumanSSv4)
}
