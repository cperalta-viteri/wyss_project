
# Download the files from GEO as well as the cell-type hierarchy
.download_ABAHumanTENx <- function(dir) {
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/dend.json',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/tsne.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv',dir)
  # .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/medians.csv',dir)
  # .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/trimmed_means.csv',dir)
  # .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/64/6d/646d3592-aff6-4364-8c3f-9e64b902638a/human_dendrogram.rds',dir)
  # .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0c/0c/0c0c882d-1c31-40a9-8039-3bf2706a77cd/sample-exp_component_mapping_human_10x_apr2020.zip',dir)
  # .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/70/32/70326830-e306-4743-a02c-a8da5bf9eb56/readme-m1-10.txt',dir)
  # .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/38/05/3805d61e-5c4d-446d-9331-cb61eb181198/taxonomy.txt',dir)
}


#' @import SingleCellExperiment
#' @import S4Vectors
.load_ABAHumanTENx <- function(dir) {
  #dir <- "../../notebooks/01_data_collection/allen_brain_atlas/data/GSE185862/"
  files <- list(
    metadata = file.path(dir,"metadata.csv"),
    tsne = file.path(dir,"tsne.csv"),
    h5 = file.path(dir,"matrix.csv"),
    medians = file.path(dir,"medians.csv"),
    trimmed_means = file.path(dir,"trimmed_means.csv")
  )

  # Process cell-type hierarchy
  dend <- aba_read_dend_json(file.path(dir,"dend.json"))

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
  # medians <- data.matrix(read.csv(files$medians,check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  # trimmed_means <- data.matrix(read.csv(files$trimmed_means,check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  # stopifnot(identical(dimnames(medians),dimnames(trimmed_means)))

  # check data alignment
  # stopifnot(all(rownames(M) == rownames(medians)))
  # stopifnot(all(rownames(M) == rownames(trimmed_means)))
  stopifnot(all(colnames(M) == rownames(meta)))

  # Build the object
  e <- SingleCellExperiment::SingleCellExperiment(
    list(counts = M),
    colData = meta,
    rowData = S4Vectors::DataFrame(symbol=rownames(M)),
    reducedDims = list(UMAP=tsne),
    metadata = list(taxonomy=dend)
  )
  e
}

#' @title ABAHumanTENx
#' @description Imports atlas data into R.
#' @param outdir Path to directory from where to read the data files.
#' @export
ABAHumanTENx <- function(outdir="ABAHumanTENx") {
  SingleCellExperimentCachedDataset(outdir,.download_ABAHumanTENx, .load_ABAHumanTENx)
}
