library(HDF5Array)
library(SingleCellExperiment)
library(DelayedArray)


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

read_sparseMatrix_from_csv <- function(file,batch_size=1024L) {
  con <- gzfile(file,"r")
  
  # Create the initial Matrix as an empty row, and set colnames
  M <- read.csv(con,nrows=batch_size,header=TRUE,row.names=1L,check.names = FALSE) |> 
    data.matrix()  |> as("RleArray")
  M <- list(M)
  
  # iterate over rows of the CSV
  while(nrow(M[[length(M)]])==batch_size) {
    message(sprintf("%d rows processed",sum(sapply(M,nrow))))
    M[[length(M)+1L]] <- read.csv(con,nrows=batch_size,header=FALSE,row.names=1L,check.names = FALSE) |> 
      data.matrix() |> as("RleArray")
  }
  M <- do.call(rbind,M)
  M
}


.download_url <- function(url,dir) {
  local_path <- file.path(dir,basename(url))
  download.file(url, local_path)
  local_path
}


.add_qc_metrics <- function(e) {
  e <- scuttle::computeLibraryFactors(e,BPPARAM=BiocParallel::MulticoreParam(6,progressbar = TRUE))
  e <- scuttle::addPerCellQCMetrics(e,BPPARAM=BiocParallel::MulticoreParam(6,progressbar = TRUE))
  e <- scuttle::addPerFeatureQCMetrics(e,BPPARAM=BiocParallel::MulticoreParam(6,progressbar = TRUE))
  e
}


# Allen Brain Atlas Mouse Cortex
download_ABAHumanTENx <- function(dir) {
  # Download the files from GEO as well as the cell-type hierarchy
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/matrix.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/metadata.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/dend.json',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/tsne.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/medians.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_m1_10x/trimmed_means.csv',dir)
  .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/64/6d/646d3592-aff6-4364-8c3f-9e64b902638a/human_dendrogram.rds',dir)
  .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/0c/0c/0c0c882d-1c31-40a9-8039-3bf2706a77cd/sample-exp_component_mapping_human_10x_apr2020.zip',dir)
  .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/70/32/70326830-e306-4743-a02c-a8da5bf9eb56/readme-m1-10.txt',dir)
  .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/38/05/3805d61e-5c4d-446d-9331-cb61eb181198/taxonomy.txt',dir)
}

download_ABAHumanSSv4 <- function(dir) {
  # Download the files from GEO as well as the cell-type hierarchy
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/matrix.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/intron.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/exon.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/trimmed_means.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/medians.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/metadata.csv',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/dend.json',dir)
  .download_url('https://idk-etl-prod-download-bucket.s3.amazonaws.com/aibs_human_ctx_smart-seq/tsne.csv',dir)
  .download_url('https://transcriptomic-viewer-downloads.s3-us-west-2.amazonaws.com/human/taxonomy.zip',dir)
  .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/30/80/30803016-9515-45fd-a6af-8595902c176f/readme_human_ctx.txt',dir)
  .download_url('https://brainmapportal-live-4cc80a57cd6e400d854-f7fdcae.divio-media.net/filer_public/b1/9f/b19f5e10-c6b2-4cf4-80fc-182e69a4511f/sample-exp_component_mapping_human_smart-seq_oct_2019.zip',dir)
}

#.download_url('ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/GSE185862_CTX-HPF-73538sample_metadata.xlsx',dir)

load_ABAHumanTENx <- function(dir) {
  #dir <- "../../notebooks/01_data_collection/allen_brain_atlas/data/GSE185862/"
  files <- list(
    metadata = file.path(dir,"metadata.csv"),
    tsne = file.path(dir,"tsne.csv"),
    h5 = file.path(dir,"matrix.csv"),
    medians = file.path(dir,"medians.csv"),
    trimmed_means = file.path(dir,"trimmed_means.csv")
  )
  
  # Process cell-type hierarchy
  dend <- read_dend_json(file.path(dir,"dend.json"))
  
  # Load count matrix
  M <- read_sparseMatrix_from_csv(files$h5,"data/counts")
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
  trimmed_means <- data.matrix(read.csv(files$trimmes_means,check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  stopifnot(identical(dimnames(medians),dimnames(trimmed_means)))
  
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


load_ABAHumanSSv4 <- function(dir) {
  #dir <- "../../notebooks/01_data_collection/allen_brain_atlas/data/GSE185862/"
  files <- list(
    metadata = file.path(dir,"metadata.csv"),
    tsne = file.path(dir,"tsne.csv"),
    h5 = file.path(dir,"matrix.csv"),
    medians = file.path(dir,"medians.csv"),
    trimmed_means = file.path(dir,"trimmed_means.csv")
  )

  # Process cell-type hierarchy
  dend <- read_dend_json(file.path(dir,"dend.json"))

  # Load count matrix
  M <- read_sparseMatrix_from_csv(files$h5,"data/counts")
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
  trimmed_means <- data.matrix(read.csv(files$trimmes_means,check.names = FALSE,stringsAsFactors = FALSE,row.names="feature"))
  stopifnot(identical(dimnames(medians),dimnames(trimmed_means)))
  
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
  
  # Compute librarySizeFactors
  #S <- .qc_metrics(e)
  e <- scuttle::computeLibraryFactors(e,BPPARAM=BiocParallel::MulticoreParam(6,progressbar = TRUE))
  e <- scuttle::addPerCellQCMetrics(e,BPPARAM=BiocParallel::MulticoreParam(6,progressbar = TRUE))
  e <- scuttle::addPerFeatureQCMetrics(e,BPPARAM=BiocParallel::MulticoreParam(6,progressbar = TRUE))
  e
}

load_ABAMouseSSv4 <- function(dir) {
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

  # Compute librarySizeFactors
  #S <- .qc_metrics(e)
  e <- scuttle::computeLibraryFactors(e,BPPARAM=BiocParallel::MulticoreParam(6,progressbar = TRUE))
  e <- scuttle::addPerCellQCMetrics(e,BPPARAM=BiocParallel::MulticoreParam(6,progressbar = TRUE))
  e <- scuttle::addPerFeatureQCMetrics(e,BPPARAM=BiocParallel::MulticoreParam(6,progressbar = TRUE))
  e
}



dir <- "../../notebooks/01_data_collection/allen_brain_atlas/data/GSE185862"
#download_ABAMouseTENx(dir)
#download_ABAMouseSSv4(dir)
e1 <- .add_qc_metrics(load_ABAMouseTENx(dir))
e2 <- .add_qc_metrics(load_ABAMouseSSv4(dir))




# Loading the data into R.
saveRDS(rowData(e1), "ABAMouseTENx_rowData.rds")
saveRDS(colData(e1), "ABAMouseTENx_colData.rds")
saveRDS(reducedDims(e1), "ABAMouseTENx_reducedDims.rds")
options(DelayedArray.block.size=1e9) # 1GB block size.
writeHDF5Array(
  counts(e1),
  file = "ABAMouseTENx_counts.h5",
  name = "counts"
)

saveRDS(rowData(e2), "ABAMouseSSv4_rowData.rds")
saveRDS(colData(e2), "ABAMouseSSv4_colData.rds")
saveRDS(reducedDims(e2), "ABAMouseSSv4_reducedDims.rds")
options(DelayedArray.block.size=1e9) # 1GB block size.
writeHDF5Array(
  counts(e2),
  file = "ABAMouseSSv4_counts.h5",
  name = "counts"
)


