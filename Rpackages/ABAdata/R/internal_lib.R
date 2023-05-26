

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

#' Read json file storing the allen taxonomy and convert into a data.frame
#' @param path to the json file containing the taxonomy in allen20 format
#' @return a data.frame
.read_dend_json <- function(json_file) {
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

.read_sparseMatrix_from_csv <- function(file,batch_size=1024L) {
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
  dend <- .read_dend_json(file.path(dir,"dend.json"))
  
  # Load count matrix
  M <- .read_sparseMatrix_from_csv(files$h5,"data/counts")
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


.load_ABAHumanSSv4 <- function(dir) {
  #dir <- "../../notebooks/01_data_collection/allen_brain_atlas/data/GSE185862/"
  files <- list(
    metadata = file.path(dir,"metadata.csv"),
    tsne = file.path(dir,"tsne.csv"),
    h5 = file.path(dir,"matrix.csv"),
    medians = file.path(dir,"medians.csv"),
    trimmed_means = file.path(dir,"trimmed_means.csv")
  )
  
  # Process cell-type hierarchy
  dend <- .read_dend_json(file.path(dir,"dend.json"))
  
  # Load count matrix
  M <- .read_sparseMatrix_from_csv(files$h5,"data/counts")
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
