


#' Read json file storing the allen taxonomy and convert into a data.frame
#' @param path to the json file containing the taxonomy in allen20 format
#' @return a data.frame
#' @importFrom rjson fromJSON
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

#' @import scuttle
#' @import BiocParallel
.add_qc_metrics <- function(e) {
  BPPARAM <- BiocParallel::MulticoreParam(6,progressbar = TRUE)
  e <- scuttle::computeLibraryFactors(e,BPPARAM=BPPARAM)
  e <- scuttle::addPerCellQCMetrics(e,BPPARAM=BPPARAM)
  e <- scuttle::addPerFeatureQCMetrics(e,BPPARAM=BPPARAM)
  e
}

#' Save a SingleCellExperiment object into a folder with generic formats (H5 and TSV only)
#' Note: only counts() assay is saved, raise a warning if more assays available
#' @import HDF5Array
#' @import SingleCellExperiment
.write_SCE <- function(e,dir) {
  stopifnot(!anyDuplicated(reducedDimNames(e)))
  stopifnot("counts" %in% assayNames(e))
  if (!identical(assayNames(e),"counts")) warning("only counts assay will be saved.")
  write.table(rowData(e), sep="\t", file.path(dir,"rowData.tsv"))
  write.table(colData(e), sep="\t", file.path(dir,"colData.tsv"))
  for(rdn in reducedDimNames(e)) {
    write.table(reducedDim(e,rdn), sep="\t", file.path(dir,sprintf("reducedDim_%s.tsv",rdn)))
  }
  writeHDF5Array(
    counts(e),
    file = file.path(dir,"counts.h5"),
    name = "counts",
    chunkdim = NULL
  )
}


#' Read a SingleCellExperiment from a directory savec with .save_SCE
#' Note: only counts() assay is saved, raise a warning if more assays available
#' @import HDF5Array
#' @import SingleCellExperiment
.read_SCE <- function(e,dir) {
  rd <- read.table(sep="\t", file.path(dir,"rowData.tsv"))
  cd <- read.table(sep="\t", file.path(dir,"colData.tsv"))
  # for(rdn in reducedDimNames(e)) {
  #   write.table(reducedDim(e,rdn), sep="\t", file.path(dir,sprintf("reducedDim_%s.tsv",rdn)))
  # }
  n <- HDF5Array(file = file.path(dir,"_counts.h5"),name = "counts")
  SingleCellExperiment(list(counts=n),rowData=rd,colData(cd))
}







