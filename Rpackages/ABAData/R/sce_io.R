


#' Save a SingleCellExperiment object into a folder with generic formats (H5 and TSV only)
#'
#' Important: only counts() assay is saved, raise a warning if more assays available
#' @param e a SingleCellExperiment object to write
#' @param dir directory to write dataset in
#' @import HDF5Array
#' @import SingleCellExperiment
#' @export
writeSingleCellExperiment <- function(e,dir) {
  stopifnot(!anyDuplicated(reducedDimNames(e)))
  stopifnot("counts" %in% assayNames(e))
  if (!identical(assayNames(e),"counts")) warning("only counts assay will be saved.")
  dir.create(dir)

  # write rowData
  write.table(rowData(e), sep="\t", gzfile(file.path(dir,"rowData.tsv.gz")))

  # write colData
  write.table(colData(e), sep="\t", gzfile(file.path(dir,"colData.tsv.gz")))

  # write reducedDims
  for(rdn in reducedDimNames(e)) {
    write.table(reducedDim(e,rdn), sep="\t", gzfile(file.path(dir,sprintf("reducedDim_%s.tsv.gz",rdn))))
  }
  closeAllConnections()

  # write counts
  writeHDF5Array(
    counts(e),
    file = file.path(dir,"counts.h5"),
    name = "counts",
    chunkdim = NULL
  )
}


#' Read a SingleCellExperiment from a directory saved with .save_SCE
#'
#' @param dir a length-one character containing the directory to load data from
#' @import HDF5Array
#' @import SingleCellExperiment
#' @export
readSingleCellExperiment <- function(dir) {
  # rowData
  rd <- read.table(sep="\t", file.path(dir,"rowData.tsv.gz"), check.names = FALSE, comment.char = "")

  # colData
  cd <- read.table(sep="\t", file.path(dir,"colData.tsv.gz"), check.names = FALSE, comment.char = "")

  # reducedDims
  red <- list.files(dir,"^reducedDim_") |>
    sub(pattern="^reducedDim_([^\\.]*)\\.tsv.*$",replacement="\\1")
  red <- lapply(red,function(r) {
    z <- read.table(file.path(dir,sprintf("reducedDim_%s.tsv.gz",r)),sep="\t") |>
      data.matrix()
    colnames(z) <- NULL
    z
  }) |> setNames(red)

  # Counts matrix
  n <- HDF5Array(file = file.path(dir,"counts.h5"),name = "counts")

  # Return SCE object
  SingleCellExperiment(list(counts=n),rowData=rd,colData=cd,reducedDims = red)
}



#' Generic function to create a cached SingleCellExperiment dataset
#'
#' The function automaticaly handle download and caching of the data
#' @param download_fun function used to download the data into a given directory
#' @param load_fun function used to load the downloaded data from the given directory into a SingleCellExperiment object
#' @import SingleCellExperiment
#' @export
SingleCellExperimentCachedDataset <- function(outdir,download_fun,load_fun) {
  # Download data and prepare QC
  if (!file.exists(file.path(outdir,"DONE"))) {
    writeLines("Start cache building")
    tmpdir <- file.path(outdir,"download")
    if (!file.exists(file.path(tmpdir,"DONE"))) {
      dir.create(tmpdir,recursive=TRUE)
      writeLines("Start data downloading")
      download_fun(tmpdir)
      writeLines("DONE",file.path(tmpdir,"DONE"))
    } else {
      writeLines("Data already downloaded")
    }
    writeLines("Load downloaded data")
    e <- load_fun(tmpdir)
    writeLines("Compute QC")
    e <- .add_qc_metrics(e)
    writeLines("Write SingleCellExperiment on disk")
    writeSingleCellExperiment(e,outdir)
    writeLines("DONE",file.path(outdir,"DONE"))
  } else {
    writeLines("Cache already built")
  }

  # Load cached data
  writeLines("Read SingleCellExperiment from disk")
  readSingleCellExperiment(outdir)
}

