


#' Save a SingleCellExperiment object into a folder with generic formats (H5 and TSV only)
#'
#' Important: only counts() assay is saved, raise a warning if more assays available
#' @param e a SingleCellExperiment object to write
#' @param dir directory to write dataset in
#' @importFrom HDF5Array writeHDF5Array
#' @importFrom SummarizedExperiment assayNames 
#' @importFrom SingleCellExperiment reducedDimNames rowData colData reducedDim counts
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
    filepath = file.path(dir,"counts.h5"),
    name = "counts",
    chunkdim = NULL
  )
}


#' Read a SingleCellExperiment from a directory saved with writeSingleCellExperiment
#'
#' @param dir a length-one character containing the directory to load data from
#' @importFrom HDF5Array HDF5Array
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @import utils
#' @export
readSingleCellExperiment <- function(dir) {
  # rowData
  rd <- utils::read.table(sep="\t", file.path(dir,"rowData.tsv.gz"), check.names = FALSE, comment.char = "")

  # colData
  cd <- utils::read.table(sep="\t", file.path(dir,"colData.tsv.gz"), check.names = FALSE, comment.char = "")

  # reducedDims
  red <- list.files(dir,"^reducedDim_") |>
    sub(pattern="^reducedDim_([^\\.]*)\\.tsv.*$",replacement="\\1")
  red <- lapply(red,function(r) {
    z <- utils::read.table(file.path(dir,sprintf("reducedDim_%s.tsv.gz",r)),sep="\t") |>
      data.matrix()
    colnames(z) <- NULL
    z
  }) |> stats::setNames(red)

  # Counts matrix
  n <- HDF5Array(filepath = file.path(dir,"counts.h5"),name = "counts")

  # Return SCE object
  SingleCellExperiment(list(counts=n),rowData=rd,colData=cd,reducedDims = red)
}



#' Generic function to create a cached SingleCellExperiment dataset
#'
#' The function automaticaly handle download and caching of the data
#' @param outdir path to directory used for I/O
#' @param download_fun function used to download the data into a given directory
#' @param load_fun function used to load the downloaded data from the given directory into a SingleCellExperiment object
#' @importFrom HDF5Array saveHDF5SummarizedExperiment
#' @importFrom HDF5Array loadHDF5SummarizedExperiment
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
    saveHDF5SummarizedExperiment(e,file.path(outdir,"sce"),chunkdim=c(nrow(e),1L),replace=TRUE)
    writeLines("DONE",file.path(outdir,"DONE"))
  } else {
    writeLines("Cache already built")
  }

  # Load cached data
  writeLines("Read SingleCellExperiment from disk")
  loadHDF5SummarizedExperiment(file.path(outdir,"sce"))
}

