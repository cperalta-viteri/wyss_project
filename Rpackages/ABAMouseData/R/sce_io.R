


#' Save a SingleCellExperiment object into a folder with generic formats (H5 and TSV only)
#' Note: only counts() assay is saved, raise a warning if more assays available
#' @import HDF5Array
#' @import SingleCellExperiment
#' @export
writeSingleCellExperiment <- function(e,dir) {
  stopifnot(!anyDuplicated(reducedDimNames(e)))
  stopifnot("counts" %in% assayNames(e))
  if (!identical(assayNames(e),"counts")) warning("only counts assay will be saved.")
  write.table(rowData(e), sep="\t", gzfile(file.path(dir,"rowData.tsv.gz")))
  write.table(colData(e), sep="\t", gzfile(file.path(dir,"colData.tsv.gz")))
  for(rdn in reducedDimNames(e)) {
    write.table(reducedDim(e,rdn), sep="\t", gzfile(file.path(dir,sprintf("reducedDim_%s.tsv.gz",rdn))))
  }
  closeAllConnections()
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
#' @export
readSingleCellExperiment <- function(dir) {
  rd <- read.table(sep="\t", file.path(dir,"rowData.tsv.gz"))
  cd <- read.table(sep="\t", file.path(dir,"colData.tsv.gz"))
  red <- list.files(dir,"^reducedDim_") |>
    sub(pattern="^reducedDim_([^\\.]*)\\.tsv.*$",replacement="\\1")
  red <- lapply(red,function(r) {
    read.table(file.path(dir,sprintf("reducedDim_%s.tsv.gz",r)),sep="\t") |>
      data.matrix()
  }) |> setNames(red)
  n <- HDF5Array(file = file.path(dir,"counts.h5"),name = "counts")
  SingleCellExperiment(list(counts=n),rowData=rd,colData=cd,reducedDims = red)
}



#' Generic function to create a SingleCellExperiment dataset
#' The function handle automatic download and caching of the data
#' @import
#' @export
SingleCellExperimentCachedDataset <- function(outdir,download_fun,load_fun) {
  # Download data and prepare QC
  if (!file.exists(file.path(outdir,"DONE"))) {
    writeLines("Start cache building")
    tmpdir <- file.path(outdir,"download")
    if (!file.exists(file.path(tmpdir,"DONE"))) {
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

