

#' @title ABAMouse_get_ISH
#' @description Fetches mouse ISH data from Allen Brain Atlas and save it into a H5 file.
#' @param dir Directory where to save H5 file.
#' @import cocoframer
#' @export
ABAMouse_get_ISH <- function(dir="ABAMouseISH"){
  if(!require(cocoframer)) devtools::install_github("AllenInstitute/cocoframer")
  # require(cocoframer)
  if(!dir.exists(dir)) dir.create(dir)

  # List of available genes
  G <- cocoframer::get_exp_gene_relationships() |>
    subset(gene_symbol!="")

  V <- with(
    G,
    lapply(id, cocoframer::get_aba_ish_data) |>
      setNames(paste0(gene_symbol,"/",id))
  )
  V <- simplify2array(V)
  HDF5Array::writeHDF5Array(V,file.path(dir,"ABAMouseISH.h5","/ish"),with.dimnames = TRUE)

}
