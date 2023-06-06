


# Build expression array of available genes
#' @import cocoframer
.download_ABAMouseISH <- function(dir) {
  
  # List of available genes
  G <- cocoframer::get_exp_gene_relationships() |>
    subset(gene_symbol!="")
  
  # Retrieve ISH 3D array for all available genes
  V <- BiocParallel::bplapply(
    BPPARAM = BiocParallel::SerialParam(progressbar = TRUE),
    G$id,cocoframer::get_aba_ish_data,values="energy"
  ) |> setNames(paste0(G$gene_symbol,"/",G$id))
  
  # Simplify into a 4D dataset
  V <- simplify2array(V)
  
  HDF5Array::writeHDF5Array(V,file.path(dir,"mouse_ish.h5"),"/energy",with.dimnames = TRUE)
}




#' Allen Brain ISH Mouse Atlas
#'
#' @param outdir name of the directory to store the dataset
#' @export
ABAMouseISH <- function(outdir="ABAMouseISH") {
  
}







