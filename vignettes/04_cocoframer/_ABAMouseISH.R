


# Build expression array of available genes
#' @importFrom cocoframer get_exp_gene_relationships get_aba_ish_data get_ccf_grid_annotation get_mba_ontology flatten_mba_ontology generate_taxons
#' @importFrom HDF5Array writeHDF5Array
#' @importFrom BiocParallel bplapply SerialParam
.download_ABAMouseISH <- function(dir) {
  
  # List of available genes
  G <- cocoframer::get_exp_gene_relationships() |>
    subset(gene_symbol!="")
  
  # Retrieve ISH 3D array for all available genes
  V <- bplapply(
    BPPARAM = SerialParam(progressbar = TRUE),
    G$id,cocoframer::get_aba_ish_data,values="energy"
  ) |> setNames(paste0(G$gene_symbol,"/",G$id))
  
  # Simplify into a 4D array
  V <- simplify2array(V)
  
  # Save array into a H5
  writeHDF5Array(V,file.path(dir,"ABAMouseISH.h5"),"/energy",with.dimnames = TRUE)
  
  # Retrieve and save CCF grid annotation
  A <- cocoframer::get_ccf_grid_annotation()
  writeHDF5Array(A,file.path(dir,"ABAMouseISH.h5"),"/ccf_grid_annotation",with.dimnames = TRUE)
  
  # retrieve and save mouse brain structure ontology
  o <- cocoframer::get_mba_ontology() |>
    cocoframer::flatten_mba_ontology() |>
    cocoframer::generate_taxons()
  write.table(o,gzfile(file.path(dir,"ABAMouseISH_ontology.tsv.gz")),sep="\t",row.names = FALSE)
  closeAllConnections()
}




#' Allen Brain ISH Mouse Atlas
#'
#' @param outdir name of the directory to store the dataset
#' @importFrom HDF5Array HDF5Array
#' @export
ABAMouseISH <- function(outdir="ABAMouseISH") {
  if (!file.exists(file.path(outdir,"ABAMouseISH.h5"))) {
    .download_ABAMouseISH(outdir)
  }
  HDF5Array::HDF5Array(file.path(outdir,"ABAMouseISH.h5"),"/energy")
}







