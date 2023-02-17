


# Build expression array of available genes
library(cocoframer)

# List of available genes
G <- get_exp_gene_relationships() |>
  subset(gene_symbol!="")

V <- with(
  G,
  lapply(id,get_aba_ish_data) |> 
    setNames(paste0(gene_symbol,"/",id))
)
V <- simplify2array(V)
HDF5Array::writeHDF5Array(V,"../out/mouse_ish.h5","/ish",with.dimnames = TRUE)











