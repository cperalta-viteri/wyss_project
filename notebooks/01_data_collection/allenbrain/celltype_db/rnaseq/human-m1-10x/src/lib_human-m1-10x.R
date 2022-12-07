
library(SummarizedExperiment)
library(HDF5Array)


human_m1_10x <- function(h5_file="out/mat.h5") {
  M <- TENxMatrix(h5_file,"/matrix")
  rd <- DataFrame(symbol=rownames(M),row.names=rownames(M))
  rd$medians <- HDF5Array(h5_file,"/medians")
  rd$trimmed_means <- HDF5Array(h5_file,"/trimmed_means")
  SummarizedExperiment(
    M,
    rowData = rd,
    colData = h5read(h5_file,"/cell_meta"),
    metadata = List(dendrogram=h5read(h5_file,"/dendrogram"))
  )  
}


