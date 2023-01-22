
library(rhdf5)
library(SingleCellExperiment)
library(HDF5Array)


saveHDF5SingleCellExperiment <- function(x,dir,prefix="main_") {
  
  .h5_save_array_list <- function(l,h5_file,group,as.sparse=NA) {
    h5createGroup(h5_file,group) |> stopifnot()
    if (length(l)>0L && is.null(names(l))) stop("array list must be named")
    h5write(as.character(names(l)),h5_file,paste0(group,"/names"))
    for(i in seq_along(l)) {
      elt_path <- paste0(group,"/element_",i)
      if ((!is.na(as.sparse) & as.sparse) | (is.na(as.sparse) & is_sparse(l[[i]]))) {
        l[[i]] <- writeTENxMatrix(l[[i]],h5_file,elt_path)
      } else {
        l[[i]] <- writeHDF5Array(l[[i]],h5_file,elt_path)
      }
    }
    l
  }
  
  .save_sce_list <- function(l,h5_file,group,prefix) {
    h5createGroup(h5_file,group) |> stopifnot()
    h5write(as.character(names(l)),h5_file,paste0(group,"/names"))
    for(i in seq_along(l)) {
      saveHDF5SingleCellExperiment(l[[i]],dir=dir,prefix=paste0(prefix,i,"_"))
    }
  }

  stopifnot(is(x,"SummarizedExperiment"))
  dir.create(dir,showWarnings = FALSE)
  
  h5_file <- file.path(dir,paste0(prefix,"arrays.h5"))
  h5createFile(h5_file) |> stopifnot()
  
  .h5_save_array_list(assays(x),h5_file,'/assays')
  assays(x) <- list()
  
  if (is(x,"SingleCellExperiment")) {
    .h5_save_array_list(reducedDims(x),h5_file,'/reducedDims',as.sparse = FALSE)
    reducedDims(x) <- list()
    .save_sce_list(altExps(x),h5_file,'/altExps','altExps')
    altExps(x) <- list()
  }
  
  rds_file <- file.path(dir,paste0(prefix,"sce.rds"))
  saveRDS(x,rds_file)
}



loadHDF5SingleCellExperiment <- function(dir,prefix="main_") {
  
  .h5_load_array_list <- function(h5_file,group) {
    elt_names <- as.character(h5read(h5_file,paste0(group,"/names")))
    l <- list()
    for(i in seq_along(elt_names)) {
      elt_path <- paste0(group,"/element_",i)
      l[[i]] <- tryCatch(TENxMatrix(h5_file,elt_path),error=function(e){HDF5Array(h5_file,elt_path)})
    }
    names(l) <- elt_names
    l
  }
  
  .load_sce_list <- function(h5_file,group,prefix) {
    elt_names <- as.character(h5read(h5_file,paste0(group,"/names")))
    l <- list()
    for(i in seq_along(elt_names)) {
      l[[i]] <- loadHDF5SingleCellExperiment(dir,paste0(prefix,i,"_"))
    }
    names(l) <- elt_names
    l
  }
  
  h5_file <- file.path(dir,paste0(prefix,"arrays.h5"))
  rds_file <- file.path(dir,paste0(prefix,"sce.rds"))
  x <- readRDS(rds_file)
  assays(x) <- .h5_load_array_list(h5_file,"/assays")
  if (is(x,"SingleCellExperiment")) {
    reducedDims(x) <- .h5_load_array_list(h5_file,"/reducedDims")
    altExps(x) <- .load_sce_list(h5_file,"/altExps","altExps")
  }
  
  x
}



E <- DropletUtils::read10xCounts("../../../10x/filtered_feature_bc_matrix.h5")
E <- E[,1:10]
E <- splitAltExps(E,mcols(E)$Type,ref = "Gene Expression")
saveHDF5SingleCellExperiment(E,"out/test2")
h5ls("out/test2/main_arrays.h5")
loadHDF5SingleCellExperiment("out/test2")



