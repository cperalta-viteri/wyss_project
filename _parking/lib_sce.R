
library(rhdf5)
library(SingleCellExperiment)
library(HDF5Array)


#' Save a SingleCellExperiment object within a H5
#' assays and reducedDims are saved as dense or sparseMatrix
#' altExps are saved with recursive calls
#' @param x a SingleCellExperiment (or SummarizedExperiment) object to that into a h5 file
#' @param h5_file name of the h5 file to create or update
#' @param group name of the group to create within the h5
#' @param force a logical if FALSE an error is raised whrn the h5 file already exists
#' @param sparse.assays name of the assays to store as sparse in the H5. If TRUE, all assays are stored as sparse; 
#'        if FALSE all assay are stored as dense; if NA sparsity is determined with a call to is_sparse()
#' @return undefined
saveHDF5SingleCellExperiment <- function(x,h5_file,group="/main",force=FALSE,sparse.assays=c("counts","logcounts")) {
  
  # Helper which can save a list of array in H5
  .h5_save_array_list <- function(l,h5_file,group,as.sparse=NA) {
    h5createGroup(h5_file,group) |> stopifnot()
    if (length(l)>0L && is.null(names(l))) stop("array list must be named")
    h5write(as.character(names(l)),h5_file,paste0(group,"/names"))
    
    if (is.character(as.sparse)) as.sparse <- names(l) %in% as.sparse
    as.sparse <- rep_len(as.sparse,length(l))
    as.sparse[is.na(as.sparse)] <- as.logical(sapply(l[is.na(as.sparse)],is_sparse))
    for(i in seq_along(l)) {
      elt_path <- paste0(group,"/element_",i)
      if (as.sparse[i]) {
        l[[i]] <- writeTENxMatrix(l[[i]],h5_file,elt_path)
      } else {
        l[[i]] <- writeHDF5Array(l[[i]],h5_file,elt_path,with.dimnames = TRUE)
      }
    }
    l
  }
  
  # Helper which can save a list of SummarizedExperiment with a recursive call
  .save_sce_list <- function(l,h5_file,group,sparse.assays) {
    h5createGroup(h5_file,group) |> stopifnot()
    h5write(as.character(names(l)),h5_file,paste0(group,"/names"))
    for(i in seq_along(l)) {
      saveHDF5SingleCellExperiment(l[[i]],h5_file,paste0(group,"/element_",i),force=TRUE,sparse.assays=sparse.assays)
    }
  }

  # Check params, and initialize H5 file
  stopifnot(is(x,"SummarizedExperiment"))
  if (!file.exists(h5_file)) {
    h5createFile(h5_file)
  } else {
    if (!force) stop("output file already exists, use force=TRUE to update")
  }
  h5createGroup(h5_file,group) || stop("cannot create main h5 group")
  
  # Save assays to H5 and remove them from the SCE object
  .h5_save_array_list(assays(x),h5_file,paste0(group,'/assays'),as.sparse = sparse.assays)
  assays(x) <- list()
  
  if (is(x,"SingleCellExperiment")) {
    # Save reducedDims to H5 and remove them from the SCE object
    .h5_save_array_list(reducedDims(x),h5_file,paste0(group,'/reducedDims'),as.sparse = FALSE)
    reducedDims(x) <- list()
    
    # Save altExps to H5 and remove from the SCE object
    .save_sce_list(altExps(x),h5_file,paste0(group,'/altExps'),sparse.assays)
    altExps(x) <- list()
  }
  
  # Serialize the remaining object within the H5
  h5write(serialize(x,NULL),h5_file,paste0(group,"/rds"))
}


#' Load a SingleCellExperiment object from a H5
#' @param h5_file name of the h5 file to load
#' @param group name of the group within the h5 containing the data
#' @return the SingleCellExperiment object
loadHDF5SingleCellExperiment <- function(h5_file,group="/main") {
  
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
  
  .load_sce_list <- function(h5_file,group) {
    elt_names <- as.character(h5read(h5_file,paste0(group,"/names")))
    l <- list()
    for(i in seq_along(elt_names)) {
      l[[i]] <- loadHDF5SingleCellExperiment(h5_file,paste0(group,"/element_",i))
    }
    names(l) <- elt_names
    l
  }
  
  x <- unserialize(h5read(h5_file,paste0(group,"/rds")))
  assays(x) <- .h5_load_array_list(h5_file,paste0(group,"/assays"))
  if (is(x,"SingleCellExperiment")) {
    reducedDims(x) <- .h5_load_array_list(h5_file,paste0(group,"/reducedDims"))
    altExps(x) <- .load_sce_list(h5_file,paste0(group,"/altExps"))
  }
  
  x
}

# Unit test
function() {
  E <- DropletUtils::read10xCounts("../../../10x/filtered_feature_bc_matrix.h5")
  E <- E[,1:30]
  reducedDim(E,"PCA") <- matrix(runif(ncol(E)*5),ncol(E))
  reducedDim(E,"UMAP") <- matrix(runif(ncol(E)*2),ncol(E))
  E <- splitAltExps(E,mcols(E)$Type,ref = "Gene Expression")
  saveHDF5SingleCellExperiment(E,"out/test2.h5")
  h5ls("out/test2.h5")
  
  e <- loadHDF5SingleCellExperiment("out/test2.h5")
  plot(reducedDim(e,"UMAP"))
}







