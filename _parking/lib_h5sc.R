
library(rhdf5)
library(HDF5Array)
library(SingleCellExperiment)



readH5 <- function(h5file,path="/") {
  
  handle_elt_ <- function(h5file,path,content) {
    message(sprintf("start processing element: %s",path))
    
    # Get attributes of the current node
    a <- h5readAttributes(h5file,path)
    message(sprintf("%s",a$"encoding-type"))
    
    # Call element handler
    switch(a$"encoding-type",
           "csr_matrix" = {
             H5SparseMatrix(h5file,path)
           },
           "dense_array" = {
             HDF5Array(h5file,path)
           },
           "factor" = {
             lev <- h5read(h5file,paste0(sub("/$","",path),"/levels"))
             i <- h5read(h5file,paste0(sub("/$","",path),"/codes"))
             factor(lev,lev)[i+1L]
           },
           "list" = {
             L <- setNames(content$path,content$name)[content$group %in% path] |> 
               lapply(handle_elt_,h5file=h5file,content=content)
           },
           "data_frame" = {
             D <- setNames(content$path,content$name)[content$group %in% path] |> 
               lapply(handle_elt_,h5file=h5file,content=content) |>
               as("DataFrame")
             rownames(D) <- as.character(D[[a$"_index"]])
             D <- D[a$"column-order"]
             D
           },
           "SingleCellExperiment" = {
             L <- setNames(content$path,content$name)[content$group %in% path] |> 
               lapply(handle_elt_,h5file=h5file,content=content)
             E <- SingleCellExperiment(list(X=L$X),colData=L$obs,rowData=L$var,reducedDims = lapply(L$obsm,t))
             L$X <- L$obs <- L$var <- L$obsm <- NULL
             metadata(E) <- L
             E
           },
           stop(sprintf("When processing %s cannot handle %s",path,handler_id))
    )
  }
  
  content <- h5ls(h5file)
  content$path <- paste0(sub("/?$","/",content$group),content$name)
  handle_elt_(h5file,path,content)  
}


writeH5 <- function(e,h5file,path="/x") {
  if (!file.exists(h5file)) h5createFile(h5file)

  if (is(e,"array") || is(e,"dgeMatrix")) {
    message(sprintf("write 'array' to %s",path))
    writeHDF5Array(e,h5file,path,with.dimnames = TRUE)
  } else if (is(e,"dgCMatrix")) {
    message(sprintf("write 'dgCMatrix' to %s",path))
    writeTENxMatrix(e,h5file,path,with.dimnames = TRUE)
  } else if (is(e,"list")) {
    message(sprintf("write 'list' to %s",path))
    h5createGroup(h5file,path)
    if (!is.null(names(e))) h5write(names(e),paste0(path,"/names"))
    for(i in seq_along(e)) writeH5(e[[i]],h5file,paste0(path,"/",i))
  } else if (is(e,"data.frame") || is(e,"DataFrame")) {
    
  } else {
    stop(sprintf("Cannot write '%s' to '%s'",class(e),path))
  }

  #h5f <- H5Fopen(h5file)
  #on.exit(H5Fclose(h5f))

  
  
  # } else if (is(e,"vector")) {
  # } else if (is(e,"factor")) {
  # } else if (is(e,"Rle")) {
  # } else if (is(e,"GRanges")) {
  #   #seqinfo
  # } else if (is(e,"DelayedArray")) {
  #   #writeHDF5Array()
  #   #writeTENxMatrix()
  
  # } else if (is(e,"SingleCellExperiment")) {    
}

e <- matrix(1:100,20,5)
h5file <- tempfile(fileext = ".h5")
path <- "/mat"
e <- list(e,e)
writeH5(e,h5file)
h5ls(h5file)

