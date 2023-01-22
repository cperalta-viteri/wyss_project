
library(rhdf5)
library(HDF5Array)
library(SingleCellExperiment)


readH5AD <- function(h5file,path="/") {
  
  handle_h5ad_elt <- function(h5file,path,content) {
    message(sprintf("start processing element: %s",path))
    
    # Get attributes of the current node
    a <- h5readAttributes(h5file,path)
    message(sprintf("%s[v%s]",a$"encoding-type",a$"encoding-version"))

    # Call element handler
    handler_id <- sprintf("%s=%s",a$"encoding-type",a$"encoding-version")
    switch(handler_id,
      "csr_matrix=0.1.0" = {
         H5SparseMatrix(h5file,path)
      },
      "string=0.2.0" = {
        as.character(h5read(h5file,path))
      },
      "string-array=0.2.0" = {
        A <- h5read(h5file,path)
        if (length(dim(A))==1L) A <- as.vector(A)
        A
      },
      "array=0.2.0" = {
        A <- h5read(h5file,path)
        if (length(dim(A))==1L) A <- as.vector(A)
        A
      },
      "categorical=0.2.0" = {
        lev <- h5read(h5file,paste0(sub("/$","",path),"/categories"))
        i <- h5read(h5file,paste0(sub("/$","",path),"/codes"))
        factor(lev,lev)[i+1L]
      },
      "dict=0.1.0" = {
         L <- setNames(content$path,content$name)[content$group %in% path] |> 
            lapply(handle_h5ad_elt,h5file=h5file,content=content)
      },
      "dataframe=0.2.0" = {
         D <- setNames(content$path,content$name)[content$group %in% path] |> 
           lapply(handle_h5ad_elt,h5file=h5file,content=content) |>
           as("DataFrame")
         rownames(D) <- as.character(D[[a$"_index"]])
         D <- D[a$"column-order"]
         D
      },
      "anndata=0.1.0" = {
        L <- setNames(content$path,content$name)[content$group %in% path] |> 
          lapply(handle_h5ad_elt,h5file=h5file,content=content)
        E <- SingleCellExperiment(list(X=L$X),colData=L$obs,rowData=L$var)
        L$X <- L$obs <- L$var <- NULL
        metadata(E) <- L
        E
      },
      "raw=0.1.0" = {
         L <- setNames(content$path,content$name)[content$group %in% path] |> 
           lapply(handle_h5ad_elt,h5file=h5file,content=content)
         E <- SummarizedExperiment(list(X=L$X),rowData=L$var)
         L$X <- L$var <- NULL
         metadata(E) <- L
         E
      },
      stop(sprintf("When processing %s cannot handle %s",path,handler_id))
    )
  }
  
  content <- h5ls(h5file)
  content$path <- paste0(sub("/?$","/",content$group),content$name)
  handle_h5ad_elt(h5file,path,content)
}


writeH5AD <- function(e,h5file) {
  stop("not implemented yet")
}





