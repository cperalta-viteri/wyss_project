
library(rhdf5)
library(HDF5Array)
library(SingleCellExperiment)


readH5AD <- function(h5file="alz_excitatory.h5ad",path="/") {
  message(sprintf("start processing element: %s",path))
  
  # Get attributes
  a <- h5readAttributes(h5file,path)
  message(sprintf("  %s [v%s]",a$"encoding-type",a$"encoding-version"))
  
  # List sub-elements
  subelts <- h5ls(h5file) |> subset(group %in% path)
  subelts$path <- paste0(sub("/?$","/",subelts$group),subelts$name)
  
  # Define Handlers and Call appropirate one
  handlers <- alist(
    "anndata[0.1.0]" = {
      L <- lapply(subelts$path,readH5AD,h5file=h5file) |> 
        setNames(subelts$name)
      E <- SingleCellExperiment(list(X=L$X),colData=L$obs,rowData=L$var)
      L$X <- L$obs <- L$var <- NULL
      metadata(E) <- L
      E
    },
    "csr_matrix[0.1.0]" = {
      H5SparseMatrix(h5file,path)
    },
    "dict[0.1.0]" = {
      L <- lapply(subelts$path,readH5AD,h5file=h5file) |> 
        setNames(subelts$name)
      L
    },
    "dataframe[0.2.0]" = {
      D <- lapply(subelts$path,readH5AD,h5file=h5file) |> 
        setNames(subelts$name) |>
        as("DataFrame")
      rownames(D) <- as.character(D[[a$"_index"]])
      D <- D[a$"column-order"]
      D
    },
    "categorical[0.2.0]" = {
      lev <- h5read(h5file,paste0(sub("/$","",path),"/categories"))
      i <- h5read(h5file,paste0(sub("/$","",path),"/codes"))
      factor(lev,lev)[i+1L]
    },
    "string[0.2.0]" = {
      as.character(h5read(h5file,path))
    },
    "string-array[0.2.0]" = {
      A <- h5read(h5file,path)
      if (length(dim(A))==1L) A <- as.vector(A)
      A
    },
    "array[0.2.0]" = {
      A <- h5read(h5file,path)
      if (length(dim(A))==1L) A <- as.vector(A)
      A
    },
    "raw[0.1.0]" = {
      L <- lapply(subelts$path,readH5AD,h5file=h5file) |> 
        setNames(subelts$name)
      SummarizedExperiment(list(X=L$X),rowData=L$var,metadata=L$uns)
    }
  )
  
  # Check attributes
  tryCatch({
    # Call element handler
    handler_id <- sprintf("%s[%s]",a$"encoding-type",a$"encoding-version")
    stopifnot("Cannot handle object" = handler_id %in% names(handlers))
    eval(handlers[[handler_id]])
  },
  error = function(e) {
    print(a)
    stop(e$message)
  }
  )
}





