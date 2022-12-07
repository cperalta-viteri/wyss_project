library(Matrix)

#' Read json file storing the allen taxonomy and convert into a data.frame
#' @param path to the json file containing the taxonomy in allen20 format
#' @return a data.frame
read_dend_json <- function(json_file) {
  x <- rjson::fromJSON(file=json_file,simplify = FALSE)
  
  # recursive function converting into dcf format
  to_dcf <- function(n,path="") {
    # Check if the node is either a leaf or an internal node
    stopifnot(xor(setequal(names(n),"leaf_attributes"),setequal(names(n),c("children","node_attributes"))))
    
    # retrieve attributes
    att <- n[[grep("_attributes",names(n))]][[1]]
    stopifnot(hasName(att,"cell_set_accession"))
    att$path <- paste0(path,"/",att$cell_set_accession) # add path
    
    l <- mapply(to_dcf,n$children,att$path)
    c(paste(names(att),att,sep=":",collapse="\n"),unlist(l))
  }
  
  dcf <- to_dcf(x)
  dcf <- read.dcf(textConnection(paste(dcf,collapse="\n\n")))
  dcf <- as.data.frame(dcf)

  # write and reparse the data.frame to convert numeric values
  local({
    con <- textConnection(NULL,"w")
    write.table(dcf,file=con,sep="\t",quote=FALSE,na = "")
    read.table(text=textConnectionValue(con),sep="\t",quote="",comment="",check.names = FALSE,stringsAsFactors = FALSE)
  })
}


#' Read rows of a count matrix from a connection to a CSV
#'
#' The first column should be a character with the rowname.
#' The CSV is expected to contain integers, and the function must be adapted for numeric values.
#'
#' @param con a connection to the CSV lines to reads
#' @param nrow desired number of row to read
#' @param ncol number of column (must be specified)
#' @return an integer matrix of size (nrow x ncol), possibly truncated if end of stream is reached. 
#'         The returned matrix has 0-rows if the function is further called after after end-of-stream.
scan_csv_matrix <- function(con,nrow,ncol) {
  m <- matrix(0L,nrow,ncol)
  rownames(m) <- seq_len(nrow)
  i <- 0L
  while(i<nrow) {
    row_name <- scan(con,character(),n = 1L,sep = ",",quiet = TRUE) # read the rowname
    if (length(row_name)==0) break # stop if rowname cannnot be read because of end-of-stream
    i <- i + 1L;
    rownames(m)[i] <- row_name
    
    n <- scan(con,integer(),n = ncol,sep = ",",nlines = 1L,quiet = TRUE)  # read the integer values of the row
    m[i,] <- n
  }
  m[seq_len(i),,drop=FALSE]
}




read_sparseMatrix_from_csv <- function(file,batch_size=1024L) {
  con <- gzfile(file)
  open(con,"r")
  
  # Create the initial Matrix as an empty row, and set colnames
  M <- local({
    # read header line of the CSV (remove first element)
    col_names <- scan(con,sep = ",",nlines=1L,what = character())[-1]
    Matrix(integer(),nrow=0,ncol=length(col_names),dimnames = list(NULL,col_names))
  })
  M <- list(M)
  
  # iterate over rows of the CSV
  while(TRUE) {
    m <- scan_csv_matrix(con,batch_size,ncol(M[[1]]))
    if (nrow(m)==0L) break # stop if no data were read because of the end-of-stream
    M[[length(M)+1L]] <- Matrix(m)
    print(sum(sapply(M,nrow))) # display progression
  }
  M <- do.call(rbind,M)
  M
}

