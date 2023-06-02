


#' Read json file storing the allen taxonomy and convert into a data.frame
#' @param json_file to the json file containing the taxonomy in allen20 format
#' @return a data.frame
#' @importFrom rjson fromJSON
#' @import utils
#' @import methods
aba_read_dend_json <- function(json_file) {
  x <- rjson::fromJSON(file=json_file,simplify = FALSE)

  # recursive function converting into dcf format
  to_dcf <- function(n,path="") {
    # Check if the node is either a leaf or an internal node
    stopifnot(xor(setequal(names(n),"leaf_attributes"),setequal(names(n),c("children","node_attributes"))))

    # retrieve attributes
    att <- n[[grep("_attributes",names(n))]][[1]]
    stopifnot(utils::hasName(att,"cell_set_accession"))
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


.download_url <- function(url,dir) {
  local_path <- file.path(dir,basename(url))
  download.file(url, local_path, method = "wget")
  local_path
}

#' @import scuttle
#' @importFrom BiocParallel MulticoreParam
.add_qc_metrics <- function(e) {
  BPPARAM <- BiocParallel::MulticoreParam(6,progressbar = TRUE)
  e <- scuttle::computeLibraryFactors(e,BPPARAM=BPPARAM)
  e <- scuttle::addPerCellQCMetrics(e,BPPARAM=BPPARAM)
  e <- scuttle::addPerFeatureQCMetrics(e,BPPARAM=BPPARAM)
  e
}


.read_sparseMatrix_from_csv <- function(file,batch_size=1024L) {
  con <- gzfile(file,"r")

  # Create the initial Matrix as an empty row, and set colnames
  M <- read.csv(con,nrows=batch_size,header=TRUE,row.names=1L,check.names = FALSE) |>
    data.matrix()  |> methods::as("RleArray")
  M <- list(M)

  # iterate over rows of the CSV
  while(nrow(M[[length(M)]])==batch_size) {
    message(sprintf("%d rows processed",sum(sapply(M,nrow))))
    M[[length(M)+1L]] <- read.csv(con,nrows=batch_size,header=FALSE,row.names=1L,check.names = FALSE) |>
      data.matrix() |> as("RleArray")
  }
  M <- do.call(rbind,M)
  M
}


