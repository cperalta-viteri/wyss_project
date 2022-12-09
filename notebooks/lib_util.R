

#' group rows of given matrix m, and compute average per group
rowmean <- function(x,group,...) {
    rowsum(x,group,...) / rowsum(1L-is.na(x),group,...)
}


