

#' group rows of given matrix m, and compute average per group
rowmean <- function(x,group,...) {
    rowsum(x,group,...) / rowsum(0L + !is.na(x),group,...)
}


order.col <- function(m) {
  col(m)[order(row(m),m)] |>
    matrix(nrow = nrow(m),byrow = TRUE)
}
