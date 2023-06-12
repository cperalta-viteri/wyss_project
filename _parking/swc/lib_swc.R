


# Read a SWC file
read_swc <- function(f) {
  x <- read.table(f,sep=" ",comment="#",col.names = c("v1","sample_id","structure_id","x","y","z","radius","parent_id"))
  lev <- c("undefined","soma","axon","basal_dendrite","apical_dendrite","fork_point","end_point","custom")
  x$structure_label <- factor(lev,lev)[pmin(x$structure_id,7L) + 1L]
  x
}

# Compute density binary matrix of a SWC file
# - morpho is projected on y-axis
# - the z is assumed to encode cortical depth (centered on pia, z=0) and with maximum depth of 1000
# - the x is assumed to range from -500 to +500
swc_density <- function(swc,n=50L) {
  xy <- list(x=swc$x,y=swc$z)
  gx <- seq(-500,500,length.out=n+1L)
  gy <- seq(-1000,0,length.out=n+1L)
  z <- Matrix::sparseMatrix(
    findInterval(xy$x,gx,all.inside = TRUE),
    findInterval(xy$y,gy,all.inside = TRUE),
    dims = c(n,n)
  )
  list(x = (gx[-1] + gx[-length(gx)]) / 2,y=(gy[-1] + gy[-length(gy)]) / 2,z=as.matrix(z))
}
