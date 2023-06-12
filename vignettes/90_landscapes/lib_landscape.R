
library(IRanges)
library(ggplot2)
library(DelayedArray)
library(SummarizedExperiment)



rowmean <- function(M,f,na.rm=TRUE) {
  f <- as.factor(f)
  M <- rowsum(M,f)
  M[levels(f),,drop=TRUE] / as.vector(table(f))
}

#' return an integer array of nearest neighbors
knn_array <- function(x,y,k=11L,gx = seq(min(x),max(x),length.out=32L),gy = seq(min(y),max(y),length.out=32L)) {
  g <- matrix(NA,length(gy),length(gx))
  knn <- FNN::get.knnx(cbind(x,y),cbind(gx[col(g)],gy[row(g)]),k=k)$nn.index
  knn <- array(knn,c(dim(g),k))
  attr(knn,"x") <- gx
  attr(knn,"y") <- gy
  knn
}


#' Split rows into batches of 100 genes
#' Compute average expression of each gene in m according to index structure found in i
#'@param m an expression matrix (row=gene,col=cell)
#'@param i an array of column indexes, to be reduced over the last dimension
knn_rowMeans <- function(m,i) {
  B <- local({
    batches <- split(seq(nrow(m)),seq(nrow(m)) %/% 100L)
    B <- lapply(batches,function(b){
      M <- m[b,as.vector(i)]
      M <- array(M,c(length(b),dim(i)))
      DelayedArray(rowMeans(M,dims=length(dim(i))))
    })
    do.call(rbind,B)
  })       
}


# plot the average batch of gene
gg_landscape <- function(B,x=seq(dim(B)[3]),y=seq(dim(B)[2])) {
  D <- reshape2::melt(as.array(B))
  D$x <- x[D$Var3]
  D$y <- y[D$Var2]
  ggplot(D) + facet_wrap(~Var1) +
    geom_tile(aes(x=x,y=y,fill=value)) +
    scale_fill_gradientn(colors = c("#35978f", "#80cdc1", "#c7eae5", "#f5f5f5","#f6e8c3", "#dfc27d", "#bf812d", "#8c510a","#543005", "#330000")) +
    #scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11,"RdYlBu"))) +
    coord_equal() +
    #guides(fill="none") +
    xlab("") + ylab("") + theme(panel.grid = element_blank())
}

gg_landscape_image <- function(m,x=seq(dim(m)[2]),y=seq(dim(m)[1])) {
  D <- reshape2::melt(m)
  D$x <- x[D$Var2]
  D$y <- y[D$Var1]
  ggplot(D) + 
    geom_tile(aes(x=x,y=y,fill=value)) +
    scale_fill_gradientn(colors = c("#35978f", "#80cdc1", "#c7eae5", "#f5f5f5","#f6e8c3", "#dfc27d", "#bf812d", "#8c510a","#543005", "#330000")) +
    #scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(11,"RdYlBu"))) +
    coord_equal() +
    theme_void() + theme(legend.position = "none")
}


# generate a random dataset
random_SE <- function(ngene=5000,ncell=1000) {
  x <- runif(ncell)
  y <- runif(ncell)
  m <- matrix(runif(ngene*ncell),ngene,ncell)
  SummarizedExperiment(m,colData=DataFrame(x=x,y=y))
}






