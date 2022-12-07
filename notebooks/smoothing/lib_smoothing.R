



knn_smooth_2d <- function(x,y,z,k=0) {
    x <- as.vector(x)
    y <- as.vector(y)
    z <- as.vector(z)
    if (k>0) {
        knn <- FNN::knnx.index(cbind(x,y),cbind(x,y),k=k)
        z <- rowMeans(cbind(z,array(z[knn],dim(knn))))
    }
    data.frame(x=x,y=y,z=z)
}


