


#' Identify contribution of every components to a multiplexed profile 
#' 
#' @param x (n x m) expression profiles of m components
#' @param y a vector with a multiplexed expression profile of length n
#' @return a vector of m positive coefficients summing to one.
#' 
#' @example 
#' x <- matrix(runif(1000),ncol=5)
#' alpha_true <- c(0.3,0.1,0.1,0.4,0.1)
#' y <- x %*% alpha_true + rnorm(nrow(x),0,0.03)
#' plot(x %*% alpha_true,y,xlim=c(0,1),ylim=c(0,1));abline(0,1)
#' alpha_pred <- demux(x,y)
#' plot(x %*% alpha_pred,y,xlim=c(0,1),ylim=c(0,1));abline(0,1)
demux <- function(x,y) {
    stop("this version seems to not work, please fix")
    x <- as.matrix(x)
    y <- as.vector(y)
    stopifnot(identical(nrow(x),length(y)))
    w <- LowRankQP::LowRankQP(
        Vmat = diag(colSums(x)^2),
        dvec = -colSums(y * x),
        Amat = matrix(1,1,ncol(x)),
        bvec = 1,
        uvec = rep(1,ncol(x))
    )
    as.vector(w$alpha)
}

#min_a sum_i (a.x_i - y_i)^2 = min_a sum_i (a.x_i)^2 - 2.a.x_i.y_i + y_i^2
#                            = min_a sum_i (a.x_i)^2 - 2.a.x_i.y_i







