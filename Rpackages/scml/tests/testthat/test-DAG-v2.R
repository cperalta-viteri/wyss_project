

library(tidyverse)
library(tidygraph)
library(ggraph)
library(igraph)
library(torch)
library(luz)
#library(scml)



#-#-#-#-#-#-#-#-#-#-#-#-#
# Definition of the Loss
#-#-#-#-#-#-#-#-#-#-#-#-#

# Convert a sparse binary matrix M into a tensor
as_tensor.Matrix <- function(M) {
  i <- Matrix::which(M,arr.ind=TRUE)
  torch_sparse_coo_tensor(indices = t(i),value = M[i],size = dim(M))
}

# Extract ancestor matrix from graph h
adj_ancestors <- function(h) {
  A <- igraph::connect(h,+Inf,mode="out")
  A[from=V(A),to=V(A)] <- TRUE
  as_tensor.Matrix(igraph::as_adj(A,names = FALSE)!=0)
}

# Extract child -> parent pairs from graph h
adj_matrix <- function(h) {
  as_tensor.Matrix(igraph::as_adj(h,names = FALSE)!=0)
}


#' @param input a (M,N) float tensor tensor
#' @param target a 1D long tensor
#' @param P Parent matrix: a square boolean 2D tensor
#' @param A Ancestor matrix: a square boolean 2D tensor
nnf_dag_multi_margin_loss <- function(input,target,P,A) {
  Y <- A$index_select(1L,target)

  # Make a 3D tensor
  P <- P$to_dense()
  Y <- Y$to_dense()
  POS <- P$unsqueeze(1L) & Y$unsqueeze(-1L) & Y$unsqueeze(2L)
  NEG <- P$unsqueeze(1L) & (!Y)$unsqueeze(-1L) & Y$unsqueeze(2L)

  pos <- (input$unsqueeze(-1L) * POS)
  pos[!POS] <- +Inf
  pos <- pos$amin(dim=2L)

  neg <- (input$unsqueeze(-1L) * NEG)
  neg[!NEG] <- -Inf
  neg <- neg$amax(dim=2L)

  loss <- (neg - pos + 1)$maximum(0)
  # TODO: scale losses according to tree depth
  # TODO: scale losses according to childs cardinality (balancing)

  loss <- loss$amax(dim=2L) # reduce over tree
  loss <- loss$mean() # reduce over samples
  loss
}




#-#-#-#-#-#-#-#-#-#-#-#-#
# Main
#-#-#-#-#-#-#-#-#-#-#-#-#

# Define label hierarchy
h <- create_tree(16,3,mode = "in")

# Display the hierarchy
mutate(h,name=row_number()) |>
  ggraph("tree",mode="in") +
  geom_edge_link(arrow = grid::arrow(10,type = "closed")) +
  geom_node_label(aes(label=name))

# Test loss
y_true <- torch_tensor(c(6L,16L,12L),torch_long())
y_pred <- torch_rand(length(y_true),igraph::vcount(h),requires_grad = TRUE)
nnf_dag_multi_margin_loss(y_pred,y_true,adj_matrix(h),adj_ancestors(h))




