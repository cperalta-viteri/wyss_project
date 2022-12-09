


library(ggplot2)
library(ggraph)
library(igraph)
library(IRanges)
library(tidygraph)



Rcpp::cppFunction(plugins="cpp11","
NumericVector hclust_reduce(const List& h,const NumericVector& v,char f='+') {
  IntegerMatrix m = h[\"merge\"];
  if (v.length() != m.nrow()+1) Rcpp::stop(\"vector length mismatch\");
  NumericVector s(m.nrow(),0.0);
  for(int i=0;i<m.nrow();i++) {
    auto v0 = (m(i,0)<0)?v[-m(i,0)-1]:s[m(i,0) - 1];
    auto v1 = (m(i,1)<0)?v[-m(i,1)-1]:s[m(i,1) - 1];
    if (f=='+') { //sum
        s[i] = v0 + v1; 
    } else if (f=='a') { //average
        s[i] = (v0 + v1)/2;
    } else if (f=='h') { //height
        s[i] = std::max(v0,v1) + 1;
    }
  }
  return s;
}")


hclust_ancestors <- function(h) {
  A <- Matrix::sparseMatrix(
    h$merge[h$merge>0],
    row(h$merge)[h$merge>0],
    dims = c(nrow(h$merge),nrow(h$merge))
  )
  Matrix::diag(A) <- TRUE
  for(k in seq(nrow(A))) {
    B <- (A %*% A) | A
    if (identical(A,B)) break;
    A <- B
  }
  A
}

hclust_leafs <- function(h) {
  leafs <- data.frame(
    id = seq_along(h$order),
    parent_id = NA_integer_,
    label = NA_character_,
    height = 0,
    order = NA_integer_
  )
  leafs$parent_id[-h$merge[h$merge<0]] <- row(h$merge)[h$merge<0]
  leafs$order[h$order] <- seq_along(h$order)
  leafs
}

hclust_nodes <- function(h) {
  nodes <- data.frame(
    id = seq(nrow(h$merge)),
    parent_id = NA_integer_,
    label = NA_character_,
    height = h$height,
    child1 = h$merge[,1],
    child2 = h$merge[,2]
  )
  nodes$parent_id[h$merge[h$merge>0]] <- row(h$merge)[h$merge>0]
  nodes
}


hclust_cut <- function(h,k) {
  # extract leafs/nodes basic metadata
  leafs <- hclust_leafs(h)
  nodes <- hclust_nodes(h)
  
  # get ancestors matrix
  A <- hclust_ancestors(h)
  a <- which(A[leafs$parent_id,],arr.ind = TRUE)
  nodes$depth <- rowSums(A)-1L
  
  # Assign clusters to nodes
  leafs$cluster <- cutree(h,k)
  nodes$clusters <- as.matrix(Matrix::sparseMatrix(a[,2],leafs$cluster[a[,1]],x=1))
  nodes$cluster <- max.col(nodes$clusters)
  is.na(nodes$cluster) <- rowSums(nodes$clusters>0)!=1
  
  nodes$num_leaf <- rowSums(nodes$clusters)
  
  
  # Assign a node type
  nodes$type <- ifelse(rowSums(nodes$clusters>0)>1,"mix","unif")
  leafs$type <- "unif"
  leafs$type[nodes$type[leafs$parent_id] %in% "mix"] <- "mix_child"
  nodes$type[(nodes$type %in% "unif") & (nodes$type[nodes$parent_id] %in% "mix")] <- "unif_root"
  
  
  # Compute tree layout
  leafs$x <- leafs$order
  leafs$y <- leafs$height
  nodes$x <- hclust_reduce(h,leafs$x,'a')
  nodes$y <- nodes$height
  leafs$parent_x <- nodes$x[leafs$parent_id]
  leafs$parent_y <- nodes$y[leafs$parent_id]
  nodes$parent_x <- nodes$x[nodes$parent_id]
  nodes$parent_y <- nodes$y[nodes$parent_id]
  
  list(leafs=leafs,nodes=nodes,A=A)
}


# Exampel usage
function() {
  # Compute hierarchical clustering
  h <- hclust(dist(data.matrix(iris[1:4])))
  
  # extract leafs/nodes basic metadata
  leafs <- hclust_leafs(h)
  nodes <- hclust_nodes(h)
  
  # get ancestors matrix
  A <- hclust_ancestors(h)
  a <- which(A[leafs$parent_id,],arr.ind = TRUE)
  nodes$depth <- rowSums(A)-1L
  
  # Assign clusters to nodes
  leafs$cluster <- cutree(h,6)
  nodes$clusters <- as.matrix(Matrix::sparseMatrix(a[,2],leafs$cluster[a[,1]],x=1))
  nodes$cluster <- max.col(nodes$clusters)
  is.na(nodes$cluster) <- rowSums(nodes$clusters>0)!=1
  
  # Compute tree layout
  leafs$x <- leafs$order
  leafs$y <- leafs$height
  nodes$x <- hclust_reduce(h,leafs$x,'a')
  nodes$y <- nodes$height
  
  # Display tree layout
  ggplot() + 
    #geom_point(aes(x=x,y=y),data=nodes) +
    geom_segment(aes(x=x,y=y,xend=x,yend=y[parent_id]),data=nodes) +
    geom_segment(aes(x=x,y=y[parent_id],xend=x[parent_id],yend=y[parent_id]),data=nodes) +
    geom_segment(aes(x=x,y=y,xend=x,yend=nodes$y[parent_id]),data=leafs) +
    geom_segment(aes(x=x,y=nodes$y[parent_id],xend=nodes$x[parent_id],yend=nodes$y[parent_id]),data=leafs) +
    geom_point(aes(x=x,y=y,color=factor(cluster)),data=leafs) +
    geom_point(aes(x=x,y=y,color=factor(cluster)),data=nodes)
}




# 
# 
# # Compute size after each merge of a hierachical community
# Rcpp::cppFunction(plugins="cpp11",includes='#include <stack>','
# DataFrame hclust_layout(const IntegerMatrix& m,int vcount) {
# 
#     // Iterate over merge from leaf nodes
#     IntegerVector n(vcount + m.nrow(),1);
#     IntegerVector h(vcount + m.nrow(),0);
#     IntegerVector p(vcount + m.nrow(),NA_INTEGER);
#     for(int i=0;i<m.nrow();i++) {
#       n[vcount + i] = n[m(i,0) - 1] + n[m(i,1) - 1];
#       h[vcount + i] = std::max(h[m(i,0) - 1],h[m(i,1) - 1]) + 1;
#       p[m(i,0) - 1] = p[m(i,1) - 1] = vcount + i + 1;
#     }
#     
#     // Iterate over merge from root node
#     IntegerVector lhs(vcount + m.nrow(),0);
#     IntegerVector d(vcount + m.nrow(),0);
#     for(int i=m.nrow()-1;i>=0;i--) {
#       lhs[m(i,0) - 1] += lhs[vcount + i];
#       lhs[vcount + i] += n[m(i,0) - 1];
#       lhs[m(i,1) - 1] += lhs[vcount + i];
#       d[m(i,0) - 1] = d[m(i,1) - 1] = d[vcount + i] + 1;
#     }
#     
#     // DFS walk from root node
#     IntegerVector o(vcount + m.nrow(),0);
#     std::stack<int> q;
#     q.push(vcount + m.nrow() - 1);
#     for(int k=1;!q.empty();k++) {
#       auto e = q.top();
#       q.pop();
#       o[e] = k;
#       if (e>=vcount) {
#         q.push(m(e-vcount,0)-1);
#         q.push(m(e-vcount,1)-1);
#       }
#     }
#     
#     
#     return DataFrame::create(
#       Named("size") = n,
#       Named("height") = h,
#       Named("lhs") = lhs,
#       Named("parent") = p,
#       Named("depth") = d,
#       Named("order") = o
#     );
# }')




