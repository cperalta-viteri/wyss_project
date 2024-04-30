

library(tidygraph)
library(ggraph)
library(tidyverse)
library(torch)
library(luz)
library(scml)



#-#-#-#-#-#-#-#-#-#-#-#
# Example 1
#-#-#-#-#-#-#-#-#-#-#-#
e <- matrix(ncol = 2,byrow = TRUE,c(
  4L,5L,
  5L,1L,
  5L,2L,
  4L,3L
)) |> torch_tensor()

dag <- as_tbl_graph(as_array(e))
dag |>
  mutate(name=row_number()) |>
  ggraph() +
  geom_edge_link(arrow = arrow(angle = 10,type = "closed")) +
  geom_node_label(aes(label=name))

# ancestors matrix
a <- t(as.matrix(igraph::as_adjacency_matrix(igraph::connect(dag,+Inf,mode = "out"))))>0
diag(a) <- TRUE



m <- nn_linear |>
  setup(loss=nn_dag_loss(e,torch_ones(5)),optimizer = optim_adam) |>
  set_hparams(in_features=4,out_features=5) |>
  set_opt_hparams(lr=0.01) |>
  fit(data=list(data.matrix(iris[1:4]),a[iris$Species,]),accelerator=accelerator(cpu=TRUE),epoch=550,dataloader_options=list(batch_size=200,drop_last=FALSE))

p <- predict(m,data.matrix(iris[1:4]))
image(as_array(p))
table(iris$Species,max.col(as_array(p)))

table(iris$Species,max.col(as_array(p)[,1:2]))

