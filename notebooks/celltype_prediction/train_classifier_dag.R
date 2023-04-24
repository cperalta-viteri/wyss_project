
library(torch)
library(luz)
library(ggraph)
library(tidygraph)
source("notebooks/celltype_prediction/lib_cell_classifier.R")
source("notebooks/01_data_collection/lib_sce.R")





label_hierarchy <- function(lev) {
  pat <- "(((.*) --- .*) --- .*)"
  e <- rbind(
    cbind("root",sub(pat,"\\3",lev)),
    cbind(sub(pat,"\\3",lev),sub(pat,"\\2",lev))
    #cbind(sub(pat,"\\2",lev),sub(pat,"\\1",lev))
  ) |> unique()
  v <- unique(as.vector(e))
  e <- array(match(e,v),dim(e))
  Matrix::sparseMatrix(i=e[,2L],j=e[,1L],dims = c(length(v),length(v)),dimnames = list(v,v))
}

# extract ancestors matrix
ancestors <- function(P) {
  for(i in seq(log2(ncol(P))+1L)) {
    P <- (P%*%P) | P
  }
  P
}

# allenSSq <- loadHDF5SingleCellExperiment("notebooks/01_data_collection/allen_brain_atlas/out/mouse_ctx_ssv4.h5");allenSSq$source <- "ssv4"
# allenSSq$y_true <- paste(allenSSq$class_label,allenSSq$neighborhood_label,allenSSq$subclass_label,allenSSq$cluster_label,sep=" --- ")
# allenSSq$y_true <- factor(allenSSq$y_true)
# fn <- rownames(allenSSq)
# lev <- levels(allenSSq$y_true)
# trn_ds <- sc_dataset(t(unname(counts(allenSSq))),y = as.integer(allenSSq$y_true))

trn_ds <- readRDS("notebooks/celltype_prediction/data/allen_subclass_trn.rds")
fn <- colnames(trn_ds$x)
lev <- levels(trn_ds$y)

# extract parent matrix from labels
P <- label_hierarchy(lev)
#P["GABAergic --- MGE --- Sst","Glutamatergic --- L2/3 IT"] <- TRUE

# extract ancestors matrix
A <- ancestors(P)




# Display hierarchy
as_tbl_graph(as.matrix(t(P))+0) |>
  mutate(label = sub(".* --- ","",name)) |>
  mutate(dfs_rank = dfs_rank(),dfs_dist=dfs_dist()) |>
  mutate(leaf = centrality_degree(mode="out")==0L) |>
  ggraph("manual",x=dfs_dist,y=dfs_rank) +
    geom_edge_diagonal(flipped = TRUE) +
    geom_node_label(aes(label=label,filter=!leaf),size=3) +
    geom_node_label(aes(label=label,filter=leaf),fill="green",size=3,hjust=0) +
    scale_y_reverse() + 
    expand_limits(x=c(1,4))


P <- torch_tensor(as.matrix(P),torch_float())
A <- torch_tensor(as.matrix(A),torch_float())


# Generate fake predition and targets
y_pred <- -(torch_tensor(1:5)[,NULL] + torch_tensor((seq(ncol(P))-1)*5)[NULL,])
y_true <- matrix(c(1:5,5L,c(9L,10L,7L,11L,7L,11L)),byrow=TRUE,nrow=2)
y_true <- torch_sparse_coo_tensor(y_true,c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE),size = dim(y_pred))


# propagate y_true labels to ancestors
y_true <- (y_true$to(torch_float())$matmul(A) + y_true)$to(torch_bool())



cst <- P$nonzero()             # Voici la liste des contraintes
cst_pred <- y_pred[,cst[,1L]]  # Recuperer les valeurs associées a chaque contrainte
cst_dir <- y_true[,cst[,1L]]   # Determiner la direction de chaque element

# reduce
cst_pos <- cst_pred$clone()
cst_neg <- cst_pred$clone()
cst_pos[!cst_dir] <- +Inf
cst_neg[cst_dir] <- -Inf
cst_pos <- torch_full(c(nrow(cst_pos),cst[,2L]$max()$item()),+Inf)$scatter_reduce(dim=2L,cst[,2L]$tile(c(5L,1L)),cst_pos,reduce="amin")
cst_neg <- torch_full(c(nrow(cst_neg),cst[,2L]$max()$item()),-Inf)$scatter_reduce(dim=2L,cst[,2L]$tile(c(5L,1L)),cst_neg,reduce="amax")



y_pred
y_true
cst_pos
cst_neg

o <- torch_arange(1,5)$unsqueeze(1) * torch_ones(10,5)
i <- torch_randint(1,3,5,dtype = torch_long())$reshape(c(1,5))$`repeat`(c(10,1))
torch_zeros_like(o)$scatter_reduce(2,i,src = o,reduce="amax")







# W <- P*colSums(Y) # and deduce node weights
# W[W>0] <- 1/W[W>0]


input <- torch_eye(5,3)
target <- matrix(c(
  2,3,-1,
  2,1,-1,
  2,1,-1,
  3,-1,-1,
  1,2,3
),5,byrow = TRUE) |> torch_tensor(torch_long())

nnf_multilabel_margin_loss(input,target,reduction = "none")



trn_ds <- tensor_dataset(
  x = torch_tensor(as.matrix(trn_ds$x)),
  y = torch_tensor(trn_ds$y,dtype = torch_long())
)



  






model <- advanced_cell_classifier |>
  setup(metrics = list(acc=luz_metric_accuracy())) |>
  set_hparams(feature_names=fn,class_names=lev) |>
  set_opt_hparams(lr1=1e-5) |>
  fit(trn_ds, epochs = 100L,valid_data = trn_ds,dataloader_options=list(batch_size=256L,shuffle=FALSE,num_workers=0L,worker_packages="HDF5Array"),callbacks = list(
    #luz_callback_model_checkpoint("notebooks/celltype_prediction/models/model_subclass_best_valid_acc.luz",monitor="valid_acc",mode="max",save_best_only=TRUE),
    #luz_callback_model_checkpoint("notebooks/celltype_prediction/models/model_subclass_best_train_acc.luz",monitor="train_acc",mode="max",save_best_only=TRUE),
    #luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"notebooks/celltype_prediction/models/model_subclass_best_valid_acc.luz")})()
  ))
luz_load_model_weights(model,"notebooks/celltype_prediction/models/model_subclass_best_valid_acc.luz")
msk <- model$model$compute_prune_mask()


model_pruned <- advanced_cell_classifier |>
  setup(metrics = list(acc=luz_metric_accuracy())) |>
  set_hparams(feature_names=model$model$feature_names,class_names=model$model$class_names,prune_mask=msk) |>
  set_opt_hparams(lr1=1e-5) |>
  fit(trn_ds, epochs = 20L,valid_data = trn_ds,dataloader_options=list(batch_size=256L,shuffle=TRUE),callbacks = list(
    luz_callback_model_checkpoint("notebooks/celltype_prediction/models/model_subclass_pruned_best_train_acc.luz",monitor="train_acc",mode="max",save_best_only=TRUE),
    luz_callback_model_checkpoint("notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.luz",monitor="valid_acc",mode="max",save_best_only=TRUE),
    luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"notebooks/celltype_prediction/models/model_subclass_best_valid_acc.luz")})()
    #luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"notebooks/celltype_prediction/models/model_subclass_pruned_best_train_acc.luz")})()
  ))

luz_load_model_weights(model_pruned,"notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.luz")
model_pruned$model(trn_ds[1:5]$x)
torch_save(model_pruned$model,"notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.torch")








