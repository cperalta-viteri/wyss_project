#'
#'
#'
#'
#' #' @param T a 3-columns integer matrix
#' #' @param C a numeric k x k cost matrix: C[i,j] is the loss for mispredicting a sample of class i into class j
#' nnf_dag_loss <- function(input,target,triplets,C) {
#'   P <- torch_rand(100,100)>0.7 # parent matrix
#'   C <- torch_rand(100,100) # cost matrix
#'   stopifnot(identical(nrow(C),ncol(C)))
#'   #A: Ancestors matrix
#'
#' }
#'
#'
#'
#'
#' # extract ancestors matrix
#' ancestors <- function(P) {
#'   for(i in seq(log2(ncol(P))+1L)) {
#'     P <- (P%*%P) | P
#'   }
#'   diag(P) <- TRUE
#'   P
#' }
#'
#' triplets <- function(P) {
#'   CST <- which(P,arr.ind = TRUE)
#'   CST <- merge(CST,CST,by=2)
#'   CST <- subset(CST,row.x<row.y)
#'   data.matrix(CST)
#' }
#'
#' nn_dag_loss <- nn_module(
#'   "nn_dag_loss",
#'
#'   #' @param P A logical square parent matrix: P[i,j] is TRUE when node j is the parent of node i
#'   #' @param C A numeric square cost matrix: C[i,j] is the loss for mispredicting a sample of class i into class j
#'   initialize = function(P,C) {
#'     # TODO: accept torch_tensor as input instead of Matrix
#'     stopifnot(identical(nrow(P),ncol(P)))
#'     stopifnot(identical(dim(P),dim(C)))
#'     self$P <- P
#'     self$C <- torch_tensor(C,torch_float())
#'     self$A <- torch_tensor(as.matrix(ancestors(P)),torch_bool()) # Extract ancestor matrix
#'     self$TRI <- torch_tensor(triplets(P)) # And triplets
#'     self$TRI_weight <- torch_gather(self$C[self$TRI[,1L],],2L,self$TRI[,2L][,NULL])$squeeze(-1L)
#'   },
#'   forward = function(y_pred, y_true) {
#'     # propagate y_true labels to their ancestors
#'     A <- self$A[y_true,]$to(torch_long())
#'
#'     # Extract triplet data for each input
#'     TRI_target1 <- A[..,self$TRI[,2L]]
#'     TRI_target2 <- A[..,self$TRI[,3L]]
#'     TRI_val1 <- y_pred[..,self$TRI[,2L]]
#'     TRI_val2 <- y_pred[..,self$TRI[,3L]]
#'     TRI_dir <- sign(TRI_target1 - TRI_target2)
#'     TRI_loss <- nnf_margin_ranking_loss(TRI_val1,TRI_val2,TRI_dir,margin=1,reduction="none") * self$TRI_weight
#'     TRI_loss[TRI_dir==0L] <- 0
#'
#'     # Reduce triplets losses
#'     #torch_full(c(length(y_true),ncol(self$P)),-Inf)$scatter_reduce(dim=-1L,self$TRI[,1L]$tile(c(length(y_true),1L)),TRI_loss,"amax")
#'     TRI_loss$amax(-1L)$mean()
#'   }
#' )
#'
#'
#'
#' trn_ds <- readRDS("notebooks/celltype_prediction/data/allen_subclass_trn.rds")
#'
#' # extract parent matrix from labels
#' P <- label_hierarchy(levels(trn_ds$y))
#'
#' # Display hierarchy
#' as_tbl_graph(as.matrix(t(P))+0) |>
#'   mutate(label = sub(".* --- ","",name)) |>
#'   mutate(dfs_rank = dfs_rank(),dfs_dist=dfs_dist()) |>
#'   mutate(leaf = centrality_degree(mode="out")==0L) |>
#'   ggraph("manual",x=dfs_dist,y=dfs_rank) +
#'     geom_edge_diagonal(flipped = TRUE) +
#'     geom_node_label(aes(label=label,filter=!leaf),size=3) +
#'     geom_node_label(aes(label=label,filter=leaf),fill="green",size=3,hjust=0) +
#'     scale_y_reverse() +
#'     expand_limits(x=c(1,4))
#'
#'
#'
#'
#'
#'
#' label_hierarchy <- function(lev) {
#'   pat <- "(((.*) --- .*) --- .*)"
#'   e <- rbind(
#'     cbind("root",sub(pat,"\\3",lev)),
#'     cbind(sub(pat,"\\3",lev),sub(pat,"\\2",lev))
#'     #cbind(sub(pat,"\\2",lev),sub(pat,"\\1",lev))
#'   ) |> unique()
#'   v <- unique(as.vector(e))
#'   e <- array(match(e,v),dim(e))
#'   Matrix::sparseMatrix(i=e[,2L],j=e[,1L],dims = c(length(v),length(v)),dimnames = list(v,v))
#' }
#'
#'
#' dag_labeller <- nn_module(
#'   "DAG_labeller",
#'   initialize = function(P,feature_names) {
#'     self$loss_fn <- nn_dag_loss(P,C=1-torch_eye(nrow(P)))
#'     self$layers <- nn_sequential(
#'       nn_dropout(),
#'       nn_linear(length(feature_names),4096L),nn_tanh(),nn_dropout(0.3),
#'       nn_linear(4096L,1024L),nn_relu(),nn_dropout(0.2),
#'       nn_linear(1024L,ncol(P))
#'     )
#'   },
#'   forward = function(x) {
#'     x |> self$layers()
#'   },
#'   loss = function(...) {self$loss_fn(...)}
#' )
#'
#'
#'
#' fm <- dag_labeller |>
#'   setup(optimizer = optim_adam) |>
#'   set_opt_hparams(lr=1e-3,weight_decay=0) |>
#'   set_hparams(P = P) |>
#'   fit(sc_dataset(trn_ds$x,trn_ds$y),
#'       epochs = 200L,
#'       dataloader_options = list(batch_size=32L,drop_last=FALSE),
#'       accelerator = accelerator(cpu = TRUE)
#'   )
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
