#Sys.setenv(PYTORCH_ENABLE_MPS_FALLBACK=1)
#Sys.getenv("PYTORCH_ENABLE_MPS_FALLBACK")

library(torch)
library(luz)
source("notebooks/celltype_prediction/lib_cell_classifier.R")

ds <- readRDS("~/Downloads/ml_dataset.rds")
fn <- colnames(ds$x)
ml_ds <- tensor_dataset(
  x = torch_tensor(ds$x),
  y = torch_tensor(ds$y+1L,dtype = torch_long())
)

model <- advanced_cell_classifier |>
  setup(metrics = list(acc=luz_metric_accuracy())) |>
  set_hparams(feature_names=fn,class_names=c("uninfected","infected")) |>
  set_opt_hparams(lr1=1e-5) |>
  fit(ml_ds, epochs = 50L,valid_data = ml_ds,dataloader_options=list(batch_size=128L),callbacks = list(
    #luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"model_best_valid_acc.luz")})(),
    luz_callback_model_checkpoint("model_best_valid_acc.luz",monitor="valid_acc",mode="max",save_best_only=TRUE),
    luz_callback_model_checkpoint("model_best_train_acc.luz",monitor="train_acc",mode="max",save_best_only=TRUE)
  ))
luz_load_model_weights(model,"model_best_train_acc.luz")
msk <- model$model$compute_prune_mask()

model_pruned <- advanced_cell_classifier |>
  setup(metrics = list(acc=luz_metric_accuracy())) |>
  set_hparams(feature_names=fn,class_names=c("uninfected","infected"),prune_mask=msk) |>
  set_opt_hparams(lr1=1e-6) |>
  fit(ml_ds, epochs = 50L,valid_data = ml_ds,dataloader_options=list(batch_size=128L),callbacks = list(
    luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"model_pruned_best_train_acc.luz")})(),
    luz_callback_model_checkpoint("model_pruned_best_train_acc.luz",monitor="train_acc",mode="max",save_best_only=TRUE),
    luz_callback_model_checkpoint("model_pruned_best_valid_acc.luz",monitor="valid_acc",mode="max",save_best_only=TRUE)
  ))

luz_load_model_weights(model_pruned,"model_pruned_best_train_acc.luz")
torch_save(model_pruned$model,"model_pruned_best_train_acc.torch")

y_true <- ml_ds$tensors$y |> as_array()  
y_pred <- model$model(ml_ds$tensors$x)$argmax(2L) |> as_array()
table(y_true,y_pred)


library(SingleCellExperiment)
x <- HDF5Array::loadHDF5SummarizedExperiment("~/Downloads/HC") |>
  as("SingleCellExperiment")
x$sum <- x$nCount_RNA
x$detected <- x$nFeature_RNA
x$subsets_mt_sum <- colSums(assay(x[grep("^mt-",rownames(x),value=TRUE)]))
x$subsets_rp_sum <- colSums(assay(x[grep("^Rp",rownames(x),value=TRUE)]))

x$infected <- NA
x$infected[grepl("POS",x$orig.ident)] <- TRUE
x$infected[grepl("NEG",x$orig.ident)] <- FALSE
x$infected[x$orig.ident %in% "NCCHC"] <- FALSE


x$qc <- "pass"
x$qc[x$nCount_RNA<2000] <- "low_read_count"
x$qc[x$nFeature_RNA<1000] <- "low_feature_count"
x$qc[x$mt_count/x$nCount_RNA>0.1] <- "high_mt"
x$qc[x$rp_count/x$nCount_RNA>0.1] <- "high_rp"



x$fold <- "tst"
x$fold[colnames(x) %in% rownames(ds$x)] <- "ml"
p <- torch_predict_cell_type(assay(x,"data"),"model_pruned_best_train_acc.torch")
x$pred <- p[,"infected"] - p[,"uninfected"]

table(x$pred_groups,x$qc=="pass")

table(cut(x$pred,c(-Inf,-1,0,1,+Inf),c("---","-","+","+++")),x$qc=="pass")


ggcells(x) + 
  geom_bar(aes(y=orig.ident,fill=qc),color="black")


#x <- x[,x$qc=="pass"]
# QC
ggcells(x) + 
  geom_violin(aes(y=orig.ident,x=sum,fill=infected),scale="width") + 
  geom_vline(xintercept=2000,linetype=2) + 
  scale_x_log10() + 
  xlab("read count")
ggcells(x) + 
  geom_violin(aes(y=orig.ident,x=detected,fill=infected),scale="width") + 
  geom_vline(xintercept=1000,linetype=2) + 
  xlab("num. detected gene")
ggcells(x) +
  geom_point(aes(x=sum,y=detected,color=qc)) + 
  geom_vline(xintercept=2000,linetype=2) +
  geom_hline(yintercept=1000,linetype=2) + 
  scale_x_log10()
ggcells(x) + 
  stat_summary(aes(x=orig.ident,y=subsets_mt_sum/sum,fill=infected)) + 
  #geom_hline(yintercept=0.1,linetype=2) + 
  scale_y_continuous(label=scales::percent) +
  ylab("% mitochondrial") +
  coord_flip()
ggcells(x) + 
  stat_summary(aes(x=orig.ident,y=subsets_rp_sum/sum,fill=infected)) + 
  #geom_hline(yintercept=0.1,linetype=2) + 
  scale_y_continuous(label=scales::percent) +
  ylab("% ribosomal protein") +
  coord_flip()




ggcells(x) + 
  facet_wrap(~fold) + 
  geom_violin(aes(y=orig.ident,x=pred,fill=infected),scale="width") + 
  geom_vline(xintercept=0)

ggcells(x[,x$fold=="tst"]) + 
  facet_wrap(~pred_groups) + 
  geom_violin(aes(y=orig.ident,x=pred,fill=infected),scale="width") + 
  geom_vline(xintercept=0)

table(x$pred_groups,(x$pred>0) == x$infected,x$fold,useNA="i")
table(x$infected,cut(x$pred,c(-Inf,0,+Inf)),x$fold,useNA="i")

table(x$orig.ident,cut(x$pred,c(-Inf,-1,0,1,+Inf),c("---","-","+","+++")),x$fold,useNA="i") |>
  as.data.frame() |> 
  subset(Var3 == "tst") |> 
  #subset(Var2 %in% c("---","+++")) |> 
  ggplot() + 
    geom_col(aes(y=Var1,x=Freq,fill=Var2),position = "fill",color="black") + 
    geom_vline(xintercept=0.3,linetype=2) + 
    scale_x_continuous(label=scales::percent) + 
    xlab("proportion of cell") + labs(fill="prediction")

table(x$orig.ident,cut(x$pred,c(-Inf,-1,0,1,+Inf),c("---","-","+","+++")),x$fold,x$pred_groups,useNA="i") |>
  as.data.frame() |> 
  subset(Var3 == "tst") |> 
  #subset(Var2 %in% c("---","+++")) |> 
  ggplot() + facet_wrap(~Var4) +
  geom_col(aes(y=Var1,x=Freq,fill=Var2),position = "fill",color="black") + 
  geom_vline(xintercept=0.3,linetype=2) + 
  scale_x_continuous(label=scales::percent) + 
  xlab("proportion of cell") + labs(fill="prediction")










