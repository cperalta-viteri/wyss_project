#Sys.setenv(PYTORCH_ENABLE_MPS_FALLBACK=1)
#Sys.getenv("PYTORCH_ENABLE_MPS_FALLBACK")

library(torch)
library(luz)
source("notebooks/celltype_prediction/lib_cell_classifier.R")




#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#
#
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
trn_ds <- readRDS("notebooks/celltype_prediction/data/allen_subclass_trn.rds")
fn <- colnames(trn_ds$x)
lev <- levels(trn_ds$y) 
trn_ds <- tensor_dataset(
  x = torch_tensor(as.matrix(trn_ds$x)),
  y = torch_tensor(trn_ds$y,dtype = torch_long())
)  

model <- advanced_cell_classifier |>
  setup(metrics = list(acc=luz_metric_accuracy())) |>
  set_hparams(feature_names=fn,class_names=lev) |>
  set_opt_hparams(lr1=1e-5) |>
  fit(trn_ds, epochs = 100L,valid_data = trn_ds,dataloader_options=list(batch_size=256L,shuffle=TRUE),callbacks = list(
    luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"notebooks/celltype_prediction/models/model_subclass_best_valid_acc.luz")})(),
    luz_callback_model_checkpoint("notebooks/celltype_prediction/models/model_subclass_best_valid_acc.luz",monitor="valid_acc",mode="max",save_best_only=TRUE),
    luz_callback_model_checkpoint("notebooks/celltype_prediction/models/model_subclass_best_train_acc.luz",monitor="train_acc",mode="max",save_best_only=TRUE)
  ))
luz_load_model_weights(model,"notebooks/celltype_prediction/models/model_subclass_best_valid_acc.luz")
msk <- model$model$compute_prune_mask()


model_pruned <- advanced_cell_classifier |>
  setup(metrics = list(acc=luz_metric_accuracy())) |>
  set_hparams(feature_names=model$model$feature_names,class_names=model$model$class_names,prune_mask=msk) |>
  set_opt_hparams(lr1=1e-5) |>
  fit(trn_ds, epochs = 20L,valid_data = trn_ds,dataloader_options=list(batch_size=256L,shuffle=TRUE),callbacks = list(
    luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"notebooks/celltype_prediction/models/model_subclass_best_valid_acc.luz")})(),
    #luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"notebooks/celltype_prediction/models/model_subclass_pruned_best_train_acc.luz")})(),
    luz_callback_model_checkpoint("notebooks/celltype_prediction/models/model_subclass_pruned_best_train_acc.luz",monitor="train_acc",mode="max",save_best_only=TRUE),
    luz_callback_model_checkpoint("notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.luz",monitor="valid_acc",mode="max",save_best_only=TRUE)
  ))

luz_load_model_weights(model_pruned,"notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.luz")
model_pruned$model(trn_ds[1:5]$x)
torch_save(model_pruned$model,"notebooks/celltype_prediction/models/model_subclass_pruned_best_valid_acc.torch")








