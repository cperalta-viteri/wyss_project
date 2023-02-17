#Sys.setenv(PYTORCH_ENABLE_MPS_FALLBACK=1)
#Sys.getenv("PYTORCH_ENABLE_MPS_FALLBACK")

library(torch)
library(luz)
source("notebooks/celltype_prediction/lib_cell_classifier.R")

ds <- readRDS("notebooks/celltype_prediction/dataset.rds")
fn <- colnames(ds$x)
ml_ds <- tensor_dataset(
  x = torch_tensor(ds$x),
  y = torch_tensor(ds$y,dtype = torch_long())
)

model <- advanced_cell_classifier |>
  setup(metrics = list(acc=luz_metric_accuracy())) |>
  set_hparams(feature_names=fn,class_names=levels(ds$y)) |>
  set_opt_hparams(lr1=1e-4) |>
  fit(ml_ds, epochs = 100L,valid_data = ml_ds,dataloader_options=list(batch_size=128L),callbacks = list(
    #luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"model_best_valid_acc.luz")})(),
    luz_callback_model_checkpoint("model_best_valid_acc.luz",monitor="valid_acc",mode="max",save_best_only=TRUE),
    luz_callback_model_checkpoint("model_best_train_acc.luz",monitor="train_acc",mode="max",save_best_only=TRUE)
  ))
luz_load_model_weights(model,"model_best_train_acc.luz")
msk <- model$model$compute_prune_mask()

model_pruned <- advanced_cell_classifier |>
  setup(metrics = list(acc=luz_metric_accuracy())) |>
  set_hparams(feature_names=model$model$feature_names,class_names=model$model$class_names,prune_mask=msk) |>
  set_opt_hparams(lr1=1e-4) |>
  fit(ml_ds, epochs = 50L,valid_data = ml_ds,dataloader_options=list(batch_size=128L),callbacks = list(
    #luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"model_pruned_best_train_acc.luz")})(),
    luz_callback_model_checkpoint("model_pruned_best_train_acc.luz",monitor="train_acc",mode="max",save_best_only=TRUE),
    luz_callback_model_checkpoint("model_pruned_best_valid_acc.luz",monitor="valid_acc",mode="max",save_best_only=TRUE)
  ))

luz_load_model_weights(model_pruned,"model_pruned_best_train_acc.luz")
torch_save(model_pruned$model,"model_pruned_best_train_acc.torch")








