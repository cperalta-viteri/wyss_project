library(torch)
library(luz)



ds <- readRDS("~/Downloads/ml_dataset.rds")
fn <- colnames(ds$x)
set.seed(12345)
fold <- sample(c("trn","val"),size = length(ds$y),replace=TRUE,prob = c(0.75,0.25))

trn_ds <- tensor_dataset(
  x = torch_tensor(ds$x)[fold=="trn",],
  y = torch_tensor(ds$y+1L,dtype = torch_long())[fold=="trn"]
)
trn_dl <- dataloader(trn_ds,batch_size = 64L)

val_ds <- tensor_dataset(
  x = torch_tensor(ds$x)[fold=="val",],
  y = torch_tensor(ds$y+1L,dtype = torch_long())[fold=="val"]
)
val_dl <- dataloader(val_ds,batch_size = 64L)


model <- advanced_cell_classifier |>
  setup(
    optimizer = optim_adam
  ) |>
  set_opt_hparams(lr=0.0001) |>
  set_hparams(feature_names=fn,class_names=c("uninfected","infected")) |>
  fit(trn_dl, epochs = 30L,valid_data = val_dl)

msk <- model$model$compute_prune_mask()

model_pruned <- advanced_cell_classifier |>
  setup(
    optimizer = optim_adam
  ) |>
  set_opt_hparams(lr=0.0001) |>
  set_hparams(feature_names=fn,class_names=c("uninfected","infected"),prune_mask=msk) |>
  fit(trn_dl, epochs = 30L,valid_data = val_dl)


y_true <- val_ds$tensors[[2]] |> as_array()  
y_pred <- model_pruned$model(val_ds$tensors[[1]])$argmax(2) |> as_array()
table(y_true,y_pred)

