
source("notebooks/celltype_prediction/lib_cell_classifier.R")

trn_ds <- delayed_dataset(
  x = data.matrix(iris[1:4]),
  y = as.integer(iris$Species)
)


C <- 1 - torch_eye(nlevels(iris$Species))
fm <- cell_scorer |>
  setup(
    loss = nn_ordinal_regression_loss(C=C),
    optimizer = optim_adam,
    metrics = list()
  ) |>
  set_opt_hparams(lr=0.001,weight_decay=0.2) |>
  set_hparams(feature_names=names(iris)[1:4],class_names="age") |>
  fit(trn_ds,
      valid_data = trn_ds,
      epochs = 300L,
      dataloader_options = list(batch_size=32L,drop_last=FALSE),
      accelerator = accelerator(cpu = TRUE),
      callbacks = list(
        luz_callback_lr_scheduler(torch::lr_step,step_size=100L,gamma=0.1)
        #luz_callback_model_checkpoint(monitor="train_loss","~/Downloads/last_full_model.luz")
        #luz_callback_auto_resume(path="~/Downloads/autoresume.pt")
      )
  )

#luz_save(fm,"~/Downloads/final_model.luz")
boxplot(as_array(predict(fm,trn_ds)$to(device="cpu"))~iris$Species)
print(plot(fm))

p <- predict(fm,trn_ds,dataloader_options = list(batch_size = 1000L),accelerator = accelerator(cpu = TRUE))

nn_ordinal_regression_loss(C)(p,torch_tensor(as.integer(iris$Species),torch_long()))

evaluate(fm,data = dataloader(trn_ds,batch_size = 1000L),accelerator = accelerator(cpu = TRUE))



