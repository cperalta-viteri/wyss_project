library(torch)
library(luz)

ds <- tensor_dataset(
  x = torch_tensor(t(as.matrix(assay(ML,"data")))),
  y = torch_tensor(ML$infected+1L,dtype = torch_long())
)
dl <- dataloader(ds,batch_size = 64)
net <- advanced_cell_classifier(rownames(ML),c("uninfected","infected"))

advanced_cell_classifier |>
  setup(
    #loss = nn_multi_margin_loss(),
    loss = nn_multi_margin_loss(),
    optimizer = optim_adam
  ) |>
  set_opt_hparams(lr=0.0001) |>
  set_hparams(feature_names=rownames(ML),c("uninfected","infected")) |>
  fit(dl, epochs = 30)
