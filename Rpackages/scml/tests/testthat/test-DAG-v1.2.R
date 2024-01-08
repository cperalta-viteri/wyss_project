

# Extract constraints:
# For each possible label, it lists pos/neg pairs involved in ancestor nodes
hierachical_cst <- function(h) {
  # Add ancestors list
  h <- mutate(h, target = row_number(),ancestors = local_members(order=+Inf,mode = "out"))

  # Ancestors table (self included)
  A <- as_tibble(h,"nodes") |> select(target,ancestors) |> unnest("ancestors")

  # Add ancestors childs (self excluded)
  Z <- A |>
    filter(target != ancestors) |>
    inner_join(as_tibble(h,"edges"),by=c("ancestors"="to"),relationship="many-to-many")

  # Check if ancestors childs are in the ancestors list
  G <- bind_rows(
    semi_join(Z,A,by=c("target","from"="ancestors")),
    anti_join(Z,A,by=c("target","from"="ancestors")),
    .id = "dir"
  ) |>
    mutate(dir=case_match(dir,"1" ~ +1,"2" ~ -1))

  # Prettify
  G |>
    relocate(target,ancestors,from,dir) |>
    arrange(target,ancestors,from,dir) |>
    rename(ancestor="ancestors")
}

G <- hierachical_cst(h)

m <- torch_sparse_coo_tensor(
  indices = torch_tensor(rbind(G$target,G$ancestor,G$from)),
  values = G$dir,
  size = c(length(h),length(h),length(h))
)
m$index_select(1L,2L)





#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# New tensor-based implementation
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#








model <- nn_sequential(
  nn_linear(4,2),nn_relu(),
  nn_linear(2,length(hlabs))
)

input <- model(data.matrix(iris[1:4]))
target <- as.integer(iris$Species)

nnf_margin_ranking_loss(2:11,1:10,rep(-1,10),margin = 1,reduction = "none")


M$to_dense()
M$matrix_power(3)
torch_bincount()

# Cost matrix ? = or determine number of label in each child













# Define the model which produces one prediction per node in the hierarchy
fm <- nn_sequential |>
  setup(
    loss = nnf_multi_margin_loss,
    optimizer = optim_adam,
    metrics = list(luz_metric_accuracy())
  ) |>
  set_hparams(nn_linear(4,2),nn_relu(),nn_linear(2,length(hlabs))) |>
  set_opt_hparams(lr=0.1) |>
  fit(list(data.matrix(iris[1:4]),as.integer(iris$Species)),epochs=30,accelerator = accelerator(cpu=TRUE))


p <- predict(fm,data.matrix(iris[1:4]))
table(iris$Species,as_array(p$max(2)[[2]]$to(device="cpu")))

# Define the loss








