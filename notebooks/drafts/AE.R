

library(torch)
library(luz)


torch_set_default_dtype(torch_float32())

x <- readRDS("notebooks/celltype_prediction/data/allen_subclass_trn.rds")
y <- x$y;x <- as.matrix(x$x)
batch <- substr(rownames(x),3,3)=="-"

# x <- data.matrix(iris[1:4])
# y <- iris$Species
ml_ds <- tensor_dataset(
  x = torch_tensor(x),
  y = torch_tensor(cbind(sign(as.integer(batch)-0.5),rowMeans(x)))
)


ae <- nn_module(
  "ae",
  
  initialize = function(feature_names,latent_dims=2L) {
    self$input_dropout <- nn_dropout()
    
    self$encoder <- nn_sequential(
      nn_linear(in_features = length(feature_names), out_features = 512L),
      nn_relu(),
      nn_linear(in_features = 512L, out_features = 512L),
      nn_relu(),
      nn_linear(in_features = 512L, out_features = latent_dims),
    )
    
    self$decoder <- nn_sequential(
      nn_linear(in_features = latent_dims+2L, out_features = 512L),
      nn_relu(),
      nn_linear(in_features = 512L, out_features = 512L),
      nn_relu(),
      nn_linear(in_features = 512L, out_features = length(feature_names))
    )
  },
  
  forward = function(x) {
    z <- x |>
      self$input_dropout() |>
      self$encoder()
    if (ctx$training) {
      z <- distr_normal(z,1)$rsample()
    }
    list(x=x,z=z)
  },
  
  loss = function(y_pred,y_true) {
    z2 <- torch_cat(list(y_pred$z,y_true),-1)
    l <- nnf_mse_loss(y_pred$x,self$decoder(z2)) + y_pred$z$norm()/100
    l
  }
)


model <- ae |>
  setup(optimizer = optim_adam) |>
  set_opt_hparams(lr=1e-4) |>
  set_hparams(colnames(x),latent_dims=3L) |>
  fit(ml_ds, epochs=300L, dataloader_options=list(batch_size=256L,shuffle=TRUE),callbacks = list(
    luz_callback(on_fit_begin=function(){luz_load_model_weights(ctx,"ae_last.luz")})(),
    luz_callback_model_checkpoint("ae_best.luz",monitor="train_loss",save_best_only=TRUE),
    luz_callback_model_checkpoint("ae_last.luz",monitor="train_loss")
  ))
luz_load_model_weights(model,"ae_last.luz")

library(ggplot2)
library(dplyr)
u <- model$model$encoder(x) |> 
  as_array()

u |> 
  as.data.frame() |>
  mutate(lev3 = y) |>
  mutate(lev2 = sub("(((.*) --- .*) --- (.*))","\\2",lev3)) |>
  mutate(lev1 = sub("(((.*) --- .*) --- (.*))","\\3",lev3)) |>
  ggplot() + #facet_wrap(~lev1) +
    geom_point(aes(x=V1,y=V2,color=lev1)) + 
    coord_equal()



rgl::clear3d()
rgl::plot3d(u,col=as.integer(factor(sub("(((.*) --- .*) --- (.*))","\\3",y))))
rgl::rglwidget()


rgl::clear3d()
rgl::plot3d(u,col=as.integer(batch)+1)
rgl::rglwidget()


rgl::clear3d()
rgl::plot3d(u,col=as.integer(batch)+1)
rgl::rglwidget()


col <- rowMeans(x)/max(rowMeans(x))
col <- rgb(1-col,1-col,1)
rgl::clear3d()
rgl::plot3d(u,col=col)
rgl::rglwidget()






