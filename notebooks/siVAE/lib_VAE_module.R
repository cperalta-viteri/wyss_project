
if(!require(pacman, quietly = TRUE))install.packages("pacman")


pacman::p_load(torch, luz)


#GGally, ggfortify, fastDummies, purrr,dplyr


# --- data prep ----
demo_geo <- demo_geo[, sample(colnames(demo_geo))]
# to reduce n of genes
demo_geo <- demo_geo[((counts(demo_geo)>0) |> rowMeans()) > 0.01]

ml_ds <- tensor_dataset(
  x = torch_tensor(t(assay(demo_geo, "logcounts"))),
  y = torch_tensor(t(assay(demo_geo, "logcounts")))
)


# ---------------------
#       Model
# ---------------------

cell_encoder <- nn_module(
  "cell-w_siVAE",
  #---- initializer -----
  # will contain mu and sigma attributes for the model
  # contains encoder-decoder layers
  initialize = function(gene_names, embedding_size=2L){
    self$gene_names <- gene_names
    self$log_scale <- nn_parameter(torch_tensor(0.0))
    # --- layers ---
    # self$layer1 <- nn_linear(in_features = 52739L, out_features = 50L)
    self$l_encoder <- nn_sequential(
      nn_linear(in_features = length(gene_names), out_features = 1024L),
      nn_relu(),
      nn_linear(in_features = 1024L, out_features = 128L),
      nn_relu()
    )
    self$layer_mu <-    nn_linear(in_features = 128L, out_features = embedding_size)
    self$layer_sigma <- nn_linear(in_features = 128L, out_features = embedding_size)
    #self$l_embedding <- nn_linear(in_features = 100L, out_features = embedding_size)
    
    # self$layer2 <- nn_linear(in_features = 2L, out_features = 50L)
    # self$layer3 <- nn_linear(in_features = 50L, out_features = 52739L)
    self$l_decoder <- nn_sequential(
      nn_linear(in_features = embedding_size, out_features = 128L),
      nn_relu(),
      nn_linear(in_features = 128L, out_features = 1024L),
      nn_relu(),
      nn_linear(in_features = 1024L, out_features = length(gene_names))
    )
  },
  
  #---- encoder ----
  # computes mu and sigma layers
  encode = function(x){
    # computes embedding's mu & std values
    list(mu=self$l_encoder(x) |> self$layer_mu(),
         std=self$l_encoder(x) |> self$layer_sigma() |> torch_exp() |> torch_sqrt()) # compute std here
  },
  
  #---- decoder ----
  # generative part
  decode = function(z){
    self$l_decoder(z)
  },
  
  forward = function(x){
    enc <- x |> 
      self$encode(x)
    
    # Define the approximated posterior probability q(z)
    q_z <- distr_normal(enc$mu, enc$std) #attribute
    
    # reparametrization trick
    ####z=self$mu+(self$std*torch_rand_like(self$std))
    q_z$rsample() |>
      self$decode()
    
  },
  
  set_optimizers= function(){
    optim_adam(self$parameters, weight_decay = 0.001)
  },
  
  step = function(){
    opt <- ctx$optimizers[[1]]
    
    x <- ctx$input
    enc <- self$encode(x)
    # q(z)
    q <- distr_normal(enc$mu, enc$std)
    # p(z)
    z <- q$rsample()
    
    ctx$pred <- z |> self$decode()
    #browser()
    # ---- loss computation ----
    loss <- elbo_loss(pred=ctx$pred, target = ctx$target, enc$mu, enc$std, q, z, self$log_scale)
    ctx$loss <- list(elbo=loss)
    
    
    if (ctx$training){
      opt$zero_grad()
      ctx$loss[[1]]$backward()
      opt$step()
    }
    ctx$loss[[1]] <- loss$detach()
  }
)



cell_wise_model <- cell_encoder |>
  setup() |>
  set_hparams(gene_names=rownames(demo_geo)) |>
  fit(ml_ds, epochs=10L, dataloader_options=list(batch_size=30L))
# 3h of training for siVAE v2




#p=distr_normal(torch_zeros_like(0), torch_ones_like(1))
#q=distr_normal(torch_zeros_like(0), torch_ones_like(2))



gaussian_likelihood <- function(preds, logscale, x){
  scale <- torch_exp(logscale)
  mean <- preds
  dist <- distr_normal(mean, scale)
  
  # computing probability p(x|z)
  log_pxz <- dist$log_prob(x)
  log_pxz <- log_pxz$sum()
}

kl_divergence <- function(q, z, mu, std){
  # 1. Getting the missing probability distribution
  #p_z<-distr_normal(torch_zeros_like(mu), torch_ones_like(std))
  p_z<-distr_normal(mu, std)
  
  # 2. Computing probabilities
  log_qzx <- q$log_prob(z)
  log_pz <- p_z$log_prob(z)
  kl_d <- log_qzx - log_pz
  kl_d <- kl_d$sum(-1)
  kl_d
}

# enc to mu
elbo_loss <- function(pred, target, mu, std, q, z, log_scale){
  recon_loss <- gaussian_likelihood(pred, log_scale, target)
  kl <- kl_divergence(q, z, mu, std)
  elbo <- kl - recon_loss
  elbo <- elbo$mean()
  elbo
}


#loss = 0.5 * ((logvar2 - logvar_1) - 1 + var_1 / var_2 + (mu_2 - mu_1)**2 / var_2)


