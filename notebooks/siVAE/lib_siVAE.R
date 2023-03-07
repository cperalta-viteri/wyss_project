
if(!require(pacman, quietly = TRUE))install.packages("pacman")


pacman::p_load(torch, luz)

#GGally, ggfortify, fastDummies, purrr,dplyr


# --- data prep ----
ml_ds <- tensor_dataset(
  x = torch_tensor(t(assay(demo_geo, "logcounts"))),
  y = torch_tensor(t(assay(demo_geo, "logcounts")))
)




cell_encoder <- nn_module(
  "cell-w_siVAE",
  #---- initializer -----
  # will contain mu and sigma attributes for the model
  # contains encoder-decoder layers
  initialize = function(){
    #self$q_z <- NULL
    #self$z <- NULL
    #self$enc <- NULL
    self$log_scale <- nn_parameter(torch_tensor(0.0))
    self$layer1 <- nn_linear(in_features = 52739L, out_features = 50L)
    self$layer_mu <-    nn_linear(in_features = 50L, out_features = 2L)
    self$layer_sigma <- nn_linear(in_features = 50L, out_features = 2L)
    self$layer2 <- nn_linear(in_features = 2L, out_features = 50L)
    self$layer3 <- nn_linear(in_features = 50L, out_features = 52739L)
  },
  
  #---- encoder ----
  # computes mu and sigma layers
  encode = function(x){
    # applies relu function on first layer
    #browser()
    r1 <- nnf_relu(self$layer1(x))
    # computes mu & std values, updating attributes
    list(mu=self$layer_mu(r1),log_var=self$layer_sigma(r1))
  },
  
  #---- decoder ----
  # generation part, layer 3
  decode = function(z){
    # applies relu function from latent space
    r2 <- nnf_relu(self$layer2(z))
    # updates layer 3, generating original data
    self$layer3(r2)
  },
  
  forward = function(x){
    enc <- x |> 
      self$encode(x)
    
    # Define the approximated posterior probability q(z)
    q_z <- distr_normal(enc$mu, enc$log_var$exp()$pow(0.5)) #attribute
    
    # reparametrization trick
    ####z=self$mu+(self$std*torch_rand_like(self$std))
    q_z$rsample() |>
      self$decode()
    
  },
  
  # gaussian_likelihood = function(preds, logscale, x){
  #   scale <- torch_exp(logscale)
  #   mean <- preds
  #   dist <- distr_normal(mean, scale)
  #   
  #   # computing probability p(x|z)
  #   log_pxz <- dist$log_prob(x)
  #   log_pxz <- log_pxz$sum()
  # },
  # 
  # kl_divergence = function(q, z, mu, log_var){
  #   # 1. Getting the missing probability distribution
  #   p_z<-distr_normal(torch_zeros_like(mu), torch_ones_like(log_var))
  #   #q=distr_normal(self$mu, self$std)
  # 
  #   # 2. Computing probabilities
  #   log_qzx <- q$log_prob(z)
  #   log_pz <- p_z$log_prob(z)
  #   kl_d <- log_qzx - log_pz
  #   kl_d <- kl_d$sum(-1)
  #   kl_d
  # },
  # 
  # elbo_loss = function(pred, target, enc, std, q, z){
  #   recon_loss <- self$gaussian_likelihood(pred, self$log_scale, target)
  #   kl <- self$kl_divergence(q, z, enc$mu, std)
  #   elbo <- kl-recon_loss
  #   elbo <- elbo$mean()
  #   elbo
  # },
  
  set_optimizers= function(){
    optim_adam(self$parameters)
  },
  
  step = function(){
    #ctx$loss <- list()
    opt <- ctx$optimizers[[1]]
    
    x <- ctx$input
    enc <- self$encode(x)
    std <- enc$log_var$exp()$pow(0.5)
    q <- distr_normal(enc$mu, std)
    z <- q$rsample()
    
    ctx$pred <- z |> self$decode()
    
    # ---- loss computation ----
    # recon_loss <- self$gaussian_likelihood(x_hat, self$log_scale, x$y)
    # kl <- self$kl_divergence(q, z, enc$mu, std)
    # elbo <- kl - recon_loss
    # loss <- elbo$mean()
    loss <- elbo_loss(pred=ctx$pred, target = ctx$target, enc, std, q, z, self$log_scale)
    ctx$loss <- list(elbo=loss)
    ###
    #ctx$loss <- ctx$loss_fn(pred=ctx$pred, target = ctx$target, enc, std, q, z, self$log_scale)
    #ctx$loss <- loss #function(pred, target, enc, std, q, z, log_scale)
    
    if (ctx$training){
      opt$zero_grad()
      ctx$loss[[1]]$backward()
      opt$step()
    }
    ctx$loss[[1]] <- loss$detach() #ctx$loss does not exist
    #browser()
  }
)



cell_wise_model <- cell_encoder |>
  setup() |>
  fit(ml_ds, epochs=100L, dataloader_options=list(batch_size=30L))





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

kl_divergence <- function(q, z, mu, log_var){
  # 1. Getting the missing probability distribution
  p_z<-distr_normal(torch_zeros_like(mu), torch_ones_like(log_var))
  #q=distr_normal(self$mu, self$std)
  
  # 2. Computing probabilities
  log_qzx <- q$log_prob(z)
  log_pz <- p_z$log_prob(z)
  kl_d <- log_qzx - log_pz
  kl_d <- kl_d$sum(-1)
  kl_d
}

elbo_loss <- function(pred, target, enc, std, q, z, log_scale){
  recon_loss <- gaussian_likelihood(pred, log_scale, target)
  kl <- kl_divergence(q, z, enc$mu, std)
  elbo <- kl-recon_loss
  elbo <- elbo$mean()
  elbo
}


