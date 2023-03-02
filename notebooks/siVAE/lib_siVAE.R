if(!require(pacman, quietly = TRUE))install.packages("pacman")
pacman::p_load(torch, luz)

#GGally, ggfortify, fastDummies, purrr,
dplyr

cell_encoder <- nn_module(
  #---- initializer -----
  # will contain mu and sigma attributes for the model
  # contains encoder-decoder layers
  initialize = function(){
    self$mu <- NULL
    self$log_var <- NULL
    self$std <- NULL
    self$q_z <- NULL
    self$z <- NULL
    self$x_hat <- NULL
    self$layer1 <- nn_linear(in_features = 52739L, out_features = 50L)
    self$layer_mu <-    nn_linear(in_features = 50L, out_features = 2L)
    self$layer_sigma <- nn_linear(in_features = 50L, out_features = 2L)
    self$layer2 <- nn_linear(in_features = 2L, out_features = 50L)
    self$layer3 <- nn_linear(in_features = 50L, out_features = 52739L)
  },
  encode = function(x){
    # applies relu function on first layer
    #browser()
    r1 = nnf_relu(self$layer1(x))
    # computes mu & std values, updating attributes
    self$mu <- self$layer_mu(r1); self$log_var <- self$layer_sigma(r1)
  },
  decode = function(z){
    # applies relu function from latent space
    r2 <- nnf_relu(self$layer2(z))
    # updates layer 3, generating original data
    self$layer3(r2)
  },
  
  forward = function(x){
    self$encode(x)
    
    ex<-self$log_var$exp()
    self$std<-ex$pow(0.5)
    
    # Define the approximated posterior probability q(z)
    self$q_z=distr_normal(self$mu, self$std)
    
    # reparametrization trick
    ####z=self$mu+(self$std*torch_rand_like(self$std))
    self$z<-self$q_z$rsample()
    
    self$decode(self$z)
    
  },
  
  gaussian_likelihood = function(x_hat, logscale, x){
    scale = torch_exp(logscale)
    mean = x_hat
    dist = distr_normal(mean, scale)
    
    # measure probability p(x|z)
    log_pxz = dist$log_prob(x)
    return(log_pxz = log_pxz$sum())
  },
  
  kl_divergence = function(){
    # 1. Getting the missing probability distribution
    p_z=distr_normal(torch_zeros_like(self$mu$loc), torch_ones_like(self$std$scale))
    #q=distr_normal(self$mu, self$std)
    
    # 2. Computing probabilities
    log_qzx = self$q_z$log_prob(self$z)
    log_pz = p_z$log_prob(self$z)
    kl = log_qzx - log_pz
    kl = kl$sum(-1)
    return(kl)
  },
  
  loss=function(pred, target){
    recon_loss = self$gaussian_likelihood(pred, nn_parameter(torch_tensor(0.0)), target)
    kl = self$kl_divergence()
    
    elbo=kl-recon_loss
    elbo=elbo$mean()
  },
  
  set_optimizers= function(){
    optim_adam(cell_encoder()$parameters)
  }
    
)


#p=distr_normal(torch_zeros_like(0), torch_ones_like(1))
#q=distr_normal(torch_zeros_like(0), torch_ones_like(2))


ml_ds <- tensor_dataset(
  x = torch_tensor(t(assay(demo_geo, "logcounts"))),
  y = torch_tensor(t(assay(demo_geo, "logcounts")))
)



cell_wise_model <- cell_encoder |>
  setup() |>
  luz::fit(ml_ds, epochs=5L, valid_data=ml_ds, dataloader_options=list(batch_size=30L))








