## Infer maps for each set of inputs

rm(list = ls())
source('../../../code.R')
source('../../../multi_epitope_code.R')
library(foreach)
library(doParallel)
dir.create('diagnostics')


fit_unconstrined_stan_MDS <- function(
  mod = '../../../Bayesian_stan/MDS_unconstrained.stan',
  observed_distances, # n_antigen x n_antibody matrix of distances
  n_antigens, # integer
  n_sera, # integer
  n_dim, # integer
  chains = 4, # Number of MCMC chains to run
  cores = 4, # For the cluster
  niter = 5000,
  antigen_coords,
  diagnostic_file = NULL,
  ...
) {
  library(rstan)
  
  
  stopifnot(nrow(observed_distances) == n_antigens)
  stopifnot(ncol(observed_distances) == n_sera)
  
  model <- stan_model(mod)
  #(print(model))
  
  model_input_data <- list(
    n_antigens = n_antigens,
    n_sera = n_sera,
    n_dim = n_dim,
    observed_distances = observed_distances
  )
  
  initfun <- function(){
    list(sigma = 1,
         strain_coords = matrix(runif(n_antigens*n_dim, -10, 10), n_antigens, n_dim),
         serum_coords =  matrix(runif(n_sera*n_dim, -10, 10), n_sera, n_dim)
    )
  }
  
  initial_fit <- sampling(model, 
                          data = model_input_data, 
                          chains = chains, 
                          init = initfun,
                          iter = niter,
                          cores = min(6, chains),
                          control = list(adapt_delta = 0.89,
                                         max_treedepth = 14),
                          diagnostic_file = diagnostic_file
  )
  
  if(all(summary(initial_fit)$summary[,'Rhat'] <= 1.02)){
    cat(sprintf('Returning initial fit'))
    return(initial_fit)
  }
  
  cat(print('Initial fit complete.\nRe-running using initial values from the best chain.\n'))
  
  initialize_with_best_chain <- function(initial_fit, nchains){
    ## Extract the summary of the best chain in terms of mean log posterior
    best_chain_index = which.max(rstan::extract(initial_fit, permute = F)[,,'lp__'] %>% colMeans())
    best_chain_summary = rstan::summary(initial_fit)$c_summary[,,best_chain_index]
    best_chain = rstan::extract(initial_fit, permute = F)[,best_chain_index,]
    
    is.antigen.coord = sapply(dimnames(best_chain_summary)$parameter, FUN = function(xx) grepl('antigen_coords', xx))
    is.serum.coord = sapply(dimnames(best_chain_summary)$parameter, FUN = function(xx) grepl('serum_coords', xx))
    
    
    get_one_list <- function(){
      ## Output an initial list using the median values from the best chain
      list(sigma = best_chain_summary['sigma', '50%'],
           antigen_coords = matrix(best_chain_summary[is.antigen.coord,'50%'], 
                                   nrow = n_antigens, 
                                   ncol = n_dim,
                                   byrow = T),
           serum_coords = matrix(best_chain_summary[is.serum.coord,'50%'], 
                                 nrow = n_sera, 
                                 ncol = n_dim,
                                 byrow = T))
    }
    
    lapply(1:nchains, function(xx){get_one_list()})
  }
  
  inits <- initialize_with_best_chain(initial_fit, 1)
  
  refit <- sampling(
    model, model_input_data, 
    chains = 1, 
    init = inits,
    iter = niter,
    diagnostic_file = diagnostic_file)
  
  if(! all(summary(refit)$summary[,'Rhat'] <= 1.02)){
    cat(sprintf('Re-doing refit to achieve Rhat < 1.1'))
    refit <- sampling(
      model, model_input_data, 
      chains = 1, 
      init = inits,
      iter = niter,
      diagnostic_file = diagnostic_file)
  }
  
  return(refit)
}




one_fit <- function(n_dim,
                    titer_map,
                    n_chains = 4,
                    n_iter= 5000,
                    idstring){
  
  nsera = unique(titer_map$serum) %>% length()
  nantigen = unique(titer_map$antigen) %>% length()
  
  cat(sprintf('Fitting %sD model', n_dim))
  fit_unconstrined_stan_MDS(mod = '../../../Bayesian_stan/MDS_unconstrained.stan',
               observed_distances = format_stan_inputs(titer_map),
               n_antigens = nsera,
               n_sera = nantigen,
               n_dim = n_dim,
               chains = n_chains,
               niter = n_iter,
               diagnostic_file = paste0('diagnostics/', idstring, '_', n_dim, 'D'))
  
}



## Fit the even immunodominance maps for 1:8 dimensions
even_inputs <- read_rds('../even_immunodominance_inputs.rds')
even_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map = even_inputs$titer_map, 
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'penalized_even')
}
names(even_fit_list) = paste0(1:8, 'D')



# ## Fit the skewed immunodominance maps for 1:8 dimensions
# skewed_inputs <- read_rds('../skewed_immunodominance_inputs.rds')
# skewed_fit_list <- foreach(this_ndim = 1:8) %do% {
#   one_fit(titer_map = skewed_inputs$titer_map, 
#           n_dim = this_ndim, 
#           n_chains = 4, 
#           n_iter = 5000,
#           idstring = 'penalized_skewed')
# }
# names(skewed_fit_list) = paste0(1:8, 'D')
# 
# 
# 
# 
## Fit the E1 complete immunodominance maps for 1:8 dimensions
E1_inputs <- read_rds('../E1_complete_dominance_inputs.rds')
E1_fit_list <- foreach(this_ndim = 2) %do% {
  one_fit(titer_map = E1_inputs$titer_map,
          n_dim = this_ndim,
          n_chains = 4,
          n_iter = 5000,
          idstring = 'penalized_E1')
}
names(E1_fit_list) = paste0(1:8, 'D')
# 
# 
# 
# 
# ## Fit the E2 complete immunodominance maps for 1:8 dimensions
# E2_inputs <- read_rds('../E2_complete_dominance_inputs.rds')
# E2_fit_list <- foreach(this_ndim = 1:8) %do% {
#   one_fit(titer_map = E2_inputs$titer_map, 
#           n_dim = this_ndim, 
#           n_chains = 4, 
#           n_iter = 5000,
#           idstring = 'penalized_E2')
# }
# names(E2_fit_list) = paste0(1:8, 'D')
# 
# 
# 
# ## Fit the E3 complete immunodominance maps for 1:8 dimensions
# E3_inputs <- read_rds('../E3_complete_dominance_inputs.rds')
# E3_fit_list <- foreach(this_ndim = 1:8) %do% {
#   one_fit(titer_map = E3_inputs$titer_map, 
#           n_dim = this_ndim, 
#           n_chains = 4, 
#           n_iter = 5000,
#           idstring = 'penalized_E3')
# }
# names(E3_fit_list) = paste0(1:8, 'D')



write_rds(even_fit_list, 'even_fit_list.rds')
# write_rds(skewed_fit_list, 'skewed_fit_list.rds')
write_rds(E1_fit_list, 'E1_fit_list.rds')
# write_rds(E2_fit_list, 'E2_fit_list.rds')
# write_rds(E3_fit_list, 'E3_fit_list.rds')