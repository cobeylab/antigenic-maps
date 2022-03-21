## Infer maps for each set of inputs

rm(list = ls())
source('../../../code.R')
source('../../../multi_epitope_code.R')
library(foreach)
library(doParallel)
library(rstan)
library(tidyverse)
if(!dir.exists('diagnostics')) dir.create('diagnostics')



fit_stan_MDS <- function(
  mod,
  observed_titers, # N-vector of observed titers
  smax, # N-vector of smax
  serum_id, # N-vector of serum id
  antigen_id, # N- vector of antigen id
  n_antigens, # integer
  n_sera, # integer
  n_dim, # integer
  N_test_set,
  smax_test_set,
  serum_id_test_set,
  antigen_id_test_set,
  observed_titers_test_set,
  chains = 3, # Number of MCMC chains to run
  cores = parallel::detectCores(logical = F), # For the cluster
  niter = 5000,
  antigen_coords,
  diagnostic_file = NULL,
  ...
) {
  
  
  model <- stan_model(mod)
  #(print(model))
  
  model_input_data <- list(
    N = length(observed_titers),
    n_sera = n_sera,
    n_antigens = n_antigens,
    n_dim = n_dim,
    observed_titers = observed_titers,
    smax = smax,
    serum_id = serum_id,
    antigen_id = antigen_id,
    N_test_set = N_test_set,
    smax_test_set = smax_test_set,
    serum_id_test_set = serum_id_test_set,
    antigen_id_test_set = antigen_id_test_set,
    observed_titers_test_set = observed_titers_test_set
  )
  
  initfun <- function(){
    list(sigma = 1,
         ag2_c1 = runif(1, 0, 10),
         strain_coords = matrix(runif((n_antigens-2)*n_dim, -10, 10), n_antigens-2, n_dim),
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
  
  if(all(summary(initial_fit)$summary[,'Rhat'] <= 1.03)){
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
           ag2_c1 = best_chain_summary['ag2_c1', '50%'],
           antigen_coords = matrix(best_chain_summary[is.antigen.coord,'50%'], 
                                   nrow = n_antigens-2, 
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
  
  if(! all(summary(refit)$summary[,'Rhat'] <= 1.03)){
    cat(sprintf('Re-doing refit to achieve Rhat < 1.03'))
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
                    titer_map_train,
                    titer_map_test,
                    n_chains = 4,
                    n_iter= 7000,
                    idstring){
  
  nsera = unique(titer_map_train$serum) %>% length()
  nantigen = unique(titer_map_train$antigen) %>% length()
  cat(sprintf('in R: n_antigens is %i; n_sera is %i \n', nantigen, nsera))
  
  stopifnot(all(titer_map_test$serum %in% titer_map_train$serum))
  stopifnot(all(titer_map_test$antigen %in% titer_map_train$antigen))
  
  fit_stan_MDS(mod = '../../../stan/MDS_predict_titer_sparse.stan',
               observed_titers = titer_map_train$logtiter, # N-vector of observed titers
               smax = (titer_map_train$serum_potency+titer_map_train$antigen_avidity)/2, # N-vector of smax
               serum_id = titer_map_train$serum, # N-vector of serum id
               antigen_id = titer_map_train$antigen, # N- vector of antigen id
               n_antigens = nantigen, # integer
               n_sera = nsera, # integer
               n_dim = n_dim, # integer
               N_test_set = nrow(titer_map_test),
               smax_test_set = (titer_map_test$serum_potency + titer_map_test$antigen_avidity)/2,
               serum_id_test_set = titer_map_test$serum,
               antigen_id_test_set = titer_map_test$antigen,
               observed_titers_test_set = titer_map_test$logtiter,
               chains = 3, # Number of MCMC chains to run
               cores = parallel::detectCores(logical = F), # For the cluster
               niter = 5000,
               diagnostic_file = NULL
  )
  
}



## Fit the even immunodominance maps for 1:8 dimensions
set.seed(11)
even_inputs <- read_rds('../even_immunodominance_inputs.rds')
split_even_map <- test_train_split(25, 5, even_inputs$titer_map)
stopifnot(1:5 %in% split_even_map$train$antigen)
stopifnot(1:5 %in% split_even_map$train$serum)
even_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map_train = split_even_map$train, 
          titer_map_test = split_even_map$test,
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'even')
}
names(even_fit_list) = paste0(1:8, 'D')
write_rds(even_fit_list, 'even_fit_list.rds')
write_rds(split_even_map, 'even_inputs_test_train_split.rds')



set.seed(11)
skewed_inputs <- read_rds('../skewed_immunodominance_inputs.rds')
split_skewed_map <- test_train_split(25, 5, skewed_inputs$titer_map)
skewed_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map_train = split_skewed_map$train,
          titer_map_test = split_skewed_map$test,
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'skewed')
}
names(skewed_fit_list) = paste0(1:8, 'D')
write_rds(skewed_fit_list, 'skewed_fit_list.rds')
write_rds(split_skewed_map, 'skewed_inputs_test_train_split.rds')


set.seed(11)
E1_inputs <- read_rds('../E1_complete_dominance_inputs.rds')
split_E1_map <- test_train_split(25, 5, E1_inputs$titer_map)
E1_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map_train = split_E1_map$train, 
          titer_map_test = split_E1_map$test,
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'E1')
}
names(E1_fit_list) = paste0(1:8, 'D')
write_rds(E1_fit_list, 'E1_fit_list.rds')
write_rds(split_E1_map, 'E1_inputs_test_train_split.rds')


set.seed(11)
E2_inputs <- read_rds('../E2_complete_dominance_inputs.rds')
split_E2_map <- test_train_split(25, 5, E2_inputs$titer_map)
E2_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map_train =  split_E2_map$train, 
          titer_map_test = split_E2_map$test,
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'E2')
}
names(E2_fit_list) = paste0(1:8, 'D')
write_rds(E2_fit_list, 'E2_fit_list.rds')
write_rds(split_E2_map, 'E2_inputs_test_train_split.rds')


set.seed(11)
E3_inputs <- read_rds('../E3_complete_dominance_inputs.rds')
split_E3_map <- test_train_split(25, 5, E3_inputs$titer_map)
E3_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map_train = split_E3_map$train, 
          titer_map_test = split_E3_map$test,
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'E3')
}
names(E3_fit_list) = paste0(1:8, 'D')
write_rds(E3_fit_list, 'E3_fit_list.rds')
write_rds(split_E3_map, 'E3_inputs_test_train_split.rds')
