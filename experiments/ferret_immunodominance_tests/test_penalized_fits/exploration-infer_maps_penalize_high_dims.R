## Infer maps for each set of inputs

rm(list = ls())
source('../../code.R')
source('../../multi_epitope_code.R')
library(foreach)
library(doParallel)
dir.create('diagnostics')


one_fit <- function(n_dim,
                    titer_map,
                    n_chains = 4,
                    n_iter= 5000,
                    idstring){
  
  nsera = unique(titer_map$serum) %>% length()
  nantigen = unique(titer_map$antigen) %>% length()
  
  fit_stan_MDS(mod = '../../Bayesian_stan/MDS_penalize_high_dims.stan',
               observed_distances = format_stan_inputs(titer_map),
               n_antigens = nsera,
               n_sera = nantigen,
               n_dim = n_dim,
               chains = n_chains,
               niter = n_iter,
               diagnostic_file = paste0('diagnostics/', idstring, '_', n_dim, 'D'))
  
}



## Fit the even immunodominance maps for 1:8 dimensions
even_inputs <- read_rds('even_immunodominance_inputs.rds')
even_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map = even_inputs$titer_map, 
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'penalized_even')
}
names(even_fit_list) = paste0(1:8, 'D')



## Fit the skewed immunodominance maps for 1:8 dimensions
skewed_inputs <- read_rds('skewed_immunodominance_inputs.rds')
skewed_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map = skewed_inputs$titer_map, 
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'penalized_skewed')
}
names(skewed_fit_list) = paste0(1:8, 'D')




## Fit the E1 complete immunodominance maps for 1:8 dimensions
E1_inputs <- read_rds('E1_complete_dominance_inputs.rds')
E1_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map = E1_inputs$titer_map, 
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'penalized_E1')
}
names(E1_fit_list) = paste0(1:8, 'D')




## Fit the E2 complete immunodominance maps for 1:8 dimensions
E2_inputs <- read_rds('E2_complete_dominance_inputs.rds')
E2_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map = E2_inputs$titer_map, 
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'penalized_E2')
}
names(E2_fit_list) = paste0(1:8, 'D')



## Fit the E3 complete immunodominance maps for 1:8 dimensions
E3_inputs <- read_rds('E3_complete_dominance_inputs.rds')
E3_fit_list <- foreach(this_ndim = 1:8) %do% {
  one_fit(titer_map = E3_inputs$titer_map, 
          n_dim = this_ndim, 
          n_chains = 4, 
          n_iter = 5000,
          idstring = 'penalized_E3')
}
names(E3_fit_list) = paste0(1:8, 'D')


if(!dir.exists('test_penalized_fits')) dir.create('test_penalized_fits')
write_rds(even_fit_list, 'test_penalized_fits/even_fit_list.rds')
write_rds(skewed_fit_list, 'test_penalized_fitsskewed_fit_list.rds')
write_rds(E1_fit_list, 'test_penalized_fitsE1_fit_list.rds')
write_rds(E2_fit_list, 'test_penalized_fitsE2_fit_list.rds')
write_rds(E3_fit_list, 'test_penalized_fitsE3_fit_list.rds')