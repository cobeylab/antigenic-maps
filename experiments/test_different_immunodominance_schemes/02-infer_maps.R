## Infer maps for each set of inputs

rm(list = ls())
source('../../code.R')
source('../../multi_epitope_code_sparse_dataset.R')
library(foreach)
library(doParallel)
library(rstan)
library(tidyverse)
if(!dir.exists('diagnostics')) dir.create('diagnostics')
filter <- dplyr::filter



one_fit <- function(n_dim,
                    titer_map_train,
                    titer_map_test,
                    n_chains = 4,
                    n_iter= 5000,
                    idstring){
  
  nsera = unique(titer_map_train$serum) %>% length()
  nantigen = unique(titer_map_train$antigen) %>% length()
  cat(sprintf('in R: n_antigens is %i; n_sera is %i \n', nantigen, nsera))
  
  stopifnot(all(titer_map_test$serum %in% titer_map_train$serum))
  stopifnot(all(titer_map_test$antigen %in% titer_map_train$antigen))
  
  fit_stan_MDS(mod = '../../stan/MDS_predict_titer_sparse.stan',
               observed_titers = titer_map_train$logtiter, # N-vector of observed titers
               smax = (titer_map_train$serum_potency+titer_map_train$antigen_avidity)/2, # N-vector of smax
               serum_id = as.numeric(as.factor(titer_map_train$serum)), # N-vector of serum id
               antigen_id = titer_map_train$antigen, # N- vector of antigen id
               n_antigens = nantigen, # integer
               n_sera = nsera, # integer
               n_dim = n_dim, # integer
               N_test_set = nrow(titer_map_test),
               smax_test_set = (titer_map_test$serum_potency + titer_map_test$antigen_avidity)/2,
               serum_id_test_set = as.numeric(as.factor(titer_map_test$serum)),
               antigen_id_test_set = titer_map_test$antigen,
               observed_titers_test_set = titer_map_test$logtiter,
               chains = 3, # Number of MCMC chains to run
               cores = parallel::detectCores(logical = F), # For the cluster
               niter = 5000,
               diagnostic_file = NULL
  )
  
}

one_fit_wrapper <- function(serum_catalog = read_rds('simulated_inputs/serum_catalog.rds'), 
                            full_titer_map = read_rds('simulated_inputs/titer_inputs.rds'), 
                            valid_serum_ids, 
                            run_name){
  this.titer.map <-  titer_map %>%
    filter(serum %in% valid_serum_ids)
  ## Fit the even immunodominance maps for 1:8 dimensions
  serum_list[[this.titer.map$valid_serum_ids]] %>%
    bind_rows() %>%
  set.seed(11)
  split_titer_map <- test_train_split(nrow(this.titer.map), floor(nrow(this.titer.map)*.2), this.titer.map)
  stopifnot(this.titer.map$antigen %in% split_titer_map$train$antigen)
  stopifnot(this.titer.map$serum %in% split_titer_map$train$serum)
  fit_list <- foreach(this_ndim = 1:8) %do% {
    cat(sprintf('Running %iD fit', this_ndim))
    one_fit(titer_map_train = split_titer_map$train, 
            titer_map_test = split_titer_map$test,
            n_dim = this_ndim, 
            n_chains = 4, 
            n_iter = 5000,
            idstring = run_name)
  }
  names(fit_list) = paste0(1:8, 'D')
  if(!dir.exists('outputs')) dir.create('outputs')
  write_rds(fit_list, paste0('outputs/', run_name, '_fit_list.rds'))
  write_rds(split_titer_map, paste0('outputs/', run_name, '_test_train_split.rds'))
}



## Generate different cases to test
titer_map <- read_rds('simulated_inputs/titer_inputs.rds')
serum_catlog <- read_rds('simulated_inputs/serum_catalog.rds')
serum_list <- read_rds('simulated_inputs/serum_list.rds')



## 1. make a list of serum panels and inputs to run on the cluster
run_list <- list(
  #E1_E1_allornothing
  list(this_catalog_selection = serum_catlog %>% 
         filter(dominant_epitope == 1, immunodominance_scheme == 'all-or-nothing', replicate <= 2),
       run_name = 'E1_E1_allornothing'),
  #E1_E2_allornothing
  list(this_catalog_selection = serum_catlog %>% 
         filter(dominant_epitope %in% c(1,2), immunodominance_scheme == 'all-or-nothing', replicate == 2),
       run_name = 'E1_E2_allornothing'),
  #E1_E1_skewed
  list(this_catalog_selection = serum_catlog %>% 
         filter(dominant_epitope == 1, immunodominance_scheme == 'skewed', replicate <= 2),
       run_name = 'E1_E1_skewed'),
  #E1_E2_skewed
  list(this_catalog_selection = serum_catlog %>% 
         filter(dominant_epitope %in% c(1,2), immunodominance_scheme == 'skewed', replicate == 2),
       run_name = 'E1_E2_skewed'),
  #even
  list(this_catalog_selection = serum_catlog %>% 
         filter(dominant_epitope %in% c(1), immunodominance_scheme == 'even', replicate <= 2),
       run_name = 'even')
  )


# plot_one_serum_panel <- function(run_list_element){
# serum_list[run_list_element$this_catalog_selection$serum_id] %>%
#   bind_rows() %>%
#   ggplot() +
#   #geom_point(aes(x = c1_Ab, y = c2_Ab), pch =3, alpha = .5)+
#   geom_density_2d_filled(aes(x = c1_Ab, y = c2_Ab)) +
#   geom_point(aes(x = c1_Ag, y = c2_Ag)) +
#   facet_wrap(.~serum_id, labeller = label_both) +
#   ggtitle(sprintf('%s: serum antibody density', run_list_element$run_name))+
#   scale_fill_viridis_d(option = 'rocket', direction = -1)
#   ggsave(paste0('plots/', run_list_element$run_name, '.png'))
# }
# 
# lapply(run_list, plot_one_serum_panel)


one_fit_wrapper(serum_catalog = read_rds('simulated_inputs/serum_catalog.rds'), 
                full_titer_map = read_rds('simulated_inputs/titer_inputs.rds'), 
                valid_serum_ids = run_list[[this_run]]$this_catalog_selection$serum_id, 
                run_name = run_list[[this_run]]$run_name)
