## Plot and analyze maps

rm(list = ls())
source('../../R/utility-functions.R')
source('../../R/stan_funs_sparse_data.R')
library(foreach)
library(doParallel)
library(tidyverse)
if(!dir.exists('plots')) dir.create('plots')
summarise <- dplyr::summarise
rename <- dplyr::rename


############# Set this number of dimensions ##############
this_ndim = 3

outdir = sprintf('plots/%s--%sD_MDS', Sys.Date(), this_ndim); outdir
output_dir_check(outdir) # Create the directory if it does not already exist


## Function to parse all outputs
analyze_one_output_set <- function(fit_list,
                                   test_train_list,
                                   immunodominance_flag,
                                   outdir){
  
  ## Posterior error
  errors = lapply(fit_list, FUN = function(ll) get_titer_error(stan_fit = ll, titer_map = test_train_list$train)) 
  whole_model_errors =  lapply(errors, function(ll) ll$model_errors) %>%
    bind_rows(.id = 'dimensions') %>%
    mutate(kind = immunodominance_flag)
  pairwise_model_errors =  lapply(errors, function(ll) ll$pairwise_errors) %>%
    bind_rows(.id = 'dimensions') %>%
    mutate(kind = immunodominance_flag)
  
  ## Posterior predictive error
  predictive_errors = lapply(fit_list, FUN = function(ll) get_predictive_errors(stan_fit = ll, test_set = test_train_list$test)) 
  whole_predictive_errors = lapply(predictive_errors, function(ll) ll$model_errors) %>%
    bind_rows(.id = 'dimensions') %>%
    mutate(kind = immunodominance_flag)
  pairwise_predictive_errors =  lapply(predictive_errors, function(ll) ll$pairwise_errors) %>%
    bind_rows(.id = 'dimensions') %>%
    mutate(kind = immunodominance_flag)
  
  pca_results <- run_pca_on_titer_map(
    convert_titer_map_to_matrix(test_train_list$train)
  )

  return(list(whole_model_errors = whole_model_errors,
              whole_predictive_errors = whole_predictive_errors,
              pairwise_predictive_errors = pairwise_predictive_errors,
              pairwise_model_errors = pairwise_model_errors,
              pca_results = pca_results))
}


plot_inferred_map <- function(fit_list, 
                              ab_ag_coords,
                              outdir, 
                              immunodominance_flag){ 
  ## Plot the inferred map
  map_plot <- lapply(fit_list, extract_summary_coords) %>%
  bind_rows(.id = 'dimensions') %>%
  select(dimensions, allele, kind, chain, matches('c\\d+$')) %>%
  ggplot() +
  geom_point(aes(x = c1, y = c2, color = chain)) +
  geom_point(aes(x = c1_Ag, y = c2_Ag), pch = 3, data = ab_ag_coords) +
  facet_wrap(.~dimensions)
map_plot
ggsave(sprintf('%s/inferred_map_%s_immunodominance.png', outdir, immunodominance_flag), width = 7, height=5, units = 'in', dpi = 200)
}



## Import saved inputs and fits
this_immunodominance = 'even'
even_fit_list <- read_rds(sprintf('outputs/%sDinputs-%s-fit_list.rds', this_ndim, this_immunodominance))
even_inputs <- read_rds(sprintf('inputs/%sD_%s_immunodominance_inputs.rds', this_ndim, this_immunodominance))
even_test_train <- read_rds(sprintf('outputs/%sDinputs-%s-test_train_split.rds', this_ndim, this_immunodominance))

this_immunodominance = 'skewed'
skewed_fit_list <- read_rds(sprintf('outputs/%sDinputs-%s-fit_list.rds', this_ndim, this_immunodominance))
skewed_inputs <- read_rds(sprintf('inputs/%sD_%s_immunodominance_inputs.rds', this_ndim, this_immunodominance))
skewed_test_train <- read_rds(sprintf('outputs/%sDinputs-%s-test_train_split.rds', this_ndim, this_immunodominance))

this_immunodominance = 'E1'
E1_fit_list <- read_rds(sprintf('outputs/%sDinputs-%s-fit_list.rds', this_ndim, this_immunodominance))
E1_inputs <- read_rds(sprintf('inputs/%sD_%s_immunodominance_inputs.rds', this_ndim, this_immunodominance))
E1_test_train <- read_rds(sprintf('outputs/%sDinputs-%s-test_train_split.rds', this_ndim, this_immunodominance))

E2_inputs = read_rds(sprintf('inputs/%sD_E2_immunodominance_inputs.rds', this_ndim))
E3_inputs = read_rds(sprintf('inputs/%sD_E3_immunodominance_inputs.rds', this_ndim))


## Parse the outputs
even_outputs <- analyze_one_output_set(fit_list = even_fit_list, 
                                       test_train_list = even_test_train, 
                                       immunodominance_flag = 'even',
                                       outdir)
skewed_outputs <- analyze_one_output_set(fit_list = skewed_fit_list, 
                                         test_train_list = skewed_test_train, 
                                         immunodominance_flag = 'skewed',
                                         outdir)
E1_outputs <- analyze_one_output_set(fit_list = E1_fit_list, 
                                         test_train_list = E1_test_train, 
                                         immunodominance_flag = 'E1',
                                         outdir)






## Plot distance error by dimension
bind_rows(list(even = even_outputs$whole_model_errors,
               skewed = skewed_outputs$whole_model_errors,
               E1 = E1_outputs$whole_model_errors),
          .id = 'kind') %>%
  ungroup() %>%
  arrange(dimensions) %>%
  mutate(x_location = as.numeric(as.factor(dimensions)) + (as.numeric(as.factor(kind))-2)/5)%>% 
  tidyr::extract(dimensions, into = 'num_dimension', regex = '(\\d+)D', convert = T, remove = F) %>%
  arrange(kind, dimensions) %>%
  mutate(dimensions = factor(num_dimension, labels = unique(dimensions))) %>%
  ggplot()+
  geom_point(aes(x = x_location, y = distance_error_med, color = kind)) +
  geom_segment(aes(x = x_location, xend = x_location, y = distance_error_0.025, yend = distance_error_0.975, color = kind)) +
  xlab('dimensions') + ylab('model distance error')
ggsave(sprintf('%s/model_error_distance.png', outdir), width = 7, height = 3.5, units = 'in', dpi = 200)


## Plot distance error by dimension, per ag-ab pair in training set
n_pairs = nrow(even_test_train$train)
bind_rows(list(even = even_outputs$whole_model_errors,
               skewed = skewed_outputs$whole_model_errors,
               E1 = E1_outputs$whole_model_errors),
          .id = 'kind') %>%
  ungroup() %>%
  arrange(dimensions) %>%
  mutate(x_location = as.numeric(as.factor(dimensions)) + (as.numeric(as.factor(kind))-2)/5)%>% 
  tidyr::extract(dimensions, into = 'num_dimension', regex = '(\\d+)D', convert = T, remove = F) %>%
  arrange(kind, dimensions) %>%
  mutate(dimensions = factor(num_dimension, labels = unique(dimensions))) %>%
  ggplot()+
  geom_point(aes(x = x_location, y = distance_error_med/n_pairs, color = kind)) +
  geom_segment(aes(x = x_location, xend = x_location, y = distance_error_0.025/n_pairs, yend = distance_error_0.975/n_pairs, color = kind)) +
  xlab('dimensions') + ylab('mean error per ag-ab pair\n(distance units)')
ggsave(sprintf('%s/model_error_distance_mean_pairwise.png', outdir), width = 7, height = 3.5, units = 'in', dpi = 200)


## Plot titer error by dimension
bind_rows(list(even = even_outputs$whole_model_errors,
               skewed = skewed_outputs$whole_model_errors,
               E1 = E1_outputs$whole_model_errors),
          .id = 'kind') %>%
  ungroup() %>%
  arrange(dimensions) %>%
  mutate(x_location = as.numeric(as.factor(dimensions)) + (as.numeric(as.factor(kind))-2)/5)%>% 
  tidyr::extract(dimensions, into = 'num_dimension', regex = '(\\d+)D', convert = T, remove = F) %>%
  arrange(kind, dimensions) %>%
  mutate(dimensions = factor(num_dimension, labels = unique(dimensions))) %>%
  ggplot()+
  geom_point(aes(x = x_location, y = titer_error_med, color = kind)) +
  geom_segment(aes(x = x_location, xend = x_location, y = titer_error_0.025, yend = titer_error_0.975, color = kind)) +
  xlab('dimensions') + ylab('model titer error')
ggsave(sprintf('%s/model_error_titer.png', ourdir), width = 7, height = 3.5, units = 'in', dpi = 200)

## Plot titer error by dimension
bind_rows(list(even = even_outputs$whole_model_errors,
               skewed = skewed_outputs$whole_model_errors,
               E1 = E1_outputs$whole_model_errors),
          .id = 'kind') %>%
  ungroup() %>%
  arrange(dimensions) %>%
  mutate(x_location = as.numeric(as.factor(dimensions)) + (as.numeric(as.factor(kind))-2)/5)%>% 
  tidyr::extract(dimensions, into = 'num_dimension', regex = '(\\d+)D', convert = T, remove = F) %>%
  arrange(kind, dimensions) %>%
  mutate(dimensions = factor(num_dimension, labels = unique(dimensions))) %>%
  ggplot()+
  geom_point(aes(x = x_location, y = titer_error_med/n_pairs, color = kind)) +
  geom_segment(aes(x = x_location, xend = x_location, y = titer_error_0.025/n_pairs, yend = titer_error_0.975/n_pairs, color = kind)) +
  xlab('dimensions') + ylab('model titer error')
ggsave(sprintf('%s/model_error_titer_mean_pairwise.png', ourdir), width = 7, height = 3.5, units = 'in', dpi = 200)



## Plot posterior predictive error by dimension
bind_rows(list(even = even_outputs$whole_predictive_errors,
               skewed = skewed_outputs$whole_predictive_errors,
               E1 = E1_outputs$whole_predictive_errors),
          .id = 'kind')  %>%
  ungroup() %>%
  arrange(dimensions) %>%
  mutate(x_location = as.numeric(as.factor(dimensions)) + (as.numeric(as.factor(kind))-2)/5)%>% 
  ggplot()+
  geom_point(aes(x = x_location, y = titer_error_med, color = kind)) +
  geom_segment(aes(x = x_location, xend = x_location, y = titer_error_0.025, yend = titer_error_0.975, color = kind)) +
  xlab('dimensions') + ylab('predictive titer error') 
ggsave(sprintf('%s/model_predictive_error.png', outdir), width = 7, height = 3.5, units = 'in', dpi = 200)

## Plot posterior predictive error by dimension
n_test_pairs = nrow(even_test_train$test)
bind_rows(list(even = even_outputs$whole_predictive_errors,
               skewed = skewed_outputs$whole_predictive_errors,
               E1 = E1_outputs$whole_predictive_errors),
          .id = 'kind')  %>%
  ungroup() %>%
  arrange(dimensions) %>%
  mutate(x_location = as.numeric(as.factor(dimensions)) + (as.numeric(as.factor(kind))-2)/5)%>% 
  ggplot()+
  geom_point(aes(x = x_location, y = titer_error_med/n_test_pairs, color = kind)) +
  geom_segment(aes(x = x_location, xend = x_location, y = titer_error_0.025/n_test_pairs, yend = titer_error_0.975/n_test_pairs, color = kind)) +
  xlab('dimensions') + ylab('predictive titer error') +
  geom_hline(aes(yintercept = 1), lty = 3)
ggsave(sprintf('%s/model_predictive_error_mean_pairwise.png', outdir), width = 7, height = 3.5, units = 'in', dpi = 200)


## Look at inferred coordinates by dimension
bind_rows(list(even = even_outputs$pairwise_model_errors,
               skewed = skewed_outputs$pairwise_model_errors,
               E1 = E1_outputs$pairwise_model_errors),
          .id = 'kind') %>%
  ggplot() +
  geom_point(aes(x = antigen, y = distance_error_med, color = dimensions), alpha = .5) +
  geom_hline(aes(yintercept = 0), lty = 2)+
  facet_grid(kind~serum) +
  ylab('median distance error')
ggsave(sprintf('%s/pairwise_distance_errors.png', outdir), width = 7, height = 3.5, units = 'in', dpi = 200)



## Plot pca analysis
cowplot::plot_grid(even_outputs$pca_results$var_explained_plot + ggtitle('even variance explained'),
                   skewed_outputs$pca_results$var_explained_plot + ggtitle('skewed variance explained'),
                   E1_outputs$pca_results$var_explained_plot + ggtitle('E1 variance explained'),
                   even_outputs$pca_results$biplot + ggtitle('even variance explained'),
                   skewed_outputs$pca_results$biplot + ggtitle('skewed variance explained'),
                   E1_outputs$pca_results$biplot + ggtitle('E1 variance explained'),
                   nrow = 2)
ggsave(sprintf('%s/pca.png', outdir), width = 10, height = 6, units = 'in')






## Calculate error between E1 predictions and all other actual titers
all_observed_titers = list(E1 = E1_inputs$titer_map,
                           E2 = E2_inputs$titer_map,
                           E3 = E3_inputs$titer_map,
                           skewed = skewed_inputs$titer_map,
                           even = even_inputs$titer_map)
E1_temp <- foreach(these_observed_titers = all_observed_titers,
                   this_observed_scheme = names(all_observed_titers)) %do% {
  get_generalized_titer_errors(stan_fit = E1_fit_list[[2]], 
                               test_set = E1_test_train$test, 
                               observed_titers = these_observed_titers, fitted_immunodominance_scheme = 'E1', 
                               observed_immunodominance_scheme = this_observed_scheme)
}
E1_generalized_overall_error = lapply(E1_temp, function(ll) ll$overall_errors) %>% bind_rows()
E1_generalized_pairwise = lapply(E1_temp, function(ll) ll$pairwise_errors) %>% bind_rows()
ag_serum_map = E1_generalized_pairwise %>% arrange(id) %>% group_by(antigen, serum, id) %>% summarise()
E1_generalized_pairwise %>%
  arrange(id) %>%
  mutate(x_location = as.numeric(id) + (as.numeric(as.factor(kind))-3)/10) %>%
ggplot() +
  geom_point(aes(x = x_location, y = titer_error_med, color = kind))+
  geom_segment(aes(x = x_location, xend = x_location, y = titer_error_0.025, yend = titer_error_0.975, color = kind)) +
  ylab('posterior predictive error') +
  xlab('') +
  geom_hline(aes(yintercept = 1), lty = 2) +
  scale_x_continuous(name = 'test set titers', breaks = 1:n_test_pairs, labels = sprintf('serum %s:\nantigen %s', ag_serum_map$serum, ag_serum_map$antigen))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf('%s/generalized_E1_error.png', outdir), width = 7, height = 4, units = 'in', dpi = 300)



## REPEAT FOR EVEN FITS
even_temp <- foreach(these_observed_titers = all_observed_titers,
                   this_observed_scheme = names(all_observed_titers)) %do% {
                     get_generalized_titer_errors(stan_fit = even_fit_list[[2]], 
                                                  test_set = even_test_train$test, 
                                                  observed_titers = these_observed_titers, fitted_immunodominance_scheme = 'even', 
                                                  observed_immunodominance_scheme = this_observed_scheme)
                   }
even_generalized_overall_error = lapply(even_temp, function(ll) ll$overall_errors) %>% bind_rows()
even_generalized_pairwise = lapply(even_temp, function(ll) ll$pairwise_errors) %>% bind_rows()
ag_serum_map = even_generalized_pairwise %>% arrange(id) %>% group_by(antigen, serum, id) %>% summarise()
even_generalized_pairwise %>%
  arrange(id) %>%
  mutate(x_location = as.numeric(id) + (as.numeric(as.factor(kind))-3)/10) %>%
  ggplot() +
  geom_point(aes(x = x_location, y = titer_error_med, color = kind))+
  geom_segment(aes(x = x_location, xend = x_location, y = titer_error_0.025, yend = titer_error_0.975, color = kind)) +
  ylab('posterior predictive error') +
  xlab('') +
  geom_hline(aes(yintercept = 1), lty = 2) +
  scale_x_continuous(name = 'test set titers', breaks = 1:n_test_pairs, labels = sprintf('serum %s:\nantigen %s', ag_serum_map$serum, ag_serum_map$antigen))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf('%s/generalized_even_error.png', outdir), width = 7, height = 4, units = 'in', dpi = 300)



# REPEAT FOR SKEWED FITS
skewed_temp <- foreach(these_observed_titers = all_observed_titers,
                   this_observed_scheme = names(all_observed_titers)) %do% {
                     get_generalized_titer_errors(stan_fit = skewed_fit_list[[2]], 
                                                  test_set = skewed_test_train$test, 
                                                  observed_titers = these_observed_titers, fitted_immunodominance_scheme = 'skewed', 
                                                  observed_immunodominance_scheme = this_observed_scheme)
                   }
skewed_generalized_overall_error = lapply(skewed_temp, function(ll) ll$overall_errors) %>% bind_rows()
skewed_generalized_pairwise = lapply(skewed_temp, function(ll) ll$pairwise_errors) %>% bind_rows()
ag_serum_map = skewed_generalized_pairwise %>% arrange(id) %>% group_by(antigen, serum, id) %>% summarise()
skewed_generalized_pairwise %>%
  arrange(id) %>%
  mutate(x_location = as.numeric(id) + (as.numeric(as.factor(kind))-3)/10) %>%
  ggplot() +
  geom_point(aes(x = x_location, y = titer_error_med, color = kind))+
  geom_segment(aes(x = x_location, xend = x_location, y = titer_error_0.025, yend = titer_error_0.975, color = kind)) +
  ylab('posterior predictive error') +
  xlab('') +
  geom_hline(aes(yintercept = 1), lty = 2) +
  scale_x_continuous(name = 'test set titers', breaks = 1:n_test_pairs, labels = sprintf('serum %s:\nantigen %s', ag_serum_map$serum, ag_serum_map$antigen))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf('%s/generalized_skewed_error.png', outdir), width = 7, height = 4, units = 'in', dpi = 300)

