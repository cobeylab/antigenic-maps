## Plot and analyze maps

rm(list = ls())
source('../../../R/code.R')
source('../../../R/multi_epitope_code_sparse_dataset.R')
library(foreach)
library(doParallel)
library(tidyverse)
summarise <- dplyr::summarise
rename <- dplyr::rename


## Get map error
even_fit_list <- read_rds('even_fit_list.rds')
even_inputs <- read_rds('../even_immunodominance_inputs.rds')
even_test_train = read_rds('even_inputs_test_train_split.rds')

skewed_fit_list <- read_rds('skewed_fit_list.rds')
skewed_inputs <- read_rds('../skewed_immunodominance_inputs.rds')
skewed_test_train <- read_rds('skewed_inputs_test_train_split.rds')

E1_fit_list <- read_rds('E1_fit_list.rds')
E1_inputs <- read_rds('../E1_complete_dominance_inputs.rds')
E1_test_train <- read_rds('E1_inputs_test_train_split.rds')

E2_inputs = read_rds('../E2_complete_dominance_inputs.rds')
E3_inputs = read_rds('../E3_complete_dominance_inputs.rds')

get_loo <- function(stan_fit){
  loo<-rstan::loo(stan_fit, pars = 'log_lik',
             cores = 4)
  as.numeric(loo$estimates['p_loo',1])
  # pareto_k_influence_values(loo)
  # plot(pareto_k_influence_values(loo))
}



analyze_one_output_set <- function(fit_list,
                                   test_train_list,
                                   immunodominance_flag){

## Posterior error
errors = lapply(fit_list, FUN = function(ll) get_titer_error(stan_fit = ll, titer_map = test_train_list$train)) 
  
whole_model_errors =  lapply(errors, function(ll) ll$model_errors) %>%
  bind_rows(.id = 'dimensions') %>%
  mutate(kind = immunodominance_flag)

whole_model_errors$loo = sapply(fit_list, get_loo)

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




map_plot <- lapply(fit_list, extract_summary_coords) %>%
  bind_rows(.id = 'dimensions') %>%
  select(dimensions, allele, kind, chain, matches('c\\d+$')) %>%
  ggplot() +
  geom_point(aes(x = c1, y = c2, color = chain)) +
  geom_point(aes(x = c1_Ag, y = c2_Ag), pch = 3, data = even_inputs$ag_ab_coords) +
  facet_wrap(.~dimensions)
if(!dir.exists('plots')) {
  dir.create('plots')
  cat('creating plots directory')
}
ggsave(sprintf('plots/inferred_map_%s_immunodominance.png', immunodominance_flag))

pca_results <- run_pca_on_titer_map(
  convert_titer_map_to_matrix(test_train_list$train)
)

return(list(whole_model_errors = whole_model_errors,
            whole_predictive_errors = whole_predictive_errors,
            pairwise_predictive_errors = pairwise_predictive_errors,
            pairwise_model_errors = pairwise_model_errors,
            map_plot = map_plot,
            pca_results = pca_results))
}

even_outputs <- analyze_one_output_set(fit_list = even_fit_list, test_train_list = even_test_train, immunodominance_flag = 'even')
skewed_outputs <- analyze_one_output_set(fit_list = skewed_fit_list, test_train_list = skewed_test_train, immunodominance_flag = 'skewed')
E1_outputs <- analyze_one_output_set(fit_list = E1_fit_list, test_train_list = E1_test_train, immunodominance_flag = 'E1')



## Plot distance error by dimension
bind_rows(list(even = even_outputs$whole_model_errors,
               skewed = skewed_outputs$whole_model_errors,
               E1 = E1_outputs$whole_model_errors),
          .id = 'kind') %>%
  ungroup() %>%
  tidyr::extract(dimensions, into = 'num_dimension', regex = '(\\d+)D', convert = T, remove = F) %>%
  arrange(kind, dimensions) %>%
  mutate(dimensions = factor(num_dimension, labels = unique(dimensions))) %>%
ggplot()+
  geom_ribbon(aes(x = num_dimension, ymin = distance_error_0.025, ymax = distance_error_0.975, fill = kind), alpha = .2) +
  geom_line(aes(x = num_dimension, y = distance_error_med, color = kind)) +
  scale_color_viridis_d(aesthetics = c('fill', 'color')) +
  xlab('dimensions') + ylab('model distance error')
ggsave('plots/model_error_distance.png', width = 7, height = 3.5, units = 'in', dpi = 200)


## Plot titer error by dimension
bind_rows(list(even = even_outputs$whole_model_errors,
               skewed = skewed_outputs$whole_model_errors,
               E1 = E1_outputs$whole_model_errors),
          .id = 'kind') %>%
  ungroup() %>%
  tidyr::extract(dimensions, into = 'num_dimension', regex = '(\\d+)D', convert = T, remove = F) %>%
  arrange(dimensions) %>%
  mutate(dimensions = factor(num_dimension, labels = unique(dimensions))) %>%
  ggplot()+
  geom_ribbon(aes(x = num_dimension, ymin = titer_error_0.025, ymax = titer_error_0.975, fill = kind), alpha = .2) +
  geom_line(aes(x = num_dimension, y = titer_error_med, color = kind)) +
  scale_color_viridis_d(aesthetics = c('fill', 'color')) +
  xlab('dimensions') + ylab('model titer error')
ggsave('plots/model_error_titer.png', width = 7, height = 3.5, units = 'in', dpi = 200)


## Plot loo by dimension
bind_rows(list(even = even_outputs$whole_model_errors,
               skewed = skewed_outputs$whole_model_errors,
               E1 = E1_outputs$whole_model_errors),
          .id = 'kind')  %>%
  tidyr::extract(dimensions, into = 'num_dimension', regex = '(\\d+)D', convert = T, remove = F) %>%
  arrange(dimensions) %>%
  mutate(dimensions = factor(num_dimension, labels = unique(dimensions))) %>%
  ggplot(aes(x = dimensions, y = loo, color = kind))+
  geom_point() +
  geom_line(aes(x = num_dimension, y = loo, color = kind)) +
  ggtitle('These may be unreliable')
ggsave('plots/model_loo.png', width = 7, height = 3.5, units = 'in', dpi = 200)


## Plot posterior predictive error by dimension
bind_rows(list(even = even_outputs$whole_predictive_errors,
               skewed = skewed_outputs$whole_predictive_errors,
               E1 = E1_outputs$whole_predictive_errors),
          .id = 'kind')  %>%
  tidyr::extract(dimensions, into = 'num_dimension', regex = '(\\d+)D', convert = T, remove = F) %>%
  ggplot() +
  geom_ribbon(aes(x = num_dimension, ymin = titer_error_0.025, ymax = titer_error_0.975, fill = kind), alpha = .2) +
  geom_line(aes(x = num_dimension, y = titer_error_med, color = kind)) +
  scale_color_viridis_d(aesthetics = c('fill', 'color')) +
  xlab('dimensions') + ylab('predictive titer error')
ggsave('plots/model_predictive_error.png', width = 7, height = 3.5, units = 'in', dpi = 200)


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
ggsave('plots/pairwise_distance_errors.png', width = 7, height = 3.5, units = 'in', dpi = 200)

## Plot pca analysis
cowplot::plot_grid(even_outputs$pca_results$var_explained_plot + ggtitle('even variance explained'),
                   skewed_outputs$pca_results$var_explained_plot + ggtitle('skewed variance explained'),
                   E1_outputs$pca_results$var_explained_plot + ggtitle('E1 variance explained'),
                   even_outputs$pca_results$biplot + ggtitle('even variance explained'),
                   skewed_outputs$pca_results$biplot + ggtitle('skewed variance explained'),
                   E1_outputs$pca_results$biplot + ggtitle('E1 variance explained'),
                   nrow = 2)
ggsave('plots/pca.png', width = 10, height = 6, units = 'in')




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
  scale_x_continuous(name = 'test set titers', breaks = 1:5, labels = sprintf('serum %s:\nantigen %s', ag_serum_map$serum, ag_serum_map$antigen))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('plots/generalized_E1_error.png', width = 7, height = 4, units = 'in', dpi = 300)



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
  scale_x_continuous(name = 'test set titers', breaks = 1:5, labels = sprintf('serum %s:\nantigen %s', ag_serum_map$serum, ag_serum_map$antigen))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('plots/generalized_even_error.png', width = 7, height = 4, units = 'in', dpi = 300)



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
  scale_x_continuous(name = 'test set titers', breaks = 1:5, labels = sprintf('serum %s:\nantigen %s', ag_serum_map$serum, ag_serum_map$antigen))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('plots/generalized_skewed_error.png', width = 7, height = 4, units = 'in', dpi = 300)

