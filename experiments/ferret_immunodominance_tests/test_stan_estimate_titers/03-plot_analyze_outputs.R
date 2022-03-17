## Plot and analyze maps

rm(list = ls())
source('../../../code.R')
source('../../../multi_epitope_code.R')
library(foreach)
library(doParallel)
library(tidyverse)
summarise <- dplyr::summarise
rename <- dplyr::rename


## Get map error
even_fit_list <- read_rds('even_fit_list.rds')
even_inputs <- read_rds('../even_immunodominance_inputs.rds')

skewed_fit_list <- read_rds('skewed_fit_list.rds')
skewed_inputs <- read_rds('../skewed_immunodominance_inputs.rds')

E1_fit_list <- read_rds('E1_fit_list.rds')
E1_inputs <- read_rds('../E1_complete_dominance_inputs.rds')


get_loo <- function(stan_fit){
  loo<-rstan::loo(stan_fit, pars = 'log_lik',
             cores = 4)
  as.numeric(loo$estimates['p_loo',1])
  # pareto_k_influence_values(loo)
  # plot(pareto_k_influence_values(loo))
}




analyze_one_output_set <- function(fit_list,
                                   input_list,
                                   immunodominance_flag){

pairwise_errors = lapply(fit_list, FUN = function(ll) get_pairwise_titer_error(stan_fit = ll, titer_map = input_list$titer_map)) %>%
  bind_rows(.id = 'dimensions') %>%
  mutate(kind = immunodominance_flag)
  
whole_model_errors =  pairwise_errors %>%
  ungroup() %>%
  group_by(dimensions) %>%
  summarise(distance_error = sum(pairwise_distance_error),
            titer_error = sum(pairwise_titer_error))
whole_model_errors$loo = sapply(fit_list, get_loo)

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
  convert_titer_map_to_matrix(input_list$titer_map)
)

return(list(pairwise_errors = pairwise_errors,
            whole_model_errors = whole_model_errors,
            map_plot = map_plot,
            pca_results = pca_results))
}

even_outputs <- analyze_one_output_set(fit_list = even_fit_list, input_list = even_inputs, immunodominance_flag = 'even')
skewed_outputs <- analyze_one_output_set(fit_list = skewed_fit_list, input_list = skewed_inputs, immunodominance_flag = 'skewed')
E1_outputs <- analyze_one_output_set(fit_list = E1_fit_list, input_list = E1_inputs, immunodominance_flag = 'E1')



## Plot error by dimension
bind_rows(list(even = even_outputs$whole_model_errors,
               skewed = skewed_outputs$whole_model_errors,
               E1 = E1_outputs$whole_model_errors),
          .id = 'kind') %>%
  ungroup() %>%
  tidyr::extract(dimensions, into = 'num_dimension', regex = '(\\d+)D', convert = T, remove = F) %>%
  arrange(dimensions) %>%
  mutate(dimensions = factor(num_dimension, labels = unique(dimensions))) %>%
ggplot(aes(x = dimensions, y = distance_error, color = kind))+
  geom_point() +
  geom_line(aes(x = as.numeric(dimensions), y = distance_error, color = kind))
ggsave('plots/model_error_distance.png')


## Plot error by dimension
bind_rows(list(even = even_outputs$whole_model_errors,
               skewed = skewed_outputs$whole_model_errors,
               E1 = E1_outputs$whole_model_errors),
          .id = 'kind') %>%
  ungroup() %>%
  tidyr::extract(dimensions, into = 'num_dimension', regex = '(\\d+)D', convert = T, remove = F) %>%
  arrange(dimensions) %>%
  mutate(dimensions = factor(num_dimension, labels = unique(dimensions))) %>%
  ggplot(aes(x = dimensions, y = titer_error, color = kind))+
  geom_point() +
  geom_line(aes(x = as.numeric(dimensions), y = titer_error, color = kind))
ggsave('plots/model_error_titer.png')


## Plot error by dimension
bind_rows(even_outputs$errors,
          skewed_outputs$errors,
          E1_outputs$errors) %>%
  ggplot(aes(x = dimensions, y = error, color = kind))+
  geom_point() +
  geom_line()
ggsave('plots/model_error_titer.png')


## Plot error by dimension
bind_rows(list(even = even_outputs$whole_model_errors,
               skewed = skewed_outputs$whole_model_errors,
               E1 = E1_outputs$whole_model_errors),
          .id = 'kind')  %>%
  ggplot(aes(x = dimensions, y = loo, color = kind))+
  geom_point() +
  geom_line() +
  ggtitle('These may be unreliable')
ggsave('plots/model_loo.png')

## Look at inferred coordinates by dimension
bind_rows(list(even = even_outputs$pairwise_errors,
               skewed = skewed_outputs$pairwise_errors,
               E1 = E1_outputs$pairwise_errors),
          .id = 'kind') %>%
  mutate(pairwise_distance_errors = mean_map_distance-titer_distance) %>%
  ggplot() +
  geom_point(aes(x = antigen, y = pairwise_distance_errors, color = dimensions), alpha = .5) +
  geom_hline(aes(yintercept = 0), lty = 2)+
  facet_grid(kind~serum)
ggsave('plots/pairwise_distance_errors.png')  

## Plot pca analysis
cowplot::plot_grid(even_outputs$pca_results$var_explained_plot + ggtitle('even variance explained'),
                   skewed_outputs$pca_results$var_explained_plot + ggtitle('skewed variance explained'),
                   E1_outputs$pca_results$var_explained_plot + ggtitle('E1 variance explained'),
                   even_outputs$pca_results$biplot + ggtitle('even variance explained'),
                   skewed_outputs$pca_results$biplot + ggtitle('skewed variance explained'),
                   E1_outputs$pca_results$biplot + ggtitle('E1 variance explained'),
                   nrow = 2)
ggsave('plots/pca.png', width = 10, height = 6, units = 'in')


