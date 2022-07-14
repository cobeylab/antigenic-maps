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

this_ndim = 2

outdir = sprintf('plots/%s-%sD_MDS', Sys.Date(), this_ndim); outdir
output_dir_check(outdir) # Create the directory if it does not already exist


even_outputs = read_reds(sprintf('outputs/2D_even_parsed_outputs.rds'))
sparse_outputs = read_reds(sprintf('outputs/2D_sparse_parsed_outputs.rds'))
E1_outputs = read_reds(sprintf('outputs/2D_E1_parsed_outputs.rds'))



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
  scale_x_continuous(name = 'test set titers', breaks = 1:5, labels = sprintf('serum %s:\nantigen %s', ag_serum_map$serum, ag_serum_map$antigen))+
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
  scale_x_continuous(name = 'test set titers', breaks = 1:5, labels = sprintf('serum %s:\nantigen %s', ag_serum_map$serum, ag_serum_map$antigen))+
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
  scale_x_continuous(name = 'test set titers', breaks = 1:nrow(skewed_test_train$test), labels = sprintf('serum %s:\nantigen %s', ag_serum_map$serum, ag_serum_map$antigen))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(sprintf('%s/generalized_skewed_error.png', outdir), width = 7, height = 4, units = 'in', dpi = 300)









