## Plot and analyze maps

rm(list = ls())
source('../../code.R')
source('../../multi_epitope_code.R')
library(foreach)
library(doParallel)


## Get map error
even_fit_list <- read_rds('even_fit_list.rds')
even_inputs <- read_rds('even_immunodominance_inputs.rds')

fits = lapply(read_rds('E1_fit_list.rds'), FUN = function(ll) extract_summary_coords(stan_fit = ll))

errors_df = lapply(even_fit_list, FUN = function(ll) get_map_error_df(stan_fit = ll, titer_map = even_inputs$titer_map)) %>%
  bind_rows(.id = 'dimensions') %>%
  mutate(kind = 'even')
  
errors = tibble(dimensions = 1:8,
                error = sapply(even_fit_list, FUN = function(ll) get_map_error(stan_fit = ll, titer_map = even_inputs$titer_map)),
                kind = 'even'
)

lapply(even_fit_list, extract_summary_coords) %>%
  bind_rows(.id = 'dimensions') %>%
  select(dimensions, allele, kind, chain, matches('c\\d+$')) %>%
  ggplot() +
  geom_point(aes(x = c1, y = c2, color = chain)) +
  geom_point(aes(x = c1_Ag, y = c2_Ag), pch = 3, data = even_inputs$ag_ab_coords) +
  facet_wrap(.~dimensions)


## Plot error by dimension
ggplot(errors)+
  geom_point(aes(x = dimensions, y = error))

## Look at inferred coordinates by dimension
errors_df %>%
  mutate(`map-observed distance` = predicted_map_distance-observed_distance) %>%
  ggplot() +
  geom_point(aes(x = antigen, y = `map-observed distance`, color = dimensions), alpha = .5) +
  geom_hline(aes(yintercept = 0), lty = 2)+
  facet_grid(.~serum)
  
  
