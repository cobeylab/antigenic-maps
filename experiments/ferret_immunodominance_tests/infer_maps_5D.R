## Infer maps for each set of inputs

rm(list = ls())
source('../../R/utility-functions.R')
source('../../R/stan_funs_sparse_data.R')
library(foreach)
library(doParallel)
library(rstan)
library(tidyverse)



## Run 3D fits for even immunodominance
fit_accross_dims(n_dim_inputs = 5, immunodominance_flag = 'even')
fit_accross_dims(n_dim_inputs = 5, immunodominance_flag = 'skewed')
fit_accross_dims(n_dim_inputs = 5, immunodominance_flag = 'E1')

