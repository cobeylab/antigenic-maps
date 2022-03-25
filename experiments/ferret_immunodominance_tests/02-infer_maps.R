## Infer maps for each set of inputs

rm(list = ls())
source('../../R/utility-functions.R')
source('../../R/stan_funs_sparse_data.R')
library(foreach)
library(doParallel)
library(rstan)
library(tidyverse)
if(!dir.exists('diagnostics')) dir.create('diagnostics')


## Run 2D fits for even immunodominance
t1 <- Sys.time()
fit_accross_dims(n_dim_inputs = 2, immunodominance_flag = 'even')
t2 <- Sys.time()
t2-t1
fit_accross_dims(n_dim_inputs = 2, immunodominance_flag = 'skewed')
t3 <- Sys.time()
t3-t2
fit_accross_dims(n_dim_inputs = 2, immunodominance_flag = 'E1')
Sys.time() - t3


## Run 3D fits for even immunodominance
t1 <- Sys.time()
fit_accross_dims(n_dim_inputs = 3, immunodominance_flag = 'even')
t2 <- Sys.time()
t2-t1
fit_accross_dims(n_dim_inputs = 3, immunodominance_flag = 'skewed')
t3 <- Sys.time()
t3-t2
fit_accross_dims(n_dim_inputs = 3, immunodominance_flag = 'E1')
Sys.time() - t3


## Run 5D fits for even immunodominance
t1 <- Sys.time()
fit_accross_dims(n_dim_inputs = 5, immunodominance_flag = 'even')
t2 <- Sys.time()
t2-t1
fit_accross_dims(n_dim_inputs = 5, immunodominance_flag = 'skewed')
t3 <- Sys.time()
t3-t2
fit_accross_dims(n_dim_inputs = 5, immunodominance_flag = 'E1')
Sys.time() - t3