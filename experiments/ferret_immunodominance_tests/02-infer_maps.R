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
cat(sprintf('\n \nstart even 2D: \n \n'))
fit_accross_dims(n_dim_inputs = 2, immunodominance_flag = 'even')
#read_rds('outputs/2Dinputs-even-fit_list.rds')
cat(sprintf('\n \nstart skewed 2D: \n \n'))
fit_accross_dims(n_dim_inputs = 2, immunodominance_flag = 'skewed')
#read_rds('outputs/2Dinputs-skewed-fit_list.rds')
cat(sprintf('\n \nstart E1 2D: \n \n'))
fit_accross_dims(n_dim_inputs = 2, immunodominance_flag = 'E1')
#read_rds('outputs/2Dinputs-E1-fit_list.rds')



## Run 3D fits for even immunodominance
cat(sprintf('\n \nstart even 3D: \n \n'))
fit_accross_dims(n_dim_inputs = 3, immunodominance_flag = 'even')
cat(sprintf('\n \nstart skewed 3D: \n \n'))
fit_accross_dims(n_dim_inputs = 3, immunodominance_flag = 'skewed')
cat(sprintf('\n \nstart E1 3D: \n \n'))
fit_accross_dims(n_dim_inputs = 3, immunodominance_flag = 'E1')



## Run 5D fits for even immunodominance
cat(sprintf('\n \nstart even 5D: \n \n'))
fit_accross_dims(n_dim_inputs = 5, immunodominance_flag = 'even')
cat(sprintf('\n \nstart skewed 5D: \n \n'))
fit_accross_dims(n_dim_inputs = 5, immunodominance_flag = 'skewed')
cat(sprintf('\n \nstart E1 5D: \n \n'))
fit_accross_dims(n_dim_inputs = 5, immunodominance_flag = 'E1')
