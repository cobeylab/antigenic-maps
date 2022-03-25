## Functions to simulate maps and test MDS

library(tidyverse)
filter <- dplyr::filter
select <- dplyr::select
extract <- tidyr::extract

zeroish <- function(x){
  abs(x) <= .Machine$double.eps^.5
}

equal_ish <- function(xx, yy = 0){
  abs(xx-yy) <= .Machine$double.eps^.5
}

ndim <- function(xx){
  ldim <- length(dim(xx)) ## Unlike dim, will print 0 if a vector
  if(ldim >= 2){
    return(ldim)
  }else{
    ifelse(length(dim) > 1, 1, 'scalar')
  }
}

test_train_split <- function(N,
                             n_test,
                             df){
  stopifnot(nrow(df)==N)
  test_inds = sample(1:N, size = n_test, replace = F)
  return(list(test = df[test_inds,],
              train = df[-test_inds,]
))
}


coord_tibble_to_matrix <- function(coord_tibble){
  coord_tibble %>%
    select(starts_with('c')) %>%
    as.matrix()
}

coord_matrix_to_tibble <- function(coord_matrix, n_antigens, n_sera){
    tibble(id = 1:(n_antigens+n_sera),
           kind = rep(c('antigen', 'serum'), c(n_antigens, n_sera))) %>%
      bind_cols(coord_matrix)
}


get_euclidean_distance<- function(v1, v2){
  stopifnot(length(v1) == length(v2))
  sqrt(sum((v1-v2)^2))
}


get_ab_ag_distances <- function(ab_coords, # matrix of Ab coords. Each column is a dimension. Each row is an Ab.
                                ag_coords # matrix of Ag coords. Each column is a dimension. Each row is an Ag.
                                ){
  distmat = matrix(NA, nrow = nrow(ag_coords), ncol = nrow(ab_coords))
  rownames(distmat) = paste0('ag', 1:nrow(ag_coords))
  colnames(distmat) = paste0('ab', 1:nrow(ab_coords))
  for(ab in 1:nrow(ab_coords)){
    for(ag in 1:nrow(ag_coords)){
      distmat[ag, ab] = get_euclidean_distance(ab_coords[ab,], ag_coords[ag,])
    }
  }
  distmat
}


output_dir_check <- function(outdir){
  if(!dir.exists(outdir)) {
    dir.create(outdir)
    cat(sprintf('creating %s directory', outdir))
  }
}
