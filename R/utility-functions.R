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



get_titer_epitope_level <- function(ab_ag_coords,
                                    this_antigen,
                                    this_serum,
                                    this_epitope = 1,
                                    alpha = .25, 
                                    r = 7){
  ## For each antigen, serum combination extract relevant coordinates
  these_inputs <- extract_titer_inputs(antigen_id = this_antigen, 
                                       serum_id = this_serum, 
                                       merged_df = ab_ag_coords)
  ## Get epitope ditances
  these_epitope_distances = get_ab_ag_distances(ab_coords = as.matrix(these_inputs$ab_position_list[[this_epitope]]), 
                                                ag_coords = matrix(these_inputs$ag_list[[this_epitope]], ncol = n_dim))
  ## Calculate titer and return
  (1-2*alpha)*sum(2^(r-these_epitope_distances))
}
