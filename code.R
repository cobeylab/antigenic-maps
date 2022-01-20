## Functions to simulate maps and test MDS

library(tidyverse)
filter <- dplyr::filter
select <- dplyr::select
extract <- tidyr::extract

generate_map_coordinates <- function(n_dim, # Dimensions in Euclidean space 
                                     n_antigens, # Number of antigens in panel (n_antigens + n_sera must be > n_dim + 1)
                                     n_sera, # Number of sera/mAbs in panel 
                                     map_range = c(0, 1) # c(min, max) of each axis in the standard basis
){
  ## Draw antigen and serum locations from a uniform distribution across each axis
  nn = (n_antigens+n_sera)
  stopifnot(nn > (n_dim + 1))
  ## Draw a vector of coordinates
  coords = runif(nn*n_dim, map_range[1], map_range[2])
  ## Recast the vector into a data frame with n_dim columns
  coords = matrix(coords, nrow = nn, 
                  ncol = n_dim, 
                  dimnames = list(NULL, paste0('c', 1:n_dim))) %>%
    as.tibble()
  ## Bind into a tibble with id columns and return
  tibble(id = 1:nn,
         kind = rep(c('antigen', 'serum'), c(n_antigens, n_sera))) %>%
    bind_cols(coords)
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

get_all_distances <- function(coord_tibble){
  pointmat <- coord_tibble_to_matrix(coord_tibble)
  dist(pointmat, diag = T) %>%
    as.matrix()
}

get_obsered_distances <- function(coord_tibble){
  ## Only antigen-sera distances are observed
  ## This fun excludes ag-ag or serum-serum distances
  all_distances = get_all_distances(coord_tibble)
  is.ag = coord_tibble$kind == 'antigen'
  all_distances[is.ag, is.ag] = NA # Exclude ag-ag distances
  all_distances[!is.ag, !is.ag] = NA # Exclude ab-ab distances
  all_distances
}
