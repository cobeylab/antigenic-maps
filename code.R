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
  tibble(id = c(1:n_antigens, 1:n_sera),
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


get_euclidean_distance<- function(v1, v2){
  stopifnot(length(v1) == length(v2))
  sqrt(sum((v1-v2)^2))
}


get_ab_ag_distances <- function(ab_coords, ag_coords){
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


generate_initial_guess <- function(n_antigens, n_sera, n_dim){
  ag_c1s <- runif(n_antigens, -10, 10); names(ag_c1s) = paste0('ag', 1:n_antigens, 'c1')
  ag_c2s <- runif(n_antigens, -10, 10); names(ag_c2s) = paste0('ag', 1:n_antigens, 'c2')
  ab_c1s <- runif(n_sera, -10, 10); names(ab_c1s) = paste0('ab', 1:n_sera, 'c1')
  ab_c2s <- runif(n_sera, -10, 10); names(ab_c2s) = paste0('ab', 1:n_sera, 'c2')
  
  c(ag_c1s, ag_c2s, ab_c1s, ab_c2s)
}


fit_MDS_least_squares <- function(
  observed_distances, # vector of delay for each individual in the data
  n_antigens,
  n_sera,
  n_dim
) {
  
  stopifnot(nrow(observed_distances) == n_antigens)
  stopifnot(ncol(observed_distances) == n_sera)
  
  loss_function <- function(pars, observed_distances, n_antigens, n_sera, n_dim){
    stopifnot(length(pars) == (n_antigens+n_sera)*n_dim)
    
    estimated_ag_coords <- matrix(pars[1:(n_antigens*n_dim)], n_antigens, n_dim)
    estimated_ab_coords <- matrix(pars[(n_antigens*n_dim+1):length(pars)], n_sera, n_dim)
    
    estimated_distances <- get_ab_ag_distances(estimated_ab_coords, estimated_ag_coords)
    sum((observed_distances - estimated_distances)^2)
  }
  
  ## Repeat the fit 10 times from different initial values to verify convergence
  fit_list <- lapply(1:10, function(x) { optim(par = generate_initial_guess(n_antigens, n_sera, n_dim),
                                               fn = loss_function,
                                               observed_distances = observed_distances,
                                               n_antigens = n_antigens, 
                                               n_sera = n_sera,
                                               n_dim = n_dim,
                                               method = 'CG',
                                               control = list(maxit = 500))}
  )
  ## Check convergence
  if(!all(sapply(fit_list, function(ll) ll$convergence) == 0)){
    warning('Not all optimizations converged! - Check individual fits')
    return(fit_list)
  }
  ## Check that if all runs converged, that they all converged to the same solution
  errfun_value <- sapply(fit_list, function(ll) ll$value)
  if(!all(abs(errfun_value - mean(errfun_value) < .5))){warning('Not all runs converged to the same solution. Check individual fits.')}
  
  return(list(best_fit = fit_list[errfun_value == min(errfun_value)][[1]],
              all_fits = fit_list))
}

align_mds_with_original <- function(mds_coords, # Nxd matrix containing all ag and ab coordinates
                                    original_coords #Nxd matrix containing all ag and ab coordinates
                                    ){
  stopifnot(all(dim(mds_coords) == dim(original_coords)))
  
  ## Shift so that the ag1 coordinate is at the origin for both maps
  shifted_mds = apply(mds_coords, 2, function(vv) vv-vv[1])
  shifted_original = apply(original_coords, 2, function(vv) vv-vv[1])
  
  ## Get the angles between corresponding coordinates in each map
  get_cosine <- function(v1, v2){
    dot_prod = sum(v1*v2) # dot product
    mag_v1 = sqrt(sum(v1^2))
    mag_v2 = sqrt(sum(v2^2))
    cos_theta = dot_prod/(mag_v1*mag_v2)
    theta = acos(cos_theta)
    theta
  }
  
  angle_original_to_mds = sapply(2:nrow(shifted_original), function(ii) get_cosine(shifted_original[ii,], shifted_mds[ii,]))
  ## If all the angles are equal, then the mds map is a rotation of the original
  ## If the angles are different, then the mds map is a reflection and a rotation
  if(!all( abs(angle_original_to_mds - mean(angle_original_to_mds)) < .1 )){# If it's not just a rotation
    ## Reflect across a vertical axis
    shifted_mds[,1] = shifted_mds[,1] - min(shifted_mds[,1]) # Shift so that all points are non-negative on the x-axis
    reflected_mds = t(matrix(c(-1, 0, 0, 1), 2, 2) %*% t(shifted_mds)) # Reflect across the y-axis
    shifted_mds = apply(reflected_mds, 2, function(vv) vv-vv[1]) # Re-shift the ag 1 point to the origin
    angle_original_to_mds = sapply(2:nrow(shifted_original), function(ii) get_cosine(shifted_original[ii,], shifted_mds[ii,])) ## Recalculate angles and check that we're now dealing with a rotation of the original map
    stopifnot(all( abs(angle_original_to_mds - mean(angle_original_to_mds)) < .1 )) ## Break if the reflected map is still not a simple rotation of the original map
  }
  
  ## Rotate to align with the original
  theta1 =  mean(angle_original_to_mds) # Get the two possible rotation angles
  theta2 =  2*pi - mean(angle_original_to_mds)
  V1 = matrix(c(cos(theta1), sin(theta1), -sin(theta1), cos(theta1)), nrow = 2, ncol = 2) # Get the rotation matrix
  rotated_mds1 = t(V1 %*% t(shifted_mds)) # Rotate by theta
  V2 = matrix(c(cos(theta2), sin(theta2), -sin(theta2), cos(theta2)), nrow = 2, ncol = 2) # Get the rotation matrix
  rotated_mds2 = t(V2 %*% t(shifted_mds)) # Rotate by theta
  ## Choose the better fit
  error1 = sum((rotated_mds1-shifted_original)^2)
  error2 = sum((rotated_mds2-shifted_original)^2)
  rotated_mds = if(error1>error2) rotated_mds2 else  rotated_mds1
  
  ## Shift to align with the original and return
  sapply(1:ncol(rotated_mds), function(cc) rotated_mds[,cc] - rotated_mds[1,cc] +  original_coords[1,cc]) 
  

  # plot(shifted_original, xlim = c(-20, 20), ylim = c(-20, 20), pch = 16)
  # points(shifted_mds, col = 'red')
  # points(reflected_mds, col = 'blue')
  # points(rotated_mds1, col = 'green')
  # points(rotated_mds2, col = 'magenta')
}