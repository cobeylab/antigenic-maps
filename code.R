## Functions to simulate maps and test MDS

library(tidyverse)
filter <- dplyr::filter
select <- dplyr::select
extract <- tidyr::extract

zeroish <- function(x){
  abs(x) <= .Machine$double.eps^.5
}

ndim <- function(xx){
  ldim <- length(dim(xx)) ## Unlike dim, will print 0 if a vector
  if(ldim >= 2){
    return(ldim)
  }else{
    ifelse(length(dim) > 1, 1, 'scalar')
  }
}

equal_ish <- function(xx, yy = 0){
  abs(xx-yy) <= 1e-6
}

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
  id_names = c(
    if(n_antigens < 1) NULL else 1:n_antigens,
    if(n_sera < 1) NULL else 1:n_sera
  )
  tibble(id = id_names,
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
  
  ag_coords <- runif(n_antigens*n_dim, -10, 10)
  names(ag_coords) = paste0('ag', rep(1:n_antigens, n_dim), 'c', rep(1:n_dim, each = n_antigens))
  
  ab_coords <- runif(n_sera*n_dim, -10, 10)
  names(ab_coords) = paste0('ab', rep(1:n_sera, n_dim), 'c', rep(1:n_dim, each = n_sera))
  
  c(ag_coords, ab_coords)
}


fit_MDS_least_squares <- function(
  observed_distances, # vector of delay for each individual in the data
  n_antigens,
  n_sera,
  n_dim,
  n_cores = (detectCores() - 1)
) {
  
  library(doParallel)
  stopifnot(nrow(observed_distances) == n_antigens)
  stopifnot(ncol(observed_distances) == n_sera)
  
  # Define the loss (error) function
  loss_function <- function(pars, observed_distances, n_antigens, n_sera, n_dim){
    stopifnot(length(pars) == (n_antigens+n_sera)*n_dim)
    
    estimated_ag_coords <- matrix(pars[1:(n_antigens*n_dim)], n_antigens, n_dim)
    estimated_ab_coords <- matrix(pars[(n_antigens*n_dim+1):length(pars)], n_sera, n_dim)
    
    estimated_distances <- get_ab_ag_distances(estimated_ab_coords, estimated_ag_coords)
    sum((observed_distances - estimated_distances)^2)
  }
  
  ## Minimize the loss function. Repeat the fit 10 times from different initial values to verify convergence
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  fit_list <- foreach(ii=1:10,
                      .export = c('generate_initial_guess', 
                                  'get_ab_ag_distances',
                                  'get_euclidean_distance')) %dopar% {
    optim(par = generate_initial_guess(n_antigens, n_sera, n_dim),
          fn = loss_function,
          observed_distances = observed_distances,
          n_antigens = n_antigens, 
          n_sera = n_sera,
          n_dim = n_dim,
          method = 'CG',
          control = list(maxit = 500))
    
  }

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



## Get the angles between corresponding coordinates in each map
get_theta <- function(v1, v2){
  dot_prod = sum(v1*v2) # dot product
  mag_v1 = sqrt(sum(v1^2))
  mag_v2 = sqrt(sum(v2^2))
  cos_theta = dot_prod/(mag_v1*mag_v2)
  theta = acos(cos_theta)
  theta
}

shift_to_origin <- function(coord_mat,   ##  a matrix of dimension N x d, where d is the number of dimensions, and N is the total number of antigens and sera in the analysis.
                            origin_row = 1 ## which row index should be shifted to the origin?
                            ){
  ## For each column, shift so that the reference row has coordinate 0
  apply(coord_mat, 2, function(vv) vv-vv[origin_row])
}

rotate_by_theta_2d <- function(theta, coords){
  stopifnot(ncol(coords) == 2)
  V = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2) # Get the rotation matrix
  t(V %*% t(coords)) ## rotate and return
}

rotate_about_y_axis_3d <- function(theta, coords){
  stopifnot(ncol(coords) == 3)
  ## Formula from equation 20 of http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf
  V = matrix(
    c(cos(theta), 0, -sin(theta),
      0,          1, 0,
      sin(theta), 0, cos(theta)),
    nrow = 3, ncol = 3, byrow = T) # Get the rotation matrix
  # cat('\nRotation matrix\n')
  # print(V)
  t(V %*% t(coords)) ## rotate and return
}

rotate_about_z_axis_3d <- function(theta, coords){
  stopifnot(ncol(coords) == 3)
  ## Formula from equation 20 of http://scipp.ucsc.edu/~haber/ph216/rotation_12.pdf
  V = matrix(
    c(cos(theta), -sin(theta), 0,
      sin(theta), cos(theta),  0,
      0,           0,          1),
    nrow = 3, ncol = 3, byrow = T) # Get the rotation matrix
  # cat('\nRotation matrix\n')
  # print(V)
  t(V %*% t(coords)) ## rotate and return
}






align_mds_with_original <- function(mds_coords, # Nxd matrix containing all ag and ab coordinates
                                    original_coords #Nxd matrix containing all ag and ab coordinates
                                    ){
  stopifnot(all(dim(mds_coords) == dim(original_coords)))
  
  ## Shift so that the ag1 coordinate is at the origin for both maps
  shifted_mds = shift_to_origin(mds_coords)
  shifted_original = shift_to_origin(original_coords)
  
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
  # theta1 =  mean(angle_original_to_mds) # Get the two possible rotation angles
  # theta2 =  2*pi - mean(angle_original_to_mds)
  # V1 = matrix(c(cos(theta1), sin(theta1), -sin(theta1), cos(theta1)), nrow = 2, ncol = 2) # Get the rotation matrix
  # rotated_mds1 = t(V1 %*% t(shifted_mds)) # Rotate by theta
  # V2 = matrix(c(cos(theta2), sin(theta2), -sin(theta2), cos(theta2)), nrow = 2, ncol = 2) # Get the rotation matrix
  # rotated_mds2 = t(V2 %*% t(shifted_mds)) # Rotate by theta
  # ## Choose the better fit
  # error1 = sum((rotated_mds1-shifted_original)^2)
  # error2 = sum((rotated_mds2-shifted_original)^2)
  # rotated_mds = if(error1>error2) rotated_mds2 else  rotated_mds1
  
  theta =  2*pi - mean(angle_original_to_mds)
  V = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2) # Get the rotation matrix
  rotated_mds = t(V %*% t(shifted_mds)) # Rotate by theta
  
  ## Shift to align with the original and return
  sapply(1:ncol(rotated_mds), function(cc) rotated_mds[,cc] - rotated_mds[1,cc] +  original_coords[1,cc]) 
  

  # plot(shifted_original, xlim = c(-20, 20), ylim = c(-20, 20), pch = 16)
  # points(shifted_mds, col = 'red')
  # points(reflected_mds, col = 'blue')
  # points(rotated_mds1, col = 'green')
  # points(rotated_mds2, col = 'magenta')
}

## Function to extact individual fits from the list of outputs
get_fits_df <- function(fit_obj){
  tibble(coordinate = fit_obj$par,
         label = names(fit_obj$par)) %>%
    tidyr::extract(label, into = c('kind','id', 'dim'), regex = '(a\\w)(\\d)c(\\d)', convert = T) %>%
    pivot_wider(names_from = dim, values_from = coordinate, names_prefix = 'c') 
}



## This function takes in a matrix of 2d or 3d coordinates, and shifts/rotates them so that:
##   1. ag1 coordinates are at the origin
##   2. ag2 coordinates fall on the x-axis
standardize_coordinates <- function(coord_mat, 
                                    ag1_row = 1,
                                    ag2_row = 2,
){
  n_dim = ncol(coord_mat) 
  stopifnot(nrow(coord_mat) >= n_dim)
  
  ## Shift the whole matrix so that ag1 is at the origin.
  shifted = shift_to_origin(coord_mat, origin_row = ag1_row)
  # cat('shifted cords\n')
  # print(shifted)
  
  if(n_dim == 3){
    ## The coordinates of the vector (ag1, ag2) are now the ag2 coordinates, because all ag1 coordinates are 0.
    ## We want to rotate the map so that ag2 falls on axis 1
    x_z_projection_vec = c(shifted[2, 1], 0, shifted[2,3])
    theta = get_theta(v1 = c(1,0,0), v2 = x_z_projection_vec)
    theta = ifelse(shifted[2,3]>0, 2*pi-theta, theta) ## Rotate in the opposite direction if z is positive
    rotated_matrix = rotate_about_y_axis_3d(theta = theta, coords = shifted)
    
    ## Now the ag2 vector is in the x-y plane.
    ## Rotate again about the z-axis to align ag2 with the x axis, and return
    theta2 = get_theta(c(rotated_matrix[2,1],0,0), rotated_matrix[2,])
    final_output = rotate_about_z_axis_3d(theta = theta2, 
                                          coords = rotated_matrix)
    
    # cat('\nafter rotation 1:\n')
    # print(rotated_matrix[2,])
    # cat('\nafter rotation 2:\n')
    stopifnot(equal_ish(final_output[2,2], 0) & equal_ish(final_output[2,3],0))
    
  }else if(n_dim == 2){
    ## Rotate about the origin to align with x and return
    theta = get_theta(shifted[ag2_row,], c(1,0))
    theta = ifelse(shifted[ag2_row,2]>0, 2*pi-theta, theta)
    # cat('theta:\n')
    # print(theta)
    # cat('original:\n')
    # print(shifted[ag2_row,])
    final_output = rotate_by_theta_2d(theta = theta, 
                                      coords = shifted)
    stopifnot(equal_ish(final_output[2,2],0))
    
  }else{
    stop('standardization only implemented for 2d or 3d maps')
  }
  return(final_output)
}
