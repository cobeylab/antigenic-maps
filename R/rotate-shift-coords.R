
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
  
  theta =  2*pi - mean(angle_original_to_mds)
  V = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2) # Get the rotation matrix
  rotated_mds = t(V %*% t(shifted_mds)) # Rotate by theta
  
  ## Shift to align with the original and return
  sapply(1:ncol(rotated_mds), function(cc) rotated_mds[,cc] - rotated_mds[1,cc] +  original_coords[1,cc]) 
  
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
                                    ag2_row = 2
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
    
    stopifnot(equal_ish(final_output[2,2], 0) & equal_ish(final_output[2,3],0))
    
  }else if(n_dim == 2){
    ## Rotate about the origin to align with x and return
    theta = get_theta(shifted[ag2_row,], c(1,0))
    theta = ifelse(shifted[ag2_row,2]>0, 2*pi-theta, theta)
    final_output = rotate_by_theta_2d(theta = theta, 
                                      coords = shifted)
    stopifnot(equal_ish(final_output[2,2],0))
    
  }else{
    stop('standardization only implemented for 2d or 3d maps')
  }
  return(final_output)
}

standardize_coordinate_df<-function(coord_df,
                                    ag1_row,
                                    ag2_row){
  # Get standarzied coords
  standardized_matrix = standardize_coordinates(coord_mat = coord_df %>% select(matches('c\\d')) %>% as.matrix())
  
  # Replace the original coordinates in the data frame and return
  colnames(standardized_matrix) = colnames(coord_df %>% select(matches('c\\d')))
  coord_df %>%
    select(-matches('c\\d')) %>%
    bind_cols(as_tibble(standardized_matrix))
}
