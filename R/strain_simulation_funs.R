generate_ag_coords_randomwalk <- function(n_antigens, n_epitopes){
  lapply(1:n_epitopes, FUN = function(xx){
    cvec = numeric(n_antigens)
    for(ii in 2:n_antigens){
      cvec[ii] = cvec[ii-1] + rnorm(1, 2, 2)
    }
    return(cvec)
  }
  ) %>%
    unlist()
}

generate_ag_coords_random <- function(n_antigens, n_epitopes){
  lapply(1:n_epitopes, FUN = function(xx){
    cvec = runif(n_antigens, 0, 10)
  }
  ) %>%
    unlist()
}



## Generate a population of Abs specific to a set of native coordinates
generate_gaussian_repertoire <- function(native_epitope_coords, ## A data frame with columns for epitope, strain id, and each coordinate
                                         n_epitopes, ## n Epitopes in simulation
                                         n_dim,
                                         n_ab = 1000, ## Total number of Abs to draw
                                         rel_immuno, ## Vector giving the relative immunodominance of each epitope. Will be normalized. The fraction of Abs specific to epitope i will be proportional to the ith entry.
                                         sigma = .1 ## sd of Ab positions around native_epitope_Coords
){
  ## Figure out how many Abs to draw for each epitope
  draw_this_many = floor(n_ab*rel_immuno/sum(rel_immuno))
  remainder = n_ab - sum(draw_this_many)
  cat(sprintf('Drawing %i Abs per individual. %i were requested. Remainder is %i.\n', sum(draw_this_many), n_ab, remainder))
  stopifnot(length(draw_this_many) == nrow(native_epitope_coords))
  stopifnot(n_dim <= sum(grepl(pattern = 'c\\d?.+', x = names(native_epitope_coords))))
  ## Drop unused antigen_coords columns
  native_epitope_coords <- select(native_epitope_coords, 'epitope', 'antigen', 'kind', paste0('c', 1:n_dim))
  
  library(foreach)
  foreach(ee = native_epitope_coords$epitope, 
          aa = native_epitope_coords$antigen, 
          this.row = 1:nrow(native_epitope_coords), 
          nn = draw_this_many) %do% {
            coord.means = native_epitope_coords %>%
              select(matches('c\\d?.+')) %>%
              slice(this.row) %>%
              as.vector
            tibble(epitope = ee,
                   antigen = aa,
                   kind = 'antibody')%>%
              bind_cols(
              map_dfc(.x = coord.means, .f = ~{rnorm(nn, .x, sd = sigma)})
              )
          } %>%
    bind_rows()
}


## Extract inputs to calculate titer 
extract_titer_inputs <- function(antigen_id, # To which antigen are we calculating titer 
                                 serum_id, # Which serum are we using?
                                 merged_df # data frame of antigen and antibody locations
                                 ){
  
  ## Extract ab_position_list
  these_serum_abs <- merged_df %>% # Subset the master data frame
    filter(serum == serum_id) %>%
    select(epitope, matches('c\\d_Ab'))
  
  ab_position_list = lapply(unique(these_serum_abs$epitope), ## Reformat into a list of matrices, one for each epitope
                            function(this.epitope){
                              these_serum_abs %>% 
                                filter(epitope == this.epitope) %>%
                                select(starts_with('c')) %>%
                                as.matrix
                            }) %>%
    set_names(paste0('E', 1:length(.)))
  
  
  ## Extract the ag coordinates as a matrix
  this_ag = merged_df %>%
    filter(antigen == antigen_id) %>%
    select(epitope, matches('c\\d_Ag')) %>%
    distinct() %>%
    select(-epitope) %>%
    as.matrix()
  ag_list = lapply(1:nrow(this_ag), function(ii){
    this_ag[ii,]
  })
  
  ## Return
  list(ab_position_list = ab_position_list, 
       ag_list = ag_list)
}


## Solve numerically for titer
## Inputs formatted by extract_titer_inputs
solve_for_titer_multi_epitope <- function(ab_position_list, # A list of E matrices whose rows represent the coordinates of Abs specific to epitope E. Each matrix should be of dimension  [n_E, d].
ag_list, # A list of length E giving the coordinate vectors of each epitope within the antigen
alpha, # Ab potency (equal and fixed across all Abs)
r # Maximum log2 titer
){ # Location of each epitope in the antigen. Vector of length E.
  
  # cat(sprintf('ag list length is %i; ab list lenght is %i', length(ag_list), length(ab_position_list)))
  stopifnot(length(ag_list)==length(ab_position_list))
  # valid_epitopes = (1:length(relative_concentrations))[relative_concentrations>0]
  # if(length(valid_epitopes) < 1) relative_concentrations = rep(1, length(ab_position_list))
  # stopifnot(length(valid_epitopes) == length(ab_position_list))
  
  
  titer_fun <- function(z, ab_position_list, x_ag, alpha, r, concentration = 1, mu = 1/2){ 
    
    # Each eptiope is a factor in Eq 13. This function calculates those factors
    one_factor <- function(ab_positions, # matrix of Ab positions to one epitope
                           x_ag = x_ag, # vector of epitope's position
                           z, # titer (we'll solve for this later)
                           alpha, # potency
                           r,
                           concentration = 1,
                           mu
    ){
      # cat(sprintf(' ID is %2.1f', this_immunodominance))
      ab_distances = apply(ab_positions, 1, FUN = function(v1){
        get_euclidean_distance(v1, v2 = x_ag)
      })
      
      ab_affinities = 2^(r-ab_distances) # Eq 6
      scaled_affinities = ab_affinities*concentration
      this_factor = (1+alpha/z*sum(scaled_affinities))/(1+1/z*sum(scaled_affinities)) # Eq 13
      return(this_factor)
      
    }
    ## This function calculates each factor and takes the product to solve for titer
    factors<-mapply(FUN = one_factor, ab_positions = ab_position_list, x_ag = ag_list, z = z, alpha = alpha, r = r)  
    # cat(print('\n'))
    # cat(print('factors are\n'))
    # print(factors)
    prod(factors) - mu
  }
  
  ## Solve for titer and return
  numerical_soln <- uniroot(f = titer_fun, 
                            interval = c(.0000000001, 10^10), 
                            ab_position_list = ab_position_list, 
                            x_ag = ag_list, 
                            alpha = alpha, 
                            r = r,
                            concentration = 1,
                            mu = 1/2,
                            tol =  .Machine$double.eps^.5)       
  
  ## Verify that the solution returns 0
  stopifnot(zeroish(titer_fun(numerical_soln$root, ab_position_list, x_ag, alpha, r, relative_concentrations)))
  
  return(numerical_soln$root)
}







## Generate serum panel, and titer map, for a given set of antigen coords and an immunodominance scheme
generate_ferret_inputs <- function(antigen_coords, # Data frame of ag coords 
                            relative_immunodominance, # vector of relative concentrations (immunodominance) of Abs to each epitope.
                            n_epitopes, # n epitiopes
                            n_antigens,
                            n_dim,
                            n_abs_per_serum = 500,
                            sigma = .1, # sd of abs around native coord
                            immunodominance_flag,
                            outdir = 'inputs',
                            alpha = .25, r = 7
                            ){ # n antigens
  stopifnot(n_dim <= sum(grepl(pattern = 'c\\d?.+', x = names(antigen_coords))))
  ## Drop unused antigen_coords columns
  antigen_coords <- select(antigen_coords, 'epitope', 'antigen', 'kind', paste0('c', 1:n_dim))
  if(sum(relative_immunodominance) != n_epitopes){warning('Relative concentrations dont sum to n_epitopes. This will rescale the total Ab concentration.')}
  stopifnot(length(unique(antigen_coords$antigen)) == n_antigens)
  stopifnot(length(unique(antigen_coords$epitope)) == n_epitopes)
  
  ## Generate ferret repertoire
  ferret_repertoires <- lapply(1:n_antigens, function(this_ag){
    
    generate_gaussian_repertoire(native_epitope_coords = antigen_coords %>% filter(antigen == this_ag), 
                               n_epitopes = n_epitopes,
                               n_dim = n_dim,
                               n_ab = n_abs_per_serum, ## n antibodies per ferret
                               rel_immuno = relative_immunodominance,
                               sigma = sigma)
  }) %>%
    bind_rows(.id = 'serum')
  
  ## Merge the antigen coordinates with the Ab coordinates to calculate titers
  merged_df <- merge(antigen_coords %>%
                       select(epitope, antigen, matches('c\\d?.+')),
                     ferret_repertoires %>%
                       select(serum, epitope, antigen, matches('c\\d?.+')),
                     by = c('epitope', 'antigen'), 
                     suffixes = c('_Ag', '_Ab'))
  
  ## Calculate the titer panel
  ## Wrapper to calculate the titer for a given serum and antigen
  
  titer_wrapper <- function(this.serum, 
                            this.antigen){
    titer = with(extract_titer_inputs(antigen_id = this.antigen, 
                                      serum_id = this.serum, 
                                      merged_df = merged_df),
                 solve_for_titer_multi_epitope(ab_position_list = ab_position_list, 
                                               ag_list = ag_list, 
                                               alpha = alpha,
                                               r = r)
    )
    c('serum' = this.serum, 'antigen' = this.antigen, 'titer' = titer)
  }
  
  titer_map <- with(expand.grid(antigen = 1:n_antigens, serum = 1:n_antigens),
                    mapply(FUN = titer_wrapper, 
                           this.serum = serum,
                           this.antigen = antigen)) %>% 
    t() %>%
    as_tibble() %>%
    set_names(c('serum', 'antigen', 'titer')) %>%
    mutate(logtiter = log2(titer/10)) %>%
    group_by(serum) %>%
    mutate(serum_potency = max(logtiter)) %>%
    ungroup() %>% group_by(antigen) %>%
    mutate(antigen_avidity = max(logtiter)) %>%
    ungroup() %>%
    mutate(titer_distance = (serum_potency+antigen_avidity)/2 - logtiter)
  
  if(!dir.exists(outdir)) dir.create(outdir)
  outfilename = sprintf('%s/%sD_%s_immunodominance_inputs.rds', outdir, n_dim, immunodominance_flag)
  write_rds(
  list(ag_ab_coords = merged_df,
              titer_map = titer_map),
  file = outfilename
  )
  cat(sprintf('saved outputs in %s', outfilename))
}
  


plot_antigen_coords <- function(antigen_coords){
  antigen_coords %>%
    ggplot() +
    geom_line(aes(x=c1, y = c2, group = epitope), lty = 2, lwd = .1)+
    geom_point(aes(x = c1, y = c2, shape = antigen, color = epitope))+
    facet_wrap(.~epitope, labeller = label_both)
}


plot_Ab_Ag_map <- function(ab_ag_df,
                           antigen_coords){
  antisera_to_strain = unique(ab_ag_df$antigen)
  ab_ag_df %>%
    ggplot() +
    geom_text(aes(x = c1, y = c2, color = antigen, label = epitope), data = antigen_coords, show.legend = F) + ## All epitope coords
    geom_point(aes(x = c1_Ab, y = c2_Ab, color = antigen), pch = 3, alpha = .4) + # Ab coords
    geom_text(aes(x = c1, y = c2, label = epitope), data = antigen_coords %>% filter(antigen == antisera_to_strain)) + # Highlight target antigen coords in black text and a circle
   # geom_point(aes(x = c1, y = c2), pch = 1, data = antigen_coords %>% filter(antigen == antisera_to_strain)) +
    guides(color=guide_legend(title="Abs to antigen"))+
    xlab('c1')+
    ylab('c2')+
    ggtitle(sprintf('Antiserum to strain %s', antisera_to_strain))
}



plot_inferred_original_map <- function(antigen_coords,
                                       stan_fit,
                                       flip_these_chains_over_x = NULL){
  reformatted_df <- antigen_coords %>%
    mutate(kind = paste0('E', epitope)) %>%
    rename(allele = antigen) %>%
    select(allele, kind, starts_with('c')) %>%
    mutate(kind = factor(kind, levels = c('antigen', 'serum', 'fixed_ag1', 'E1', 'E2', 'E3'),
                         labels = c('antigen (estimated)', 'serum (estimated)', 'antigen 1 position (fixed)', 'epitope 1 (truth)', 'epitope 2 (truth)', 'epitope 3 (truth)')), ordered = T)
  
  epitope_strain_map <- extract_summary_coords(stan_fit) %>%
    mutate(across(starts_with('c2'), ~ifelse(chain %in% flip_these_chains_over_x, -.x, .x))) %>%
    mutate(allele = as.factor(id)) %>%
    mutate(kind = factor(kind, levels = c('antigen', 'serum', 'fixed_ag1', 'E1', 'E2', 'E3'),
                         labels = c('antigen (estimated)', 'serum (estimated)', 'antigen 1 position (fixed)', 'epitope 1 (truth)', 'epitope 2 (truth)', 'epitope 3 (truth)')), ordered = T) %>%
    ggplot() +
    geom_point(aes(x = c1, y = c2, shape = kind), data = reformatted_df) +
    geom_polygon(aes(x = c1, y = c2, fill = allele), data = reformatted_df, alpha = .2)+
    geom_point(aes(x = c1, y = c2, shape = kind, color = allele)) +
    facet_grid(.~chain, labeller = label_both) +
    scale_shape_manual(values = c(4, 14, 21, 22, 23, 3))
  
  epitope_strain_map
}





generate_OAS_serum <- function(existing_Ab_coords, ## ab_ag_coords data frame with columns: `epitope`, c1, c2, ... cn
                               strain_coords, ## data frame with columns: epitope, c1, c2, ... cn
                               total_ab = 1000, ## Total Abs to return
                               preserve_epitope_immunodominance = FALSE, # if TRUE, return the fraction of Abs specific to each epitope as in the original serum. if FALSE, let the relative affinities of existing Abs to their cognate epitopes dictate immunodominance.
                               this_antigen,
                               sigma = .1,
                               r=7){
  
  stopifnot('epitope' %in% colnames(existing_Ab_coords))
  stopifnot('epitope' %in% colnames(strain_coords))
  stopifnot(any(grepl('c\\d?', x = colnames(strain_coords))))
  stopifnot(any(grepl('c\\d?_Ab', x = colnames(existing_Ab_coords))))
  existing_Ab_coords <- existing_Ab_coords %>% select(epitope, matches('c\\d?_Ab'))
  strain_coords <- strain_coords %>% select(epitope, matches('c\\d?'))
  n_dim = sum(grepl(pattern = 'c\\d?', colnames(strain_coords)))
  stopifnot(n_dim == sum(grepl(pattern = 'c\\d?', colnames(existing_Ab_coords))))
  existing_Ab_coords <- set_names(existing_Ab_coords, c('epitope', gsub(pattern = '(c\\d?)_Ab', replacement = '\\1', x = colnames(existing_Ab_coords)[-1])))
  
  if(preserve_epitope_immunodominance){
    ## Calculate the number of Abs per epitope in the original serum
    n_epitopes = length(unique(existing_Ab_coords$epitope))
    epitope_counts = existing_Ab_coords %>% group_by(epitope) %>% summarise(n = n()) %>% arrange(epitope) %>% pull(n)
    nAbs_per_epitope = floor(total_ab*(epitope_counts/sum(epitope_counts)))
    remainder = total_ab-sum(nAbs_per_epitope)
    nAbs_per_epitope = nAbs_per_epitope + rmultinom(n = 1, size = remainder, prob = rep(1, n_epitopes)) %>% t()
    stopifnot(sum(nAbs_per_epitope) == total_ab)
    
    library(foreach)
    ## For each epitope...
    new_Abs = foreach(this_epitope = 1:n_epitopes,
        this_nAb = nAbs_per_epitope) %do% {
          ## START MODIFYING HERE
        candidate_locations = rbind(strain_coords %>% filter(epitope == this_epitope) %>% select(matches('c\\d?')),
                                    existing_Ab_coords %>% filter(epitope == this_epitope) %>% select(matches('c\\d?')))
        distances_to_candidate_locations = apply(candidate_locations, 1, function(vv){get_euclidean_distance(vv, candidate_locations[1,])})
        affinities = 2^(r-distances_to_candidate_locations)
        ## Select locations
       #new_location_weights = rmultinom(n = 1, size = nrow(candidate_locations), prob = affinities/sum(affinities)) 
        new_location_indices =  sample(x = 1:nrow(candidate_locations), size = this_nAb, replace = T, prob = affinities) 
        new_location_coords = as.data.frame(candidate_locations)[new_location_indices, ]
        stopifnot(nrow(new_location_coords)==this_nAb)
        
        ## And for each location...
        foreach(this.row = 1:nrow(new_location_coords)) %do% {
                  coord.means = new_location_coords[this.row, ]
                  tibble(epitope = this_epitope,
                         antigen = this_antigen,
                         kind = 'antibody')%>%
                    bind_cols(
                      map_dfc(.x = coord.means, .f = ~{rnorm(1, .x, sd = sigma)})
                    )
                } %>%
          bind_rows()
        } %>%
      bind_rows()
  }else{
    candidate_locations = rbind(strain_coords,
                                existing_Ab_coords)
    candidate_epitopes = candidate_locations$epitope
    candidate_locations = select(candidate_locations, -epitope)
    distances_to_candidate_locations = apply(candidate_locations, 1, function(vv){get_euclidean_distance(vv, candidate_locations[1,])})
    affinities = 2^(r-distances_to_candidate_locations)
    new_location_indices =  sample(x = 1:nrow(candidate_locations), size = total_ab, replace = T, prob = affinities) 
    new_location_coords = as.data.frame(candidate_locations)[new_location_indices, ]
    new_epitopes = candidate_epitopes[new_location_indices]
    stopifnot(nrow(new_location_coords)==total_ab)
    ## And for each location...
    library(foreach)
    new_Abs = foreach(this.row = 1:nrow(new_location_coords),
                      this_epitope = new_epitopes) %do% {
      coord.means = new_location_coords[this.row, ]
      tibble(epitope = this_epitope,
             antigen = this_antigen,
             kind = 'antibody')%>%
        bind_cols(
          map_dfc(.x = coord.means, .f = ~{rnorm(1, .x, sd = sigma)})
        )
    } %>%
      bind_rows()
  }
  new_Abs
}



generate_one_serum <- function(native_epitope_coords, ## A data frame with columns for epitope, strain id, and each coordinate
                                n_epitopes,
                                total_ab = 1000,
                                epitope_immunodominance, # a vector of the relative immunodominances of each epitope. Will be normalized.
                                strain_dominance, # a vector of the relative dominance of each strain in the repertoire
                                sigma = .1, ## sd of Ab positions around native_epitope_Coords)
                                n_dim
){
  
  n_ab_per_strain <- floor(total_ab*strain_dominance/sum(strain_dominance))
  strains <- as.numeric(unique(native_epitope_coords$antigen))
  
  ## For each strain, draw Abs
  foreach(this.strain = strains, this.n_ab = n_ab_per_strain) %do% {
    generate_gaussian_repertoire(native_epitope_coords = filter(native_epitope_coords, antigen == this.strain),
                                 n_epitopes = n_epitopes, 
                                 n_dim,
                                 n_ab = this.n_ab, 
                                 rel_immuno = epitope_immunodominance, 
                                 sigma = sigma) 
  } %>%
    bind_rows() %>%
    merge(native_epitope_coords,
          by = c('epitope', 'antigen'), suffixes = c('_Ab', '_Ag')) %>%
    select(-starts_with('kind'))
}  





get_weighted_centroids <- function(antigen_coords,
                                   immunodominance_weights){
  immunodominance_weights = immunodominance_weights/sum(immunodominance_weights)
  antigen_coords %>%
    group_by(antigen) %>%
    summarise(epitope = NA,
              antigen = unique(antigen), 
              kind = 'antigen',
              c1 = sum(c1*immunodominance_weights),
              c2 = sum(c2*immunodominance_weights))
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
