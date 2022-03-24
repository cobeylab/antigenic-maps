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
                                          relative_concentrations = NULL, # An optional vector that scales the immunodominance of each epitope. (Assume in eq. 18 that concentration is proportional to the cognate epitope's immunodominance)
ag_list, # A list of length E giving the coordinate vectors of each epitope within the antigen
alpha, # Ab potency (equal and fixed across all Abs)
r # Maximum log2 titer
){ # Location of each epitope in the antigen. Vector of length E.
  
  # cat(sprintf('ag list length is %i; ab list lenght is %i', length(ag_list), length(ab_position_list)))
  stopifnot(length(ag_list)==length(ab_position_list))
  valid_epitopes = (1:length(relative_concentrations))[relative_concentrations>0]
  if(length(valid_epitopes) < 1) relative_concentrations = rep(1, length(ab_position_list))
  stopifnot(length(valid_epitopes) == length(ab_position_list))
  
  
  titer_fun <- function(z, ab_position_list, x_ag, alpha, r, relative_concentrations){ 
    
    # Each eptiope is a factor in Eq 13. This function calculates those factors
    one_factor <- function(ab_positions, # matrix of Ab positions to one epitope
                           x_ag = x_ag, # vector of epitope's position
                           z, # titer (we'll solve for this later)
                           alpha, # potency
                           r,  # max titer
                           this_immunodominance # the immunodominance (relative concentration) of this epitope
    ){
      # cat(sprintf(' ID is %2.1f', this_immunodominance))
      ab_distances = apply(ab_positions, 1, FUN = function(v1){
        get_euclidean_distance(v1, v2 = x_ag)
      })
      
      ab_affinities = 2^(r-ab_distances) # Eq 6
      scaled_affinities = ab_affinities*this_immunodominance
      this_factor = (1+alpha/z*sum(scaled_affinities))/(1+1/z*sum(scaled_affinities)) # Eq 13
      return(this_factor)
      
    }
    ## This function calculates each factor and takes the product to solve for titer
    factors<-mapply(FUN = one_factor, ab_positions = ab_position_list, x_ag = ag_list, this_immunodominance = relative_concentrations, z = z, alpha = alpha, r = r)  
    # cat(print('\n'))
    # cat(print('factors are\n'))
    # print(factors)
    prod(factors) - 1/2
  }
  
  ## Solve for titer and return
  numerical_soln <- uniroot(f = titer_fun, 
                            interval = c(.0000000001, 10^10), 
                            ab_position_list = ab_position_list, 
                            x_ag = ag_list, 
                            alpha = alpha, 
                            r = r,
                            relative_concentrations = relative_concentrations,
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
                            sigma = 1, # sd of abs around native coord
                            immunodominance_flag,
                            outdir = 'inputs'
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
                                               relative_concentrations = relative_immunodominance, 
                                               ag_list = ag_list, 
                                               alpha = .25,
                                               r = 7)
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



generate_one_serum <- function(native_epitope_coords, ## A data frame with columns for epitope, strain id, and each coordinate
                                n_epitopes,
                                total_ab = 1000,
                                epitope_immunodominance, # a vector of the relative immunodominances of each epitope. Will be normalized.
                                strain_dominance, # a vector of the relative dominance of each strain in the repertoire
                                sigma = .1 ## sd of Ab positions around native_epitope_Coords)
){
  
  n_ab_per_strain <- floor(total_ab*strain_dominance/sum(strain_dominance))
  strains <- as.numeric(unique(native_epitope_coords$antigen))
  
  ## For each strain, draw Abs
  foreach(this.strain = strains, this.n_ab = n_ab_per_strain) %do% {
    generate_gaussian_repertoire(native_epitope_coords = filter(native_epitope_coords, antigen == this.strain),
                                 n_epitopes = n_epitopes, 
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


