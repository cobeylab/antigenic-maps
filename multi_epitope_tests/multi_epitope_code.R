
reset <- function(){
  rm(list = ls())
  library(tidyr)
  filter <- dplyr::filter
  extract <- tidyr::extract
  select <- dplyr::select
  source('../code.R')
}


generate_ag_coords_randomwalk <- function(n_ag, n_epitopes){
  lapply(1:n_epitopes, FUN = function(xx){
    cvec = numeric(n_ag)
    for(ii in 2:n_ag){
      cvec[ii] = cvec[ii-1] + rnorm(1, 2, 2)
    }
    return(cvec)
  }
  ) %>%
    unlist()
}

generate_ag_coords_random <- function(n_ag, n_epitopes){
  lapply(1:n_epitopes, FUN = function(xx){
    cvec = runif(n_ag, 0, 10)
  }
  ) %>%
    unlist()
}



## Generate a population of Abs specific to a set of native coordinates
generate_ferret_repertoire <- function(native_epitope_coords, ## A data frame with columns for epitope, strain id, and each coordinate
                                       n_epitopes,
                                       n_ab = 99,
                                       rel_immuno,
                                       sigma = 1 ## sd of Ab positions around native_epitope_Coords
                                       ){
  ## Figure out how many Abs to draw for each epitope
  draw_this_many = floor(n_ab*rel_immuno/sum(rel_immuno))
  remainder = n_ab - sum(draw_this_many)
  cat(sprintf('Drawing %i Abs per individual. %i were requested. Remainder is %i.\n', sum(draw_this_many), n_ab, remainder))
  
  library(foreach)
  foreach(ee = native_epitope_coords$epitope, 
          aa = native_epitope_coords$antigen, 
          c1.mean = native_epitope_coords$c1, 
          c2.mean = native_epitope_coords$c2, 
          nn = draw_this_many) %do% {
    tibble(epitope = ee,
           antigen = aa,
           kind = 'antibody',
           c1 = rnorm(nn, c1.mean, sd = sigma),
           c2 = rnorm(nn, c2.mean, sd = sigma))
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
    filter(antigen == serum_id) %>%
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
    #cat(print('\n'))
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


format_stan_inputs <- function(titer_map){
  antigens = unique(titer_map$antigen)
  sera = unique(titer_map$serum)
  
  distmat = matrix(NA, length(antigens), length(sera))
  for(aa in 1:length(antigens)){
    for(ss in 1:length(sera)){
      distmat[aa,ss] = filter(titer_map, antigen == aa & serum == ss) %>%
        pull(titer_distance)
    }
  }
  distmat
}


fit_stan_MDS <- function(
  mod = 'MDS.stan',
  observed_distances, # n_antigen x n_antibody matrix of distances
  n_antigens, # integer
  n_sera, # integer
  n_dim, # integer
  chains = 3, # Number of MCMC chains to run
  cores = parallel::detectCores(logical = F), # For the cluster
  niter = 5000,
  ...
) {
  library(rstan)
  
  stopifnot(nrow(observed_distances) == n_antigens)
  stopifnot(ncol(observed_distances) == n_sera)
  
  model <- stan_model('../Bayesian_stan/MDS.stan')
  (print(model))
  
  model_input_data <- list(
    n_strains = n_antigens,
    n_sera = n_sera,
    n_dim = n_dim,
    observed_distances = observed_distances
  )
  
  initfun <- function(){
    list(sigma = 1,
         ag2_c1 = runif(1, 0, 10),
         strain_coords = matrix(runif((n_antigens-2)*n_dim, -10, 10), n_antigens-2, n_dim),
         serum_coords =  matrix(runif(n_sera*n_dim, -10, 10), n_sera, n_dim)
    )
  }
  
  fit <- sampling(
    model, model_input_data, chains = chains, cores = cores, 
    init = initfun,
    iter = niter,
    control = list(adapt_delta = 0.89,
                   max_treedepth = 14),
    ...
  )
}




## Wrapper function to do it all
infer_titer_map <- function(antigen_coords, # Data frame of ag coords 
                            relative_concentrations, # vector of relative concentrations (immunodominance) of Abs to each epitope. WILL NOT BE NORMALIZED. SHOULD SUM TO N EPITOPES.
                            n_epitopes, # n epitiopes
                            n_antigens,
                            n_abs_per_serum = 500,
                            plotdir,
                            sigma = 1 # sd of abs around native coord
                            ){ # n antigens
  
  if(sum(relative_concentrations) != n_epitopes){warning('Relative concentrations dont sum to n_epitopes. This will rescale the total Ab concentration.')}
  stopifnot(length(unique(antigen_coords$antigen)) == n_antigens)
  stopifnot(length(unique(antigen_coords$epitope)) == n_epitopes)
  
  ## Generate ferret repertoire
  ferret_repertoires <- lapply(1:n_antigens, function(this_ag){
    generate_ferret_repertoire(native_epitope_coords = antigen_coords %>% filter(antigen == this_ag), 
                               n_epitopes = n_epitopes,
                               n_ab = n_abs_per_serum, ## n antibodies per ferret
                               rel_immuno = relative_concentrations,
                               sigma = sigma)
  }) %>%
    bind_rows(.id = 'antigen')
  
  ## Merge the antigen coordinates with the Ab coordinates to calculate titers
  merged_df <- merge(antigen_coords %>%
                       select(epitope, antigen, c1, c2),
                     ferret_repertoires %>%
                       select(epitope, antigen, c1, c2),
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
                                               relative_concentrations = relative_concentrations, 
                                               ag_list = ag_list, 
                                               alpha = .25,
                                               r = 7)
    )
    c('serum' = this.serum, 'antigen' = this.antigen, 'titer' = titer)
  }
  
  titer_map <- with(expand.grid(antigen = 1:n_ag, serum = 1:n_ag),
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
  
  ## Fit the stan model
  fits <- fit_stan_MDS(observed_distances = format_stan_inputs(titer_map), 
                       n_antigens = n_ag, 
                       n_sera = n_ag, 
                       n_dim = n_dim, 
                       chains = n_chains,
                       niter = n_iter)
  
  
  
  return(list(ag_ab_coords = merged_df,
              titer_map = titer_map,
              summary_inferred_coords = summary_coords,
              stan_fits = fits))
}




plot_fits <- function(antigen_coords,
                      fits, 
                      n_iter,
                      outdir){
  trplot <- rstan::traceplot(fits, pars = names(fits))
  
  contour_posteriors <- extract_long_coords(fits) %>% 
    filter(iter %% 5 == 0) %>% ## THIN
    ggplot()+
    geom_density_2d(aes(x = c1, y = c2, color = factor(chain))) +
    facet_grid(kind ~ id)
  
  epitope_strain_map <- plot_inferred_original_map(antigen_coords, 
                                                   fits)
  
  outdir = paste0(Sys.Date(), '/', outdir)
  if(!dir.exists(paste0(Sys.Date()))) dir.create(paste0(Sys.Date()))
  if(!dir.exists(outdir))  dir.create(outdir)
  ggsave(paste0(outdir, '/traceplot.png'), plot = trplot)
  ggsave(paste0(outdir, '/epitope_strain_map.png'), plot = epitope_strain_map)
  ggsave(paste0(outdir, '/contour_posteriors.png'), plot = contour_posteriors)
  
  ## pairs plots
  png(paste0(outdir, '/antigen_pairs.png'))
  pairs(fits, pars = names(fits)[grepl(names(fits), pattern = 'antigen.+')])
  dev.off()
  png(paste0(outdir, '/serum_pairs.png'))
  pairs(fits, pars = names(fits)[grepl(names(fits), pattern = 'serum.+')])
  dev.off()
  png(paste0(outdir, '/other_pairs.png'))
  pairs(fits, pars = names(fits)[!grepl(names(fits), pattern = 'antigen.+') & !grepl(names(fits), pattern = 'serum.+')])
  dev.off()
  
  
  cat(sprintf('plots saved in %s', outdir))
  return(epitope_strain_map)
}


plot_antigen_coords <- function(antigen_coords){
  antigen_coords %>%
    ggplot() +
    geom_line(aes(x=c1, y = c2, group = epitope), lty = 2, lwd = .1)+
    geom_point(aes(x = c1, y = c2, shape = antigen, color = epitope))+
    facet_wrap(.~epitope, labeller = label_both)
}


plot_Ab_Ag_map <- function(ab_ag_df){
  ab_ag_df %>%
    ggplot() +
    geom_point(aes(x = c1_Ab, y = c2_Ab, color = epitope), pch = 3) +
    geom_point(aes(x = c1_Ag, y = c2_Ag, fill = epitope), pch = 21) +
  #  geom_point(aes(x = c1, y = c2, fill = epitope), data = antigen_coords, pch = 22) +
    guides(color=guide_legend(title="Abs to epitope"))+
    xlab('c1')+
    ylab('c2')
}

plot_inferred_original_map <- function(antigen_coords,
                                       stan_fit){

  
  reformatted_df <- antigen_coords %>%
    mutate(kind = paste0('E', epitope)) %>%
    rename(allele = antigen) %>%
    select(allele, kind, starts_with('c')) %>%
    mutate(kind = factor(kind, levels = c('antigen', 'serum', 'fixed_ag1', 'E1', 'E2', 'E3'),
                         labels = c('antigen (estimated)', 'serum (estimated)', 'antigen 1 position (fixed)', 'epitope 1 (truth)', 'epitope 2 (truth)', 'epitope 3 (truth)')), ordered = T)
  
  epitope_strain_map <- extract_summary_coords(stan_fit) %>%
    mutate(allele = as.factor(id)) %>%
    mutate(kind = factor(kind, levels = c('antigen', 'serum', 'fixed_ag1', 'E1', 'E2', 'E3'),
                         labels = c('antigen (estimated)', 'serum (estimated)', 'antigen 1 position (fixed)', 'epitope 1 (truth)', 'epitope 2 (truth)', 'epitope 3 (truth)')), ordered = T) %>%
    ggplot() +
    geom_point(aes(x = c1, y = c2, shape = kind), data = reformatted_df) +
    geom_polygon(aes(x = c1, y = c2, fill = allele), data = reformatted_df, alpha = .2)+
    geom_point(aes(x = c1, y = c2, shape = kind, color = allele)) +
    facet_grid(allele~chain, labeller = label_both) +
    scale_shape_manual(values = c(4, 14, 21, 22, 23, 3))
  
  epitope_strain_map
}



extract_long_coords <- function(stan_fit){
  raw_fits <- rstan::extract(stan_fit)
  
  long_antigen_coords <- lapply(1:(n_ag-2), function(ll) raw_fits$antigen_coords[,ll,] %>%
                                  magrittr::set_colnames(c('c1', 'c2')) %>%
                                  as_tibble() %>%
                                  mutate(chain = rep(1:n_chains, each = n_iter/2),
                                         id = ll+2,
                                         iter = 1:nrow(.),
                                         kind = 'antigen')) %>%
    bind_rows()
  
  long_serum_coords <- lapply(1:n_ag, function(ll) raw_fits$serum_coords[,ll,] %>%
                                magrittr::set_colnames(c('c1', 'c2')) %>%
                                as_tibble() %>%
                                mutate(chain = rep(1:n_chains, each = n_iter/2),
                                       id = ll,
                                       iter = 1:nrow(.),
                                       kind = 'serum')) %>%
    bind_rows()
  
  ag2_long <- tibble(id = 2, kind = 'antigen', chain = rep(1:n_chains, each = n_iter/2), c1 = raw_fits$ag2_c1) 
  
  long_coords <- bind_rows(long_antigen_coords,
                           long_serum_coords,
                           ag2_long)
  long_coords
}

  
extract_summary_coords <- function(stan_fit){
  summary_coords <- extract_long_coords(stan_fit) %>%
    group_by(id, kind, chain) %>%
    summarise(c1.10 = quantile(c1, .1),
              c1 = median(c1),
              c1.90 = quantile(c1, .9),
              c2.10 = quantile(c2, .1, na.rm = T),
              c2 = median(c2, na.rm = T),
              c2.90 = quantile(c2, .9, na.rm = T)) %>%
    ungroup() %>%
    bind_rows(tibble(id = 1,
                     kind = 'fixed_ag1', chain = 1:3, c1.10 = NA, c1 = 0, c1.90 = NA, c2.10 = NA, c2 = 0, c2.90 = NA)) 
  summary_coords
}


flip_summary_over_x <- function(original_coords,
                        inferred_summary){
 flip.these.chains <-  merge(original_coords, 
        inferred_summary %>%
          filter(kind == 'antigen') %>%
          rename(antigen = id),
        by = 'antigen') %>%
    group_by(chain) %>%
    summarise(flip.me = sum((c2.x > 0) == (c2.y > 0)) < sum((c2.x > 0) == (-c2.y > 0)))
 
 inferred_summary %>% 
   merge(flip.these.chains) %>%
   mutate(across(starts_with("c2"), negate))
  
}

negate <- function(x){
  -x
}
