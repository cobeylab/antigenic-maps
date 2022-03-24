
reset <- function(){
  rm(list = ls())
  library(tidyr)
  filter <- dplyr::filter
  extract <- tidyr::extract
  select <- dplyr::select
  source('../code.R')
}


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
                                         n_ab = 1000, ## Total number of Abs to draw
                                         rel_immuno, ## Vector giving the relative immunodominance of each epitope. Will be normalized. The fraction of Abs specific to epitope i will be proportional to the ith entry.
                                         sigma = .1 ## sd of Ab positions around native_epitope_Coords
){
  ## Figure out how many Abs to draw for each epitope
  draw_this_many = floor(n_ab*rel_immuno/sum(rel_immuno))
  remainder = n_ab - sum(draw_this_many)
  cat(sprintf('Drawing %i Abs per individual. %i were requested. Remainder is %i.\n', sum(draw_this_many), n_ab, remainder))
  stopifnot(length(draw_this_many) == nrow(native_epitope_coords))
  
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

format_stan_titers <- function(titer_map){
  antigens = unique(titer_map$antigen)
  sera = unique(titer_map$serum)
  
  titer_matrix = matrix(NA, length(antigens), length(sera))
  for(aa in 1:length(antigens)){
    for(ss in 1:length(sera)){
      titer_matrix[aa,ss] = filter(titer_map, antigen == aa & serum == ss) %>%
        pull(logtiter)
    }
  }
  titer_matrix
}

format_stan_smax <- function(titer_map){
  antigens = unique(titer_map$antigen)
  sera = unique(titer_map$serum)
  
  outmat = matrix(NA, length(antigens), length(sera))
  for(aa in 1:length(antigens)){
    for(ss in 1:length(sera)){
      outmat[aa,ss] = filter(titer_map, antigen == aa & serum == ss) %>%
        mutate(smax = mean(c(serum_potency, antigen_avidity))) %>%
        pull(smax)
    }
  }
  outmat
}


fit_stan_MDS <- function(
  mod,
  observed_titers, # N-vector of observed titers
  smax, # N-vector of smax
  serum_id, # N-vector of serum id
  antigen_id, # N- vector of antigen id
  n_antigens, # integer
  n_sera, # integer
  n_dim, # integer
  N_test_set,
  smax_test_set,
  serum_id_test_set,
  antigen_id_test_set,
  observed_titers_test_set,
  chains = 3, # Number of MCMC chains to run
  cores = parallel::detectCores(logical = F), # For the cluster
  niter = 5000,
  antigen_coords,
  diagnostic_file = NULL,
  ...
) {
  
  
  model <- stan_model(mod)
  #(print(model))
  
  model_input_data <- list(
    N = length(observed_titers),
    n_sera = n_sera,
    n_antigens = n_antigens,
    n_dim = n_dim,
    observed_titers = observed_titers,
    smax = smax,
    serum_id = serum_id,
    antigen_id = antigen_id,
    N_test_set = N_test_set,
    smax_test_set = smax_test_set,
    serum_id_test_set = serum_id_test_set,
    antigen_id_test_set = antigen_id_test_set,
    observed_titers_test_set = observed_titers_test_set
  )
  
  initfun <- function(){
    list(sigma = 1,
         ag2_c1 = runif(1, 0, 10),
         strain_coords = matrix(runif((n_antigens-2)*n_dim, -10, 10), n_antigens-2, n_dim),
         serum_coords =  matrix(runif(n_sera*n_dim, -10, 10), n_sera, n_dim)
    )
  }
  
  initial_fit <- sampling(model, 
                          data = model_input_data, 
                          chains = chains, 
                          init = initfun,
                          iter = niter,
                          cores = min(6, chains),
                          control = list(adapt_delta = 0.89,
                                         max_treedepth = 14),
                          diagnostic_file = diagnostic_file
  )
  
  if(all(summary(initial_fit)$summary[,'Rhat'] <= 1.03)){
    cat(sprintf('Returning initial fit'))
    return(initial_fit)
  }
  
  cat(print('Initial fit complete.\nRe-running using initial values from the best chain.\n'))
  
  initialize_with_best_chain <- function(initial_fit, nchains){
    ## Extract the summary of the best chain in terms of mean log posterior
    best_chain_index = which.max(rstan::extract(initial_fit, permute = F)[,,'lp__'] %>% colMeans())
    best_chain_summary = rstan::summary(initial_fit)$c_summary[,,best_chain_index]
    best_chain = rstan::extract(initial_fit, permute = F)[,best_chain_index,]
    
    is.antigen.coord = sapply(dimnames(best_chain_summary)$parameter, FUN = function(xx) grepl('antigen_coords', xx))
    is.serum.coord = sapply(dimnames(best_chain_summary)$parameter, FUN = function(xx) grepl('serum_coords', xx))
    
    
    get_one_list <- function(){
      ## Output an initial list using the median values from the best chain
      list(sigma = best_chain_summary['sigma', '50%'],
           ag2_c1 = best_chain_summary['ag2_c1', '50%'],
           antigen_coords = matrix(best_chain_summary[is.antigen.coord,'50%'], 
                                   nrow = n_antigens-2, 
                                   ncol = n_dim,
                                   byrow = T),
           serum_coords = matrix(best_chain_summary[is.serum.coord,'50%'], 
                                 nrow = n_sera, 
                                 ncol = n_dim,
                                 byrow = T))
    }
    
    lapply(1:nchains, function(xx){get_one_list()})
  }
  
  inits <- initialize_with_best_chain(initial_fit, 1)
  
  refit <- sampling(
    model, model_input_data, 
    chains = 1, 
    init = inits,
    iter = niter,
    diagnostic_file = diagnostic_file)
  
  if(! all(summary(refit)$summary[,'Rhat'] <= 1.03)){
    cat(sprintf('Re-doing refit to achieve Rhat < 1.03'))
    refit <- sampling(
      model, model_input_data, 
      chains = 1, 
      init = inits,
      iter = niter,
      diagnostic_file = diagnostic_file)
  }
  
  return(refit)
}





## Generate serum panel, and titer map, for a given set of antigen coords and an immunodominance scheme
generate_ferret_inputs <- function(antigen_coords, # Data frame of ag coords 
                            relative_immunodominance, # vector of relative concentrations (immunodominance) of Abs to each epitope.
                            n_epitopes, # n epitiopes
                            n_antigens,
                            n_abs_per_serum = 500,
                            sigma = 1 # sd of abs around native coord
                            ){ # n antigens
  
  if(sum(relative_immunodominance) != n_epitopes){warning('Relative concentrations dont sum to n_epitopes. This will rescale the total Ab concentration.')}
  stopifnot(length(unique(antigen_coords$antigen)) == n_antigens)
  stopifnot(length(unique(antigen_coords$epitope)) == n_epitopes)
  
  ## Generate ferret repertoire
  ferret_repertoires <- lapply(1:n_antigens, function(this_ag){
    generate_gaussian_repertoire(native_epitope_coords = antigen_coords %>% filter(antigen == this_ag), 
                               n_epitopes = n_epitopes,
                               n_ab = n_abs_per_serum, ## n antibodies per ferret
                               rel_immuno = relative_immunodominance,
                               sigma = sigma)
  }) %>%
    bind_rows(.id = 'serum')
  
  ## Merge the antigen coordinates with the Ab coordinates to calculate titers
  merged_df <- merge(antigen_coords %>%
                       select(epitope, antigen, c1, c2),
                     ferret_repertoires %>%
                       select(serum, epitope, antigen, c1, c2),
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
  
  return(list(ag_ab_coords = merged_df,
              titer_map = titer_map))
  
}
  

# 
#   ## Fit the stan model
#   fits <- fit_stan_MDS(observed_distances = format_stan_inputs(titer_map), 
#                        n_antigens = n_antigens, 
#                        n_sera = n_antigens, 
#                        n_dim = n_dim, 
#                        chains = n_chains,
#                        niter = n_iter)
#   
#   
# 
# }




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



extract_long_coords <- function(stan_fit){
  raw_fits <- rstan::extract(stan_fit, permuted = F, inc_warmup = F)
  niter <- dim(raw_fits)[1]
  nchains <- dim(raw_fits)[2]
  
  parnames = dimnames(raw_fits)[3]
  parname_contains <- function(parnames, this_pattern){sapply(parnames, function(xx) grepl(pattern = this_pattern, x = xx)) %>% as.vector()}
  is_antigen_coord = parname_contains(parnames, 'antigen_coords')
  is_serum_coord = parname_contains(parnames, 'serum_coords')
  is_log_posterior = parname_contains(parnames, 'lp_')
  
  
  # ag_coords = raw_fits[,,is_antigen_coord]
  # serum_coords = raw_fits[,,is_serum_coord]
  # ag2_coords = raw_fits[,,parname_contains(parnames, 'ag2')]
  
  #cat(print('checkpoint 1\n'))
  long_antigen_coords <- lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_antigen_coord]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    mutate(kind = 'antigen') %>%
    pivot_longer(contains('antigen_coords'), names_to = c('allele', 'coordinate'), names_pattern = 'antigen_coords.(\\d+),(\\d+).', values_to = 'value') %>%
    mutate(allele = as.numeric(allele) + 2) %>%
    pivot_wider(id_cols = c(chain, iter, allele, kind), names_from = coordinate, names_prefix = 'c', values_from = value) 
  
  #cat(print('checkpoint 2\n'))
  long_serum_coords <- lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_serum_coord]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    mutate(kind = 'serum') %>%
    pivot_longer(contains('serum_coords'), names_to = c('allele', 'coordinate'), names_pattern = 'serum_coords.(\\d+),(\\d+).', values_to = 'value') %>%
    mutate(allele = as.numeric(allele)) %>%
    pivot_wider(id_cols = c(chain, iter, allele, kind), names_from = coordinate, names_prefix = 'c', values_from = value) 
  
  #cat(print('checkpoint 3\n'))
  ag2_long <- lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,parname_contains(parnames, 'ag2')]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    rename(c1 = value) %>%
    mutate(allele = 2) %>%
    mutate(kind = 'antigen') %>%
    select(chain, iter, allele, kind, c1)
  
  ## If more than 3 dimensions, pad ag2 data frame with 0s
  if(any(! (colnames(long_antigen_coords) %in% colnames(ag2_long)) )) {
    add_these_columns <- colnames(long_antigen_coords)[!(colnames(long_antigen_coords) %in% colnames(ag2_long))]
    for(cc in add_these_columns){
      ag2_long[[cc]] = 0
    }  
  }
  
  n_allele <- max(long_serum_coords$allele)
  coord_names = long_serum_coords %>% select(matches('c\\d$')) %>% colnames()
  fill_vec = numeric(length(coord_names)+1)
  names(fill_vec) = c('iter', coord_names)
  fill_vec['iter'] = NA
  
  # (print(long_antigen_coords))
  # (print(long_serum_coords))
  # (print(ag2_long))
  
  long_coords <- bind_rows(long_antigen_coords,
                           long_serum_coords,
                           ag2_long)  %>%
    complete(allele, kind, chain, fill = as.list(fill_vec)) ## Add 0 coordinates for the fixed antigen at the origin
  long_coords
}

  
extract_summary_coords <- function(stan_fit){
  summary_coords <- extract_long_coords(stan_fit) %>%
    group_by(allele, kind, chain) %>%
    summarise_at(.vars = vars(matches('c\\d+')), .funs = list(`.500` = median,
                                                              `.025` = ~quantile(.x, probs = .025),
                                                              `.100` = ~quantile(.x, probs = .025),
                                                              `.900` = ~quantile(.x, probs = .025),
                                                              `.975` = ~quantile(.x, probs = .025))) %>%
    rename_with(.fn = ~gsub("_.500", "", .x, fixed = TRUE), .cols = contains('.500'))
  if(ncol(summary_coords) == 8){
    summary_coords <- summary_coords %>%
      rename_with(.fn = ~gsub(".", "c1.", .x, fixed = TRUE), .cols = matches('\\d\\d\\d')) %>%
      rename(c1 = c1.500)
  }
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
    summarise(flipped_error = sum((c2.x + c2.y)^2, na.rm = T),
              original_error = sum((c2.x - c2.y)^2, na.rm = T),
              flip.me = flipped_error<original_error)
 cat(sprintf('chain %s flipped over x\n', flip.these.chains %>% filter(flip.me) %>% pull(chain)))
 print(flip.these.chains)
 
 inferred_summary %>% 
   merge(flip.these.chains) %>%
   mutate(across(starts_with("c2"), negate))
  
}

negate <- function(x){
  -x
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




## Test this tomorrow

infer_human_map <- function(ab_ag_df, # a long data frame of antigen coords and corresponding antibodies with an id column for distinct sera
                            n_dim,
                            n_chains = 4,
                            n_iter = 5000){ 
  
  n_antigen = length(unique(ab_ag_df$antigen))
  n_sera = length(unique(ag_ag_df$serum))
  
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
  
  titer_map <- with(expand.grid(antigen = 1:n_antigen, serum = 1:n_sera),
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
                       n_antigens = n_antigen, 
                       n_sera = n_sera, 
                       n_dim = n_dim, 
                       chains = n_chains,
                       niter = n_iter)
  
  
  
  return(list(ag_ab_coords = merged_df,
              titer_map = titer_map,
              summary_inferred_coords = extract_summary_coords(fits),
              stan_fits = fits))
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


plot_inferred_original_map2 <- function(stan_fits,
                                        weights_for_shift, # An E-vector giving the weight of each epitope for calculating the centroids
                                        antigen_coords,
                                        flip_these_chains_over_x){
  inferred_map <- extract_summary_coords(stan_fits) %>%
    mutate(id = as.factor(id)) %>%
    mutate(c2 = ifelse(chain %in% flip_these_chains_over_x, -c2, c2))
  
  ag1_ag2_centroids = get_weighted_centroids(antigen_coords, weights_for_shift)
  
  true_coords = standardize_coordinate_df(coord_df = bind_rows(ag1_ag2_centroids, antigen_coords), # Standardize to Ag1, Ag2 centroids
                                          ag1_row = 1, ag2_row = 2) %>%
    filter(!is.na(epitope)) %>% # Remove ag1, ag2 centroids
    group_by(antigen) %>%
    mutate(centroid1 = sum(c1*weights_for_shift),
           centroid2 = sum(c2*weights_for_shift))
  
  ggplot() +
    geom_polygon(aes(x = c1, y = c2, fill = antigen), data = true_coords, alpha = .05) +
    geom_point(aes(x = centroid1, y = centroid2, color = antigen), data = true_coords, pch = 1) +
    geom_point(aes(x = c1, y = c2, color = id, shape = kind), data = inferred_map)
}



compare_distances <- function(fitlist){
  nchains = (rstan::extract(fitlist$stan_fits, permuted = F) %>% dim())[2]
  fitlist$titer_map %>%
    expand_grid(chain = 1:nchains) %>%
    mutate(chain = as.numeric(chain)) %>%
    rowwise() %>%
    mutate(inferred_distance = get_one_inferred_distance(serum, antigen, chain, fitlist$summary_inferred_coords)) %>%
    filter(!is.na(inferred_distance)) %>%
    ungroup() %>%
    arrange(inferred_distance)
}

plot_compare_distances <- function(fitlist){
  compare_distances <- compare_distances(fitlist)
  sum_square_error = compare_distances %>%
    summarise(rs_error = sum((inferred_distance-titer_distance)^2))
  
  ggplot(compare_distances) +
    geom_point(aes(x = (titer_distance), y = (inferred_distance), color = as.factor(chain)), alpha = .3, size = 3) +
    geom_abline(lty = 2, lwd = .5) +
    theme(legend.position = c(.25, .75)) +
    ggtitle(sprintf('sum square error is %2.2f', sum_square_error))
}

get_one_inferred_distance <- function(serum, 
                                      antigen, 
                                      this.chain,
                                      summary_coords){
  this.serum = summary_coords %>% 
    filter(id == serum & kind == 'serum') %>%
    filter(chain == this.chain)
  this.antigen = summary_coords %>%
    filter(id == antigen & kind == 'antigen') %>%
    filter(chain == this.chain)
  # each_chain = foreach(ag.c1 = this.antigen$c1, 
  #                      ag.c2 = this.antigen$c2,
  #                      ab.c1 = this.serum$c1,
  #                      ab.c2 = this.serum$c2,
  #                      .combine = 'c') %do% {
  #                        get_euclidean_distance(v1 = c(ab.c1, ab.c2),
  #                                               v2 = c(ag.c1, ag.c2))
  #                      }
  # mean(each_chain)
  stopifnot(nrow(this.serum) <= 1)
  stopifnot(nrow(this.antigen) <= 1)
  if(any(nrow(this.serum)==0, nrow(this.antigen) == 0)){return(NA)}
  get_euclidean_distance(v1 = c(this.serum$c1, this.serum$c2),
                         v2 = c(this.antigen$c1, this.antigen$c2))
}

plot_compare_distance_tiles <- function(fitlist){
  compare_distances(fitlist) %>%
    mutate(distance_difference = titer_distance - inferred_distance) %>% 
    pivot_longer(contains('distance')) %>% 
    ggplot() + 
    geom_tile(aes(x = antigen, y = serum, fill = value))   +
    facet_wrap(.~name) +
    scale_fill_viridis_c()
}


## Write a function to compare distances that would have occurred with centroid solution
centroid_counterfactual <- function(fitlist, 
                                          centroid_weights,
                                           antigen_coords){
  
  centroid_df <- get_weighted_centroids(antigen_coords, centroid_weights)
  get_one_centroid_distance <- function(ag1, ag2){
    v1 = c(centroid_df$c1[ag1], centroid_df$c2[ag1])
    v2 = c(centroid_df$c1[ag2], centroid_df$c2[ag2])
    get_euclidean_distance(v1, v2)
  }
  
  distance_comparison <- fitlist$titer_map %>%
    rowwise() %>%
    mutate(centroid_distance = get_one_centroid_distance(ag1 = antigen, ag2 = serum)) %>%
    ungroup() %>%
    arrange(centroid_distance) %>%
    mutate(distance_difference = titer_distance - centroid_distance) 
  
  tiles = distance_comparison %>% 
    pivot_longer(contains('distance')) %>%
    ggplot() + 
    geom_tile(aes(x = antigen, y = serum, fill = value))   +
    facet_wrap(.~name) +
    scale_fill_viridis_c()
  
  sum_square_error = distance_comparison %>%
    summarise(rs_error = sum((centroid_distance-titer_distance)^2))
  
  lineplot = ggplot(distance_comparison) +
    geom_point(aes(x = (titer_distance), y = (centroid_distance), color = as.factor(antigen)), alpha = .3, size = 3) +
    geom_abline(lty = 2, lwd = .5) +
    theme(legend.position = c(.25, .75)) +
    ggtitle(sprintf('sum square error is %2.2f', sum_square_error))
  
  cowplot::plot_grid(lineplot, tiles, nrow = 2)
  
}



get_titer_error <- function(stan_fit,
                             titer_map){
  raw_fits <- rstan::extract(stan_fit, permuted = F, inc_warmup = F)
  niter <- dim(raw_fits)[1]
  nchains <- dim(raw_fits)[2]
  
  parnames = dimnames(raw_fits)[3]
  parname_contains <- function(parnames, this_pattern){sapply(parnames, function(xx) grepl(pattern = this_pattern, x = xx)) %>% as.vector()}
  is_map_distance = parname_contains(parnames, 'map')
  is_titer = parname_contains(parnames, 'estimated_titers')
  
  stopifnot(all(titer_map %>% group_by(antigen, serum) %>% summarise(n = n()) %>% pull(n) == 1 ))
  titer_map = titer_map %>% 
    mutate(id = 1:nrow(.)) 
  
full_posterior_error <- merge(
  lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_map_distance]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    pivot_longer(contains('map'), names_to = c('id'), names_pattern = 'map_distances.(\\d+).', values_to = 'map_distance'),
  lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_titer]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    pivot_longer(contains('titer'), names_to = c('id'), names_pattern = 'estimated_titers.(\\d+).', values_to = 'estimated_titer')
  )%>%
  merge(titer_map, 
        by = c('id')) %>%
  arrange(serum, antigen, chain, iter)

  pairwise_errors = full_posterior_error %>%
  ungroup() %>% 
  mutate(titer_error = (sqrt((estimated_titer-logtiter)^2)),
         distance_error = (sqrt((map_distance-titer_distance)^2))) %>%
    group_by(antigen, serum, id) %>%
  dplyr::summarise(mean_map_distance = mean(map_distance),
                   titer_distance = unique(titer_distance),
                   mean_estimated_titer = mean(estimated_titer),
                   observed_titer = unique(logtiter),
                   titer_error_med = median(titer_error),
                   titer_error_0.025 = quantile(titer_error, 0.025),
                   titer_error_0.975 = quantile(titer_error, 0.975),
                   distance_error_med = median(distance_error),
                   distance_error_0.025 = quantile(distance_error, 0.025),
                   distance_error_0.975 = quantile(distance_error, 0.975)) %>%
    ungroup()
  
  model_errors = full_posterior_error %>%
    mutate(titer_error = (sqrt((estimated_titer-logtiter)^2)),
           distance_error = (sqrt((map_distance-titer_distance)^2))) %>%
    group_by(chain, iter) %>%
    summarise(distance_error = sum(distance_error),
              titer_error = sum(titer_error)) %>%
    ungroup() %>%
    dplyr::summarise(titer_error_med = median(titer_error),
                     titer_error_0.025 = quantile(titer_error, 0.025),
                     titer_error_0.975 = quantile(titer_error, 0.975),
                     distance_error_med = median(distance_error),
                     distance_error_0.025 = quantile(distance_error, 0.025),
                     distance_error_0.975 = quantile(distance_error, 0.975)) %>%
    ungroup()
  
  return(list(model_errors = model_errors,
              pairwise_errors = pairwise_errors))
}





get_predictive_errors <- function(stan_fit, 
                                 test_set){
  raw_fits <- rstan::extract(stan_fit, permuted = F, inc_warmup = F)
  niter <- dim(raw_fits)[1]
  nchains <- dim(raw_fits)[2]
  
  parnames = dimnames(raw_fits)[3]
  parname_contains <- function(parnames, this_pattern){sapply(parnames, function(xx) grepl(pattern = this_pattern, x = xx)) %>% as.vector()}
  is_predicted_titer = parname_contains(parnames, 'predicted_titers')
  is_predictive_error = parname_contains(parnames, 'posterior_predictive_titer_error')
  
  stopifnot(all(test_set %>% group_by(antigen, serum) %>% summarise(n = n()) %>% pull(n) == 1 ))
  test_set = test_set %>% 
    mutate(id = 1:nrow(.)) 
  
full_posterior <- merge(
  ## Reformatted predicted titers
  lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_predicted_titer]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    pivot_longer(contains('predicted'), names_to = c('id'), names_pattern = 'predicted_titers.(\\d+).', values_to = 'predicted_log_titer'),
  ## Reformatted predictive error
  lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_predictive_error]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    pivot_longer(contains('predictive'), names_to = c('id'), names_pattern = 'posterior_predictive_titer_error.(\\d+).', values_to = 'predictive_error')
  ) %>%
    merge(test_set) %>%
    arrange(id, chain, iter) %>%
  ungroup()


pairwise_errors <- full_posterior %>%
  group_by(antigen, serum, id)  %>%
      dplyr::summarise(titer_error_med = median(predictive_error),
                       titer_error_0.025 = quantile(predictive_error, 0.025),
                       titer_error_0.975 = quantile(predictive_error, 0.975)) %>%
      ungroup()

model_errors = full_posterior %>%
  group_by(chain, iter) %>% 
  summarise(titer_error = sum(predictive_error)) %>% ## Sum across data points
  ungroup() %>%
  dplyr::summarise(titer_error_med = median(titer_error),
                   titer_error_0.025 = quantile(titer_error, 0.025),
                   titer_error_0.975 = quantile(titer_error, 0.975)) %>%
  ungroup()

return(list(model_errors = model_errors,
            pairwise_errors = pairwise_errors))
}





## Calculate the error between titers predicted in a model with one immunodominance scheme, against titers observed under any other immunodominance scheme
get_generalized_titer_errors <- function(stan_fit, 
                                          test_set,
                                          observed_titers,
                                          fitted_immunodominance_scheme = 'E1',
                                          observed_immunodominance_scheme = 'even'){

  
  
  raw_fits <- rstan::extract(stan_fit, permuted = F, inc_warmup = F)
  niter <- dim(raw_fits)[1]
  nchains <- dim(raw_fits)[2]
  
  parnames = dimnames(raw_fits)[3]
  parname_contains <- function(parnames, this_pattern){sapply(parnames, function(xx) grepl(pattern = this_pattern, x = xx)) %>% as.vector()}
  is_predicted_titer = parname_contains(parnames, 'predicted_titers')
  is_predictive_error = parname_contains(parnames, 'posterior_predictive_titer_error')
  
  stopifnot(all(test_set %>% group_by(antigen, serum) %>% summarise(n = n()) %>% pull(n) == 1 ))
  test_set = test_set %>% 
    mutate(id = 1:nrow(.)) 
  
  ## Take the observed titer map and extract the entries in the test set
  observed_subset = merge(test_set %>% 
                            select(id, serum, antigen),
                          observed_titers,
                          all.x = T, all.y = F)
  
  full_posterior <- merge(
    ## Reformatted predicted titers
    lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_predicted_titer]) %>% mutate(iter = 1:niter)}) %>%
      bind_rows(.id = 'chain') %>%
      pivot_longer(contains('predicted'), names_to = c('id'), names_pattern = 'predicted_titers.(\\d+).', values_to = 'predicted_log_titer'),
    ## Reformatted predictive error
    lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_predictive_error]) %>% mutate(iter = 1:niter)}) %>%
      bind_rows(.id = 'chain') %>%
      pivot_longer(contains('predictive'), names_to = c('id'), names_pattern = 'posterior_predictive_titer_error.(\\d+).', values_to = 'predictive_error')
  ) %>%
    merge(observed_subset) %>%
    arrange(id, chain, iter) %>%
    ungroup()
  
  
  ## calculate pairwise errors and return
  pairwise_errors <- full_posterior %>%
    mutate(titer_error = sqrt((predicted_log_titer - logtiter)^2)) %>%
    group_by(antigen, serum, id)  %>%
    dplyr::summarise(titer_error_med = median(titer_error),
                     titer_error_0.025 = quantile(titer_error, 0.025),
                     titer_error_0.975 = quantile(titer_error, 0.975)) %>%
    ungroup() %>%
    mutate(kind = sprintf('%s_predictions:%s_observations', fitted_immunodominance_scheme, observed_immunodominance_scheme))
  
  overall_errors = full_posterior %>%
    mutate(titer_error = sqrt((predicted_log_titer - logtiter)^2)) %>%
    group_by(chain, iter) %>% 
    summarise(titer_error = sum(titer_error)) %>% ## Sum across data points
    ungroup() %>%
    dplyr::summarise(titer_error_med = median(titer_error),
                     titer_error_0.025 = quantile(titer_error, 0.025),
                     titer_error_0.975 = quantile(titer_error, 0.975)) %>%
    ungroup() %>%
    mutate(kind = sprintf('%s_predictions:%s_observations', fitted_immunodominance_scheme, observed_immunodominance_scheme))
  

  
  return(list(overall_errors = overall_errors,
              pairwise_errors = pairwise_errors))
}




run_pca_on_titer_map <- function(titer_matrix){
  
  this_pca <- prcomp(x = titer_matrix, scale = TRUE, center = TRUE)
  pca_summary <- summary(this_pca)
  this_ndim = ncol(pca_summary$importance)
  
  biplot <- ggbiplot::ggbiplot(this_pca, labels = rownames(titer_matrix)) +
    ylim(c(-2, 2)) + xlim(c(-2,2))
  
  var_explained_plot <- tibble(component = factor(1:this_ndim, labels = paste0('PC', 1:this_ndim)),
                               variance_explained = summary(this_pca)$importance[2,],
                               cumulative_var_explained = cumsum(variance_explained)) %>%
    ggplot() +
    geom_bar(aes(x = component, y = variance_explained), stat = 'identity') +
    geom_point(aes(x = component, y = cumulative_var_explained), color = 'deepskyblue') +
    geom_line(aes(x = as.numeric(component), y = cumulative_var_explained), color = 'deepskyblue') 
  
  return(list(pca_summary = pca_summary,
              biplot = biplot,
              var_explained_plot = var_explained_plot))
}


convert_titer_map_to_matrix <- function(titer_map){
  sera = unique(titer_map$serum)
  antigens = unique(titer_map$antigen)
  stopifnot(all(titer_map %>% group_by(serum, antigen) %>% dplyr::summarise(n = n()) %>% pull(n) == 1))
  titer_map <- titer_map %>% arrange(antigen, serum)
  matrix(titer_map$logtiter, nrow = length(sera), ncol = length(antigens), byrow = F, 
         dimnames = list(paste0('serum', sera), paste0('antigen', antigens)))
}


test_train_split <- function(N, n_test, df){
  test_indices = sample(1:N, size = n_test, replace = F)
  list(train = df[-test_indices,],
       test = df[test_indices,])
}
  