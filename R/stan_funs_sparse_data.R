## ---------------------------------------------------------------------
## ----------- Functions to extract columns from 
## ----------- the long-form "titer_map" inputs 
## ----------- and format into matrices
distance_matrix_format <- function(titer_map){
  max_allele = max(c(titer_map$antigen, titer_map$serum))
  ## Extract titer distances from the titer map and format as a matrix
  titer_map %>% select(serum, antigen, titer_distance) %>%
    complete(serum = 1:max_allele, antigen = 1:max_allele) %>%
    arrange(serum, antigen) %>%
    pivot_wider(names_from = antigen, values_from = titer_distance) %>%
    select(-serum) %>%
    as.matrix()
}



titer_matrix_format <- function(titer_map){
  max_allele = max(c(titer_map$antigen, titer_map$serum))
  titer_map %>% select(serum, antigen, logtiter)  %>%
    complete(serum = 1:max_allele, antigen = 1:max_allele) %>%
    arrange(serum, antigen) %>%
    pivot_wider(names_from = antigen, values_from = logtiter) %>%
    select(-serum) %>%
    as.matrix()
}


smax_matrix_format <- function(titer_map){
  max_allele = max(c(titer_map$antigen, titer_map$serum))
  ## Extract smax estimates from the titer map and format as a matrix
  titer_map %>% 
    mutate(smax = (serum_potency+antigen_avidity)/2) %>%
    select(serum, antigen, smax)  %>%
    complete(serum = 1:max_allele, antigen = 1:max_allele) %>%
    arrange(serum, antigen) %>%
    pivot_wider(names_from = antigen, values_from = smax) %>%
    select(-serum) %>%
    as.matrix()
}
## ---------------------------------------------------------------------




## ---------------------------------------------------------------------
##  Split into a test and training set
## ---------------------------------------------------------------------
## First, write an internal function to make sure the training set
##  contains all the necessary observations for the model to run
check_validity <- function(df, test_inds){
  train = df[-test_inds, ]
  ## 1. All antigens and sera must appear at least once in the training set
  check1.1 = c(check1.antigens = all(df$antigen %in% train$antigen))
  check1.2 = c(check1.sera = all(df$serum %in% train$serum))
  ## 2. In order to calculate constraints, at least one of the entries:
  ##  [i,j] or [j,i] must be in the training set for i, j in 1:8
  testmat = distance_matrix_format(train)
  check2.vec = vector('logical', 8)
  names(check2.vec) = paste0('check2_dim', 1:8)
  for(ii in 1:8){
    vec1 = testmat[ii,1:ii]; vec1
    vec2 = testmat[ii,1:ii]; vec2
    is.na(vec1)+is.na(vec2)
    
    check2.vec[ii] = all(is.na(vec1)+is.na(vec2) < 2)
  }
  c(check1.1, check1.2, check2.vec)
}


test_train_split <- function(df){
  ## Split a dataframe into 20% test set and 80% train
  N = nrow(df)
  n_test = floor(N*.2)
  ## Sample test indices from the upper triangle of the serum-antigen matrix
  ##  This avoids instances where an antigen-serum pair are not observed at all
  upper_triangle_indices = df %>% 
    arrange(serum, antigen) %>%
    mutate(id = 1:nrow(.)) %>%
    filter(antigen > serum) %>% 
    pull(id)
  ## Downsample
  cat(sprintf('Split the %s-row data frame into %s training rows (80%%) and %s test rows (20%%)\n', N, N-n_test, n_test))
  test_inds = sample(upper_triangle_indices, size = n_test, replace = F)
  ## Check validity of sample
  stopifnot(all(check_validity(df, test_inds) == TRUE))
  return(list(test = df[test_inds,],
              train = df[-test_inds,]))
}
## ---------------------------------------------------------------------




## ---------------------------------------------------------------------
## ------------  Function to calculate constraints on coordinates
## ------------  inferred in the model
get_constrainted_coords <-   function(distmat, 
                                      n_dim,
                                      verbose = F) {
  ## Constrain the first few coordinates
  constrained_coords = matrix(0, nrow = n_dim+1, ncol = n_dim)
  constrained_coords[2,1] = mean(c(distmat[2,1], distmat[1,2]), na.rm = T)
  
  if(n_dim >= 2){
    this_fun <- function(pars, 
                         this_dim){
      distances = colMeans(rbind(distmat[this_dim, 1:(this_dim-1)], 
                                 distmat[1:(this_dim-1), this_dim]), 
                           na.rm = T)
      target = 0
      for(ii in 1:(this_dim-1)){
        target = target + abs( sqrt(sum( (pars - constrained_coords[ii,1:(this_dim-1)])^2 )) - distances[ii] )
        #cat(sprintf('target = %2.2f; this dist is %2.2f; this estimate is %2.2f\n', target, distances[ii], sqrt(sum( (pars - constrained_coords[ii,1:(this_dim-1)])^2 ) )))
      }
      target
    }
    
    for(this_dim in 3:(n_dim+1)){
      pars = vector(length = this_dim-1) + 1
      names(pars) = letters[1:this_dim-1]
      solution = optim(pars, fn = this_fun, this_dim = this_dim, method = 'BFGS', 
                       control = list(abstol = 10^-7))
      constrained_coords[this_dim, 1:(this_dim-1)] = solution$par
    }
  }
  
  estimated_distances = get_ab_ag_distances(constrained_coords, constrained_coords)   ## Check that estimated distances from constrained coords match inputs
  if(!max(abs(estimated_distances - distmat[1:(n_dim+1), 1:(n_dim+1)]), na.rm = T) < 0.1){
    warning(sprintf('max error when inferring constraints is > 0.1. Try rerunning as verbose.'))
  }
  
  if(verbose == TRUE){
    print('Estimated distances are:\n')
    print(estimated_distances)
    print('Target distances are\n:')
    print( distmat[1:(n_dim+1),1:(n_dim+1)] )
    print('\n')
  }
  
  constrained_coords
}
## ---------------------------------------------------------------------


## ---------------------------------------------------------------------
## -------------Function to generate initial values, consistent with
## ------------ constraints
initfun = function(distmat,
                   n_dim,
                   n_antigens,
                   n_sera){
  ## Generate random guesses
  initlist <- list(sigma = 1,
                   antigen_coords = matrix(runif(n_antigens*n_dim, -10, 10), n_antigens, n_dim),
                   serum_coords = matrix(runif(n_antigens*n_dim, -10, 10), n_sera, n_dim))
  constrained_coords = get_constrainted_coords(distmat, n_dim)
  initlist$antigen_coords[1:(n_dim+1),] = constrained_coords
  return(initlist)
}
## ---------------------------------------------------------------------


## ---------------------------------------------------------------------
## ---------------Main function used to fit the model in stan
fit_stan_MDS <- function(
  mod,
  titer_map_train,
  titer_map_test,
  coord_prior_sd,
  n_dim, # integer
  chains = 3, # Number of MCMC chains to run
  cores = parallel::detectCores(logical = F), # For the cluster
  niter = 5000,
  diagnostic_file = NULL,
  debug = FALSE,
  return_diagnostic_plots = FALSE,
  ...
) {
  
  beepr::beep()
  ## Check inputs
  n_sera = unique(titer_map_train$serum) %>% length()
  n_antigens = unique(titer_map_train$antigen) %>% length()
  stopifnot(all(titer_map_test$serum %in% titer_map_train$serum)) # Verify that all antigens and sera in the test set appear at least once in the training set
  stopifnot(all(titer_map_test$antigen %in% titer_map_train$antigen))
  cat(sprintf('in R: n_antigens is %i; n_sera is %i \n', n_antigens, n_sera))
  
  ## Calculate distances and coordinate constraints
  distmat = distance_matrix_format(titer_map = titer_map_train)
  constrained_coords = get_constrainted_coords(distmat = distmat, n_dim = n_dim, verbose = TRUE)
  cat(sprintf('%s dim checkpoint0: generated inputs', n_dim))
  
  ## Load model
  model <- stan_model(mod)
  #(print(model))
  
  ## Format model inputs
  model_input_data <- list(
    N = nrow(titer_map_train),
    n_sera = n_sera,
    n_antigens = n_antigens,
    n_dim = n_dim,
    observed_titers = titer_map_train$logtiter, # N-vector of observed titers,
    smax = (titer_map_train$serum_potency+titer_map_train$antigen_avidity)/2,
    serum_id = titer_map_train$serum,
    antigen_id = titer_map_train$antigen,
    antigen_priors = constrained_coords,
    coord_prior_sd = coord_prior_sd,
    N_test_set = nrow(titer_map_test),
    smax_test_set = (titer_map_test$serum_potency + titer_map_test$antigen_avidity)/2,
    serum_id_test_set = titer_map_test$serum,
    antigen_id_test_set = titer_map_test$antigen,
    observed_titers_test_set = titer_map_test$logtiter
  )
  
  initlist <- lapply(1:chains, function(xx) initfun(distmat = distmat, 
                                                      n_dim = n_dim, 
                                                      n_antigens = n_antigens, 
                                                      n_sera = n_sera))
  cat(sprintf('%s dim checkpoint1: starting initial fit', n_dim))
  ## Run initial fit
  initial_fit <- sampling(model, 
                          data = model_input_data, 
                          chains = chains, 
                          init = initlist,
                          iter = niter,
                          cores = min(6, chains),
                          control = list(adapt_delta = 0.89,
                                         max_treedepth = 14),
                          diagnostic_file = diagnostic_file
  )
  
  ## Check if all 4 chains of the initial fit have mixed. If not, we'll 
  ##  re-run the chain, 
  ##  initializing with the median values of the highest likelihood chain
  return_initial_fit =   all(summary(initial_fit)$summary[,'Rhat'] <= 1.02) # logical -- have all parameters mixed?
  
  if(return_initial_fit == TRUE & return_diagnostic_plots!=TRUE){
    cat(sprintf('returning initial fit'))
    return(initial_fit)
  }
  
  cat(sprintf('%s dim checkpoint1: making diagnostic plots', n_dim))
  if(return_diagnostic_plots == TRUE | debug == TRUE){
  ## Make diagnostic plots
    ##  Sometimes chains converge to different local maxima, leading to high r-hat values for some inferred coordinates
    ##  To diagnose these situations, extract the coordinates wtih the worst R-hat values.
    ##  Plot their traceplots alongside violin plots of the corresponding log-posterior
    ##  We know that the likelihood has multiple local maxima, and that our model constraints only make the global 
    ##   maximum slightly better than a corresponding maximum reflected across the mth axis.
    ##  Below, if the chains haven't all mixed, we assume that one or more is stuck at the local (reflected) maximum
    ##  To solve this, we initialize a new run using values from the highest likelihood chain (this gets the algorithm off the lower-likelihood false peak)
    ##  The plots here allows us to check that some of the chains really do have higher log likelihoods than others 
    ##   (i.e. we use these plots to rule out other convergence issues)
  worst_4_coords <- as.data.frame(summary(initial_fit)$summary) %>% ## extract the coordinates with the worst rhat values
    as_tibble() %>%
    mutate(parname = names(initial_fit)) %>%
    arrange(-Rhat) %>%
    filter(grepl(pattern = 'coords', parname)) %>%
    slice(1:4) %>%
    pull(parname)
  # Traceplots for the worst 4 coordinates
  worst_traceplots =  traceplot(initial_fit, pars = worst_4_coords) +
    scale_color_viridis_d()+
    ggtitle('Initial fit traceplots-\n pars with worst Rhat') +
    theme(legend.position = 'bottom')
  # violin plot of log posterior by chain
  lp_plot =   as.array(initial_fit)[,,'lp__'] %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = 'chain', values_to = 'log_posterior') %>%
    ggplot() +
    geom_violin(aes(x = chain, y = log_posterior, color = chain), draw_quantiles = .5) +
    scale_color_viridis_d()+
    ggtitle('Initial log posteriors') +
    theme(legend.position = 'bottom')
  ## Combine into a single plot
  initial_diagnostic_plots = cowplot::plot_grid(
    worst_traceplots,
    lp_plot
  )
  
  if(debug == TRUE){ ## If debug 
    return(list(initial_fit = initial_fit, initial_diagnostic_plots = initial_diagnostic_plots))
  }
  if(return_initial_fit == TRUE){ ## If initial fit mixed and converged and we want to return diagnostic plots
    cat(sprintf('returning initial fit'))
    return(list(final_fit = initial_fit, initial_diagnostic_plots = initial_diagnostic_plots, final_diagnostic_plots = NULL))
  }
  }
  
  ## If the initial fit didn't converge, re-run, initializing with values from the best chain
  cat(print('Initial fit complete.\nRe-running using initial values from the best chain.\n'))
  
  ## Extract new initial conditions
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
           antigen_coords = matrix(best_chain_summary[is.antigen.coord,'50%'], 
                                   nrow = n_antigens, 
                                   ncol = n_dim,
                                   byrow = T),
           serum_coords = matrix(best_chain_summary[is.serum.coord,'50%'], 
                                 nrow = n_sera, 
                                 ncol = n_dim,
                                 byrow = T))
    }
    
    lapply(1:nchains, function(xx){get_one_list()})
  }
  
  ## Generate inits and refit
  inits <- initialize_with_best_chain(initial_fit, 1)
  refit <- sampling(
    model, model_input_data, 
    chains = 1, 
    init = inits,
    iter = niter,
    diagnostic_file = diagnostic_file)
  
  ## Check again that the new fit has mixed
  if(! all(summary(refit)$summary[,'Rhat'] <= 1.02)){
    cat(sprintf('Re-doing refit to achieve Rhat < 1.02'))
    refit <- sampling(
      model, model_input_data, 
      chains = 1, 
      init = inits,
      iter = niter,
      diagnostic_file = diagnostic_file)
  }
  
  ## Replicate the diagnostic plots from the initial fit, using the new fit, for comparison
  if(return_diagnostic_plots == TRUE){
    ## Replot the diagnostics for the new chain
  worst_traceplots =  traceplot(refit, pars = worst_4_coords)+
    scale_color_viridis_d()+
    ggtitle('Final traceplots') +
    theme(legend.position = 'bottom')
  lp_plot =   as.array(refit)[,,'lp__'] %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = 'chain', values_to = 'log_posterior') %>%
    ggplot() +
    geom_violin(aes(x = chain, y = log_posterior), draw_quantiles = .5)   +
    scale_color_viridis_d()+
    ggtitle('Final log posteriors')+
    theme(legend.position = 'bottom')
  
  refit_diagnostic_plots = cowplot::plot_grid(
    worst_traceplots,
    lp_plot
  )
  return(list(final_fit = refit, 
              initial_diagnostic_plots = initial_diagnostic_plots, 
              final_diagnostic_plots = refit_diagnostic_plots))
  }
  # else
  return(refit)
}
## ---------------------------------------------------------------------




# 
# plot_fits <- function(antigen_coords,
#                       fits, 
#                       n_iter,
#                       outdir){
#   trplot <- rstan::traceplot(fits, pars = names(fits))
#   
#   contour_posteriors <- extract_long_coords(fits) %>% 
#     filter(iter %% 5 == 0) %>% ## THIN
#     ggplot()+
#     geom_density_2d(aes(x = c1, y = c2, color = factor(chain))) +
#     facet_grid(kind ~ id)
#   
#   epitope_strain_map <- plot_inferred_original_map(antigen_coords, 
#                                                    fits)
#   
#   outdir = paste0(Sys.Date(), '/', outdir)
#   if(!dir.exists(paste0(Sys.Date()))) dir.create(paste0(Sys.Date()))
#   if(!dir.exists(outdir))  dir.create(outdir)
#   ggsave(paste0(outdir, '/traceplot.png'), plot = trplot)
#   ggsave(paste0(outdir, '/epitope_strain_map.png'), plot = epitope_strain_map)
#   ggsave(paste0(outdir, '/contour_posteriors.png'), plot = contour_posteriors)
#   
#   ## pairs plots
#   png(paste0(outdir, '/antigen_pairs.png'))
#   pairs(fits, pars = names(fits)[grepl(names(fits), pattern = 'antigen.+')])
#   dev.off()
#   png(paste0(outdir, '/serum_pairs.png'))
#   pairs(fits, pars = names(fits)[grepl(names(fits), pattern = 'serum.+')])
#   dev.off()
#   png(paste0(outdir, '/other_pairs.png'))
#   pairs(fits, pars = names(fits)[!grepl(names(fits), pattern = 'antigen.+') & !grepl(names(fits), pattern = 'serum.+')])
#   dev.off()
#   
#   
#   cat(sprintf('plots saved in %s', outdir))
#   return(epitope_strain_map)
# }



## ---------------------------------------------------------------------
## --------------- Extract estimated coordinates from the stan fit
## --------------- object 
extract_long_coords <- function(stan_fit){
  raw_fits <- rstan::extract(stan_fit, permuted = F, inc_warmup = F)
  niter <- dim(raw_fits)[1]
  nchains <- dim(raw_fits)[2]
  
  parnames = dimnames(raw_fits)[3]
  parname_contains <- function(parnames, this_pattern){sapply(parnames, function(xx) grepl(pattern = this_pattern, x = xx)) %>% as.vector()}
  is_antigen_coord = parname_contains(parnames, 'antigen_coords')
  is_serum_coord = parname_contains(parnames, 'serum_coords')
  is_log_posterior = parname_contains(parnames, 'lp_')

  
  #cat(print('checkpoint 1\n'))
  long_antigen_coords <- lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_antigen_coord]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    mutate(kind = 'antigen') %>%
    pivot_longer(contains('antigen_coords'), names_to = c('allele', 'coordinate'), names_pattern = 'antigen_coords.(\\d+),(\\d+).', values_to = 'value') %>%
    mutate(allele = as.numeric(allele) ) %>%
    pivot_wider(id_cols = c(chain, iter, allele, kind), names_from = coordinate, names_prefix = 'c', values_from = value) 
  
  #cat(print('checkpoint 2\n'))
  long_serum_coords <- lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_serum_coord]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    mutate(kind = 'serum') %>%
    pivot_longer(contains('serum_coords'), names_to = c('allele', 'coordinate'), names_pattern = 'serum_coords.(\\d+),(\\d+).', values_to = 'value') %>%
    mutate(allele = as.numeric(allele)) %>%
    pivot_wider(id_cols = c(chain, iter, allele, kind), names_from = coordinate, names_prefix = 'c', values_from = value) 
  
  n_allele <- max(long_serum_coords$allele)
  coord_names = long_serum_coords %>% select(matches('c\\d$')) %>% colnames()
  fill_vec = numeric(length(coord_names)+1)
  names(fill_vec) = c('iter', coord_names)
  fill_vec['iter'] = NA
  
  # (print(long_antigen_coords))
  # (print(long_serum_coords))
  # (print(ag2_long))
  
  long_coords <- bind_rows(long_antigen_coords,
                           long_serum_coords)  %>%
    complete(allele, kind, chain, fill = as.list(fill_vec)) ## Add 0 coordinates for the fixed antigen at the origin
  long_coords
}


## ---------------------------------------------------------------------
## --------------- summarize coordinates extracted from the stan fit 
## --------------- object --
## ---------------- this is a wrapper for extract_long_cords()
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
## ---------------------------------------------------------------------


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



# ## Test this tomorrow
# 
# infer_human_map <- function(ab_ag_df, # a long data frame of antigen coords and corresponding antibodies with an id column for distinct sera
#                             n_dim,
#                             n_chains = 4,
#                             n_iter = 5000){ 
#   
#   n_antigen = length(unique(ab_ag_df$antigen))
#   n_sera = length(unique(ag_ag_df$serum))
#   
#   ## Calculate the titer panel
#   ## Wrapper to calculate the titer for a given serum and antigen
#   
#   titer_wrapper <- function(this.serum, 
#                             this.antigen){
#     titer = with(extract_titer_inputs(antigen_id = this.antigen, 
#                                       serum_id = this.serum, 
#                                       merged_df = merged_df),
#                  solve_for_titer_multi_epitope(ab_position_list = ab_position_list, 
#                                                relative_concentrations = relative_concentrations, 
#                                                ag_list = ag_list, 
#                                                alpha = .25,
#                                                r = 7)
#     )
#     c('serum' = this.serum, 'antigen' = this.antigen, 'titer' = titer)
#   }
#   
#   titer_map <- with(expand.grid(antigen = 1:n_antigen, serum = 1:n_sera),
#                     mapply(FUN = titer_wrapper, 
#                            this.serum = serum,
#                            this.antigen = antigen)) %>% 
#     t() %>%
#     as_tibble() %>%
#     set_names(c('serum', 'antigen', 'titer')) %>%
#     mutate(logtiter = log2(titer/10)) %>%
#     group_by(serum) %>%
#     mutate(serum_potency = max(logtiter)) %>%
#     ungroup() %>% group_by(antigen) %>%
#     mutate(antigen_avidity = max(logtiter)) %>%
#     ungroup() %>%
#     mutate(titer_distance = (serum_potency+antigen_avidity)/2 - logtiter)
#   
#   ## Fit the stan model
#   fits <- fit_stan_MDS(observed_distances = format_stan_inputs(titer_map), 
#                        n_antigens = n_antigen, 
#                        n_sera = n_sera, 
#                        n_dim = n_dim, 
#                        chains = n_chains,
#                        niter = n_iter)
#   
#   
#   
#   return(list(ag_ab_coords = merged_df,
#               titer_map = titer_map,
#               summary_inferred_coords = extract_summary_coords(fits),
#               stan_fits = fits))
# }
# 
# 
# 
# plot_inferred_original_map2 <- function(stan_fits,
#                                         weights_for_shift, # An E-vector giving the weight of each epitope for calculating the centroids
#                                         antigen_coords,
#                                         flip_these_chains_over_x){
#   inferred_map <- extract_summary_coords(stan_fits) %>%
#     mutate(id = as.factor(id)) %>%
#     mutate(c2 = ifelse(chain %in% flip_these_chains_over_x, -c2, c2))
#   
#   ag1_ag2_centroids = get_weighted_centroids(antigen_coords, weights_for_shift)
#   
#   true_coords = standardize_coordinate_df(coord_df = bind_rows(ag1_ag2_centroids, antigen_coords), # Standardize to Ag1, Ag2 centroids
#                                           ag1_row = 1, ag2_row = 2) %>%
#     filter(!is.na(epitope)) %>% # Remove ag1, ag2 centroids
#     group_by(antigen) %>%
#     mutate(centroid1 = sum(c1*weights_for_shift),
#            centroid2 = sum(c2*weights_for_shift))
#   
#   ggplot() +
#     geom_polygon(aes(x = c1, y = c2, fill = antigen), data = true_coords, alpha = .05) +
#     geom_point(aes(x = centroid1, y = centroid2, color = antigen), data = true_coords, pch = 1) +
#     geom_point(aes(x = c1, y = c2, color = id, shape = kind), data = inferred_map)
# }
# 
# 

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
                                         observed_immunodominance_scheme = 'even' 
){
  
  
  
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
    geom_text(aes(x = component, y = cumulative_var_explained+.02, label = sprintf('%2.2f', cumulative_var_explained)))+
    geom_line(aes(x = as.numeric(component), y = cumulative_var_explained), color = 'deepskyblue') +
    ylim(c(0, 1.1))
  
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





fit_accross_dims <- function(n_dim_inputs,
                             immunodominance_flag){
  if(!dir.exists('diagnostic_plots')){dir.create('diagnostic_plots')}
  ## Fit maps across 1:8 dimensions
  inputs <- read_rds(sprintf('inputs/%sD_%s_immunodominance_inputs.rds', n_dim_inputs, immunodominance_flag))
  split_inputs <- test_train_split(df = inputs$titer_map)
  # Verify that all antigens and sera appear at least once in the training set
  stopifnot(all(unique(inputs$titer_map$antigen) %in% split_inputs$train$antigen))
  stopifnot(all(unique(inputs$titer_map$serum) %in% split_inputs$train$serum))
  fit_list <- vector(mode = 'list', length = 8)
  foreach(this_ndim = 1:8) %do% {
    cat(sprintf('Fitting %s dim model %s immunodominance\n', this_ndim, immunodominance_flag))
    this_fit <- fit_stan_MDS(mod = '../../stan/MDS_predict_titer_sparse.stan', 
                             titer_map_train = split_inputs$train, 
                             titer_map_test = split_inputs$test, 
                             coord_prior_sd = 0.1, 
                             n_dim = this_ndim, 
                             chains = 4, 
                             cores = 4, 
                             niter = 5000, 
                             diagnostic_file = NULL, 
                             debug = F, 
                             return_diagnostic_plots = T)
    cat(sprintf('%s dim fit_across_dims checkpoint1', this_ndim))
    ggsave(plot = cowplot::plot_grid(this_fit$initial_diagnostic_plots,
                                     this_fit$final_diagnostic_plots,
                                     nrow = 2),
           filename = sprintf('diagnostic_plots/%sDinput_%sDmap_%s.png', n_dim_inputs, this_ndim, immunodominance_flag),
           width = 9, height = 7, units = 'in', dpi = 200)
    this_fit$final_fit
  }
 # cat(sprintf('checkpoint1'))
  print(typeof(fit_list))
  names(fit_list) = paste0(1:8, 'D')
  if(!dir.exists('outputs')) dir.create('outputs')
  write_rds(fit_list, sprintf('outputs/%sDinputs-%s-fit_list.rds', n_dim_inputs, immunodominance_flag))
  write_rds(split_inputs, sprintf('outputs/%sDinputs-%s-test_train_split.rds', n_dim_inputs, immunodominance_flag))
 # cat(sprintf('checkpoint2'))
}
