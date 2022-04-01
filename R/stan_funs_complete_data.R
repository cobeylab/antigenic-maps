distance_matrix_format <- function(titer_map){
  ## Extract titer distances from the titer map and format as a matrix
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


titer_matrix_format <- function(titer_map){
  ## Extract log titers from the titer map and format as a matrix
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


smax_matrix_format <- function(titer_map){
  ## Extract smax estimates from the titer map and format as a matrix
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


get_constrainted_coords <-   function(distmat, 
                                      n_dim,
                                      verbose = F) {
  ## Constrain the first few coordinates
  constrained_coords = matrix(0, nrow = n_dim+1, ncol = n_dim)
  constrained_coords[2,1] = mean(distmat[2,1], distmat[1,2], na.rm = T)
  
  if(n_dim >= 2){
    this_fun <- function(pars, 
                         this_dim){
      distances = colMeans(rbind(distmat[this_dim, 1:(this_dim-1)], 
                                 distmat[1:(this_dim-1), this_dim]), 
                           na.rm = T)
      target = 0
      for(ii in 1:(this_dim-1)){
        target = target + abs( sqrt(sum( (pars - constrained_coords[ii,1:(this_dim-1)])^2 )) - distances[ii] )
        if(verbose) cat(sprintf('target = %2.2f; this dist is %2.2f; this estimate is %2.2f\n', target, distances[ii], sqrt(sum( (pars - constrained_coords[ii,1:(this_dim-1)])^2 ) )))
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
  if(!all(abs(estimated_distances - distmat[1:(n_dim+1), 1:(n_dim+1)]) < 0.1)){
    warning(sprintf('max error when inferring constraints is > 0.1. Try rerunning as verbose.'))
  }
  
  print('Estimated distances are:\n')
  print(estimated_distances)
  print('Target distances are\n:')
  print( distmat[1:(n_dim+1),1:(n_dim+1)] )
  print('\n')
  
  constrained_coords
}

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



### FUNCTION TO FIT THE MODEL IN STAN
fit_stan_MDS <- function(
  mod, # path to stan model
  titer_map, ## data frame of antigen id, serum id, titer distance, titer, and smax
  n_dim, # integer
  chains = 4, # Number of MCMC chains to run. The function will require this many cores!
  niter = 5000,
  diagnostic_file = NULL,
  coord_prior_sd = 0.01, # sd of the prior used to constrain the first n_dim+1 antigen coordinates
  debug = FALSE,
  ...
) {
  library(rstan)
  
  n_antigens = length(unique(titer_map$antigen))
  n_sera = length(unique(titer_map$serum))
  
  observed_titers = titer_matrix_format(titer_map)
  observed_distances = distance_matrix_format(titer_map)
  smax = smax_matrix_format(titer_map)
  
  constrained_coords = get_constrainted_coords(distmat = observed_distances, 
                                               n_dim = n_dim)
  
  stopifnot(nrow(observed_titers) == n_antigens)
  stopifnot(ncol(observed_titers) == n_sera)
  stopifnot(nrow(smax) == n_antigens)
  stopifnot(ncol(smax) == n_sera)
  
  model <- stan_model(mod)
  (print(model))
  
  model_input_data <- list(
    n_antigens = n_antigens,
    n_sera = n_sera,
    n_dim = n_dim,
    antigen_priors = constrained_coords,
    observed_titers = observed_titers,
    smax = smax,
    coord_prior_sd = coord_prior_sd
  )
  
  initlist <- lapply(1:n_chains, function(xx) initfun(distmat = observed_distances, 
                                                      n_dim = n_dim, 
                                                      n_antigens = n_antigens, 
                                                      n_sera = n_sera))
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
  
  if(debug==TRUE){return(initial_fit)}
  
  if(all(summary(initial_fit)$summary[,'Rhat'] <= 1.02)){
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

  
  refit <- sampling(
    model, model_input_data, 
    chains = 1, 
    init = initialize_with_best_chain(initial_fit, 1),
    iter = niter,
    diagnostic_file = diagnostic_file)
  
  if(! all(summary(refit)$summary[,'Rhat'] <= 1.02)){
    cat(sprintf('Re-doing refit to achieve Rhat < 1.02'))
    refit <- sampling(
      model, model_input_data, 
      chains = 1, 
      init = initialize_with_best_chain(initial_fit, 1),
      iter = niter,
      diagnostic_file = diagnostic_file)
  }
  
  return(refit)
}














# 
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
# 
# 
# plot_antigen_coords <- function(antigen_coords){
#   antigen_coords %>%
#     ggplot() +
#     geom_line(aes(x=c1, y = c2, group = epitope), lty = 2, lwd = .1)+
#     geom_point(aes(x = c1, y = c2, shape = antigen, color = epitope))+
#     facet_wrap(.~epitope, labeller = label_both)
# }
# 
# 
# plot_Ab_Ag_map <- function(ab_ag_df,
#                            antigen_coords){
#   antisera_to_strain = unique(ab_ag_df$antigen)
#   ab_ag_df %>%
#     ggplot() +
#     geom_text(aes(x = c1, y = c2, color = antigen, label = epitope), data = antigen_coords, show.legend = F) + ## All epitope coords
#     geom_point(aes(x = c1_Ab, y = c2_Ab, color = antigen), pch = 3, alpha = .4) + # Ab coords
#     geom_text(aes(x = c1, y = c2, label = epitope), data = antigen_coords %>% filter(antigen == antisera_to_strain)) + # Highlight target antigen coords in black text and a circle
#    # geom_point(aes(x = c1, y = c2), pch = 1, data = antigen_coords %>% filter(antigen == antisera_to_strain)) +
#     guides(color=guide_legend(title="Abs to antigen"))+
#     xlab('c1')+
#     ylab('c2')+
#     ggtitle(sprintf('Antiserum to strain %s', antisera_to_strain))
# }
# 
# 

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
  
  # #cat(print('checkpoint 3\n'))
  # ag2_long <- lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,parname_contains(parnames, 'ag2')]) %>% mutate(iter = 1:niter)}) %>%
  #   bind_rows(.id = 'chain') %>%
  #   rename(c1 = value) %>%
  #   mutate(allele = 2) %>%
  #   mutate(kind = 'antigen') %>%
  #   select(chain, iter, allele, kind, c1)
  
  # ## If more than 3 dimensions, pad ag2 data frame with 0s
  # if(any(! (colnames(long_antigen_coords) %in% colnames(ag2_long)) )) {
  #   add_these_columns <- colnames(long_antigen_coords)[!(colnames(long_antigen_coords) %in% colnames(ag2_long))]
  #   for(cc in add_these_columns){
  #     ag2_long[[cc]] = 0
  #   }  
  # }
  
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
# 
# 
# flip_summary_over_x <- function(original_coords,
#                         inferred_summary){
#  flip.these.chains <-  merge(original_coords, 
#         inferred_summary %>%
#           filter(kind == 'antigen') %>%
#           rename(antigen = id),
#         by = 'antigen') %>%
#     group_by(chain) %>%
#     summarise(flipped_error = sum((c2.x + c2.y)^2, na.rm = T),
#               original_error = sum((c2.x - c2.y)^2, na.rm = T),
#               flip.me = flipped_error<original_error)
#  cat(sprintf('chain %s flipped over x\n', flip.these.chains %>% filter(flip.me) %>% pull(chain)))
#  print(flip.these.chains)
#  
#  inferred_summary %>% 
#    merge(flip.these.chains) %>%
#    mutate(across(starts_with("c2"), negate))
#   
# }
# 
# negate <- function(x){
#   -x
# }
# 
# 
# 
# 

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



get_pairwise_titer_error <- function(stan_fit,
                             titer_map){
  raw_fits <- rstan::extract(stan_fit, permuted = F, inc_warmup = F)
  niter <- dim(raw_fits)[1]
  nchains <- dim(raw_fits)[2]
  
  parnames = dimnames(raw_fits)[3]
  parname_contains <- function(parnames, this_pattern){sapply(parnames, function(xx) grepl(pattern = this_pattern, x = xx)) %>% as.vector()}
  is_map_distance = parname_contains(parnames, 'map')
  is_titer = parname_contains(parnames, 'estimated_titers')
  
merge(
  lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_map_distance]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    pivot_longer(contains('map'), names_to = c('antigen', 'serum'), names_pattern = 'map_distances.(\\d+),(\\d+).', values_to = 'map_distance'),
  lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_titer]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    pivot_longer(contains('titer'), names_to = c('antigen', 'serum'), names_pattern = 'estimated_titers.(\\d+),(\\d+).', values_to = 'estimated_titer')
  )%>%
  merge(titer_map, 
        by = c('serum', 'antigen')) %>%
  arrange(serum, antigen, chain, iter) %>%
  group_by(antigen, serum) %>%
  dplyr::summarise(mean_map_distance = mean(map_distance),
                   titer_distance = unique(titer_distance),
                   pairwise_distance_error = mean(sqrt( (map_distance-titer_distance)^2 )),
                   pairwise_titer_error = mean(sqrt((estimated_titer-logtiter)^2)),
                   mean_estimated_titer = mean(estimated_titer),
                   observed_titer = unique(logtiter)) %>%
    ungroup()
}

get_titer_error <- function(stan_fit, 
                          titer_map){
  get_pairwise_titer_error(stan_fit,
                   titer_map) %>%
    ungroup() %>%
    summarise(titer_error = sum(pairwise_titer_error)) %>%
    pull(titer_error)
}


get_distance_error <- function(stan_fit, 
                               titer_map){
  get_pairwise_titer_error(stan_fit,
                         titer_map) %>%
    ungroup() %>%
    summarise(distance_error = sum(pairwise_distance_error)) %>%
    pull(distance_error)
}



get_predictive_error <- function(stan_fit, 
                                 test_set){
  raw_fits <- rstan::extract(stan_fit, permuted = F, inc_warmup = F)
  niter <- dim(raw_fits)[1]
  nchains <- dim(raw_fits)[2]
  
  parnames = dimnames(raw_fits)[3]
  parname_contains <- function(parnames, this_pattern){sapply(parnames, function(xx) grepl(pattern = this_pattern, x = xx)) %>% as.vector()}
  is_predicted_titer = parname_contains(parnames, 'predicted_titers')
  is_predictive_error = parname_contains(parnames, 'posterior_predictive_titer_error')
  
merge(
  ## Reformatted predicted titers
  lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_predicted_titer]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    pivot_longer(contains('predicted'), names_to = c('id'), names_pattern = 'predicted_titers.(\\d+).', values_to = 'predicted_log_titer') %>%
    group_by(id),
  ## Reformatted predictive error
  lapply(1:nchains, function(xx){as.tibble(raw_fits[,xx,is_predictive_error]) %>% mutate(iter = 1:niter)}) %>%
    bind_rows(.id = 'chain') %>%
    pivot_longer(contains('predictive'), names_to = c('id'), names_pattern = 'posterior_predictive_titer_error.(\\d+).', values_to = 'predictive_error') %>%
    group_by(id)
)  %>%
    merge(test_set %>% mutate(id = 1:nrow(.))) %>%
    arrange(id, chain, iter) %>%
    ungroup()
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
  