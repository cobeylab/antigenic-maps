
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
  
  generate_initial_random_guess <- function(n_antigens, n_sera, n_dim){
    
    ag_coords <- runif(n_antigens*n_dim, -10, 10)
    names(ag_coords) = paste0('ag', rep(1:n_antigens, n_dim), 'c', rep(1:n_dim, each = n_antigens))
    
    ab_coords <- runif(n_sera*n_dim, -10, 10)
    names(ab_coords) = paste0('ab', rep(1:n_sera, n_dim), 'c', rep(1:n_dim, each = n_sera))
    
    c(ag_coords, ab_coords)
  }
  
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
                      .export = c('generate_initial_random_guess', 
                                  'get_ab_ag_distances',
                                  'get_euclidean_distance')) %dopar% {
                                    optim(par = generate_initial_random_guess(n_antigens, n_sera, n_dim),
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






