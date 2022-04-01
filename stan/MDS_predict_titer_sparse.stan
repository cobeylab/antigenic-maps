data {
  int<lower=1> N;
  int<lower=1> n_sera;
  int<lower=1> n_antigens;
  int<lower=1> n_dim;
  real observed_titers[N]; // Together these 4 columns are a long data frame
  real<lower=0> smax[N]; // used to calculate titer
  int<lower=1,upper=n_sera> serum_id[N]; // serum id -- used as index
  int<lower=1,upper=n_antigens> antigen_id[N]; // antigen id -- used as index
  // inputs to constrain antigen coordinates in the lower triangle of the coordinate matrix 1:(n_dim+1)
  vector[n_dim] antigen_priors[n_dim+1]; // Array of priors used to constrain antigen coords 1:n_dim
  real<lower=0> coord_prior_sd; // prior variance for constrained coordinates
  // test set inputs for posterior prediction
  int<lower=1> N_test_set;
  real observed_titers_test_set[N_test_set]; // Together these 4 columns are a long data frame
  real<lower=0> smax_test_set[N_test_set]; // used to calculate titer
  int<lower=1,upper=n_sera> serum_id_test_set[N_test_set]; // serum id -- used as index
  int<lower=1,upper=n_antigens> antigen_id_test_set[N_test_set]; // antigen id -- used as index
}

parameters {
  real<lower=0> sigma; // sd of normal density around observed distances
  vector[n_dim] antigen_coords[n_antigens]; // Array of vectors, one for each strain, and each containing n_dim coords
  vector[n_dim] serum_coords[n_sera];
}

transformed parameters{
  real map_distances[N]; 
  real estimated_titers[N];

  for(ii in 1:N){
    // Declare temporary variables for the row-specific ag and antigen coords
    vector[n_dim] this_ag_coords;
    vector[n_dim] this_serum_coords;

    // set serum and ag coord values
    this_serum_coords = serum_coords[serum_id[ii]];
    this_ag_coords = antigen_coords[antigen_id[ii]];

    // Calculate the map distance and estimated titer for each row in the data
    map_distances[ii] = distance(this_ag_coords, this_serum_coords);
    estimated_titers[ii] = smax[ii] - map_distances[ii];
  }
}

model {
  // prior on sigma
  sigma ~ normal(1.0, 1.0);
  for(strain in 1:(n_dim+1)){
    for(coord in 1:n_dim){
      antigen_coords[strain][coord] ~ normal(antigen_priors[strain][coord], coord_prior_sd);
    }
  }

  for(ii in 1:N){
    observed_titers[ii] ~ normal(estimated_titers[ii], sigma);
  }
 }

generated quantities{
  real log_lik[N];
  real predicted_titers[N_test_set];
  real posterior_predictive_titer_error[N_test_set];

  // Compute the log likelihood at each iteration so that we can later calculate the loo
  for(ii in 1:N){
    log_lik[ii] = normal_lpdf(observed_titers[ii] | estimated_titers[ii], sigma);
  } 

  // Compute the posterior predictive distribution, and error

  for(ii in 1:N_test_set){
    vector[n_dim] this_ag_coords;
    vector[n_dim] this_serum_coords;
    real this_map_distance;

    // set serum and ag coord values
    this_serum_coords = serum_coords[serum_id_test_set[ii]];
    this_ag_coords = antigen_coords[antigen_id_test_set[ii]];

    // Calculate the map distance 
    this_map_distance = distance(this_ag_coords, this_serum_coords);
    //print("ag coords are", this_ag_coords);
    //print("serum coords are", this_serum_coords);
    //print("this distance is", this_map_distance);

    // Generate posterior predictive draws of the estimated titers and error
    predicted_titers[ii] = normal_rng(smax_test_set[ii] - this_map_distance, sigma);
    //print("ssmax is", smax_test_set[ii]);
    //print("predicted titer is", predicted_titers[ii]);
    posterior_predictive_titer_error[ii] = sqrt( (predicted_titers[ii]-observed_titers_test_set[ii])^2 ) ;
  }
}
