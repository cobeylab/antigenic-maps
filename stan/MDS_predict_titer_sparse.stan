data {
  int<lower=1> N;
  int<lower=1> n_sera;
  int<lower=1> n_antigens;
  int<lower=1> n_dim;
  real<lower=0> observed_titers[N]; // Together these 4 columns are a long data frame
  real<lower=0> smax[N]; // used to calculate titer
  int<lower=1,upper=n_sera> serum_id[N]; // serum id -- used as index
  int<lower=1,upper=n_antigens> antigen_id[N]; // antigen id -- used as index
  // test set inputs for posterior prediction
  int<lower=1> N_test_set;
  real<lower=0> observed_titers_test_set[N_test_set]; // Together these 4 columns are a long data frame
  real<lower=0> smax_test_set[N_test_set]; // used to calculate titer
  int<lower=1,upper=n_sera> serum_id_test_set[N_test_set]; // serum id -- used as index
  int<lower=1,upper=n_antigens> antigen_id_test_set[N_test_set]; // antigen id -- used as index
}

parameters {
  // Note: ag1 coords are all 0 (ag1 falls at the origin)
  //       ag1 coord 1 is free. Ag2 coords 2, 3, etc. are 0 (ag2 is on the x-axis).
  real<lower=0> sigma; // sd of normal density around observed distances
  real<lower=0> ag2_c1; // x coordinate of antigen 2
  vector[n_dim] antigen_coords[n_antigens-2]; // Array of vectors, one for each strain, and each containing n_dim coords
  vector[n_dim] serum_coords[n_sera];
}

transformed parameters{
  real map_distances[N]; 
  real estimated_titers[N];
  {
    vector[n_dim] ag1_coords;
    vector[n_dim] ag2_coords;
    // set all ag1 coords to 0
    ag1_coords = rep_vector(0.0, n_dim);
    // ag2 c1 is a free par, and set all other coords to 0
    ag2_coords = rep_vector(0.0, n_dim);
    ag2_coords[1] = ag2_c1;

    for(ii in 1:N){
    // Declare temporary variables for the row-specific ag and antigen coords
      vector[n_dim] this_ag_coords;
      vector[n_dim] this_serum_coords;

    // set serum and ag coord values
      this_serum_coords = serum_coords[serum_id[ii]];
      if(antigen_id[ii]==1)
        this_ag_coords = ag1_coords;
      else if(antigen_id[ii]==2)
        this_ag_coords = ag2_coords;
      else
        this_ag_coords = antigen_coords[antigen_id[ii]-2];

    // Calculate the map distance and estimated titer for each row in the data
      map_distances[ii] = distance(this_ag_coords, this_serum_coords);
      estimated_titers[ii] = smax[ii] - map_distances[ii];
    }
  }
}

model {
  // prior on sigma
  sigma ~ normal(1.0, 1.0);

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
    vector[n_dim] ag2_coords;
    real this_map_distance;

    // set serum and ag coord values
    this_serum_coords = serum_coords[serum_id_test_set[ii]];

    if(antigen_id_test_set[ii]==1){
      this_ag_coords = rep_vector(0.0, n_dim);
    }
    else if(antigen_id_test_set[ii]==2){
      ag2_coords = rep_vector(0.0, n_dim);
      ag2_coords[1] = ag2_c1;
      this_ag_coords = ag2_coords;
    }
    else {
      this_ag_coords = antigen_coords[antigen_id_test_set[ii]-2];
    }

    // Calculate the map distance 
    this_map_distance = distance(this_ag_coords, this_serum_coords);

    // Generate posterior predictive draws of the estimated titers and error
    predicted_titers[ii] = normal_rng(smax_test_set[ii] - this_map_distance, sigma);
    posterior_predictive_titer_error[ii] = sqrt( (predicted_titers[ii]-observed_titers_test_set[ii])^2 ) ;
  }
}
