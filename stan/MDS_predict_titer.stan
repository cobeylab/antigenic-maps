data {
  int<lower=0> n_antigens;
  int<lower=0> n_sera;
  int<lower=0> n_dim;
  vector[n_dim] antigen_priors[n_dim+1]; // Array of priors used to constrain antigen coords 1:n_dim
  real observed_titers[n_antigens,n_sera]; // Array of observed titers
  real smax[n_antigens,n_sera]; // Used to calculate titer in the model
  real<lower=0> coord_prior_sd;
}

parameters {
  real<lower=0> sigma; // sd of normal density around observed distances
  vector[n_dim] antigen_coords[n_antigens]; // Array of vectors, one for each strain, and each containing n_dim coords
  vector[n_dim] serum_coords[n_sera];
}

transformed parameters{
  real map_distances[n_antigens, n_sera]; 
  real estimated_titers[n_antigens, n_sera];
  {
    for (serum in 1:n_sera){
      for (strain in 1:n_antigens){
        map_distances[strain,serum] = distance(antigen_coords[strain], serum_coords[serum]);
        estimated_titers[strain,serum] = smax[strain,serum] - map_distances[strain,serum];
      }
    }
  }
}

model {
  // prior on sigma
  sigma ~ normal(1.0, 1.0);
  // priors on the first n_dim+1 antigen coords to facilitate identifiability
  for(strain in 1:(n_dim+1)){
    for(coord in 1:n_dim){
      antigen_coords[strain][coord] ~ normal(antigen_priors[strain][coord], coord_prior_sd);
    }
  }

  for (serum in 1:n_sera){
    for (strain in 1:n_antigens){
      observed_titers[strain, serum] ~ normal(estimated_titers[strain,serum], sigma);
    }
  }
 }

 generated quantities{
    matrix[n_antigens, n_sera] log_lik;

    for(serum in 1:n_sera){
      for(strain in 1:n_antigens){
        log_lik[strain,serum] = normal_lpdf(observed_titers[strain,serum] | estimated_titers[strain,serum], sigma);
      }
    } 

 }
