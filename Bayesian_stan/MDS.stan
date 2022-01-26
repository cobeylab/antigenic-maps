data {
  int<lower=0> n_antigens;
  int<lower=0> n_sera;
  int<lower=0> n_dim;
  //real sigma;
  real observed_distances[n_antigens,n_sera]; // Array of distances
}
parameters {
  real<lower=0> sigma; // sd of normal density around observed distances
  vector[n_dim] antigen_coords[n_antigens]; // Array of vectors, one for each strain, and each containing n_dim coords
  vector[n_dim] serum_coords[n_sera];
}
model {
  real estimated_distance; 
  sigma ~ normal(1.0, 1.0);

  for (serum in 1:n_sera){
  	for (strain in 1:n_antigens){
      estimated_distance = distance(antigen_coords[strain], serum_coords[serum]);
  		observed_distances[strain, serum] ~ normal(estimated_distance, sigma);
  	}
  }
 }


 