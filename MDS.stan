data {
  int<lower=0> n_strains;
  int<lower=0> n_sera;
  int<lower=0> n_dim;
  matrix[n_strains,n_sera] observed_distances;
}
parameters {
  real<lower=0> sigma;
  vector[n_dim] strain_coords[n_strains]; // Array of vectors, one for each strain, and each containing n_dim coords
  vector[n_dim] serum_coords[n_sera];
}
model {
  real estimated_distance; 

  for (serum in 1:n_sera){
  	for (strain in 2:n_strains){
      estimated_distance = distance(strain_coords[strain], serum_coords[serum]);
  		observed_distances[strain, serum] ~ normal(estimated_distance, sigma);
  	}
  }
 }


 