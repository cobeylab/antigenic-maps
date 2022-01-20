data {
  int<lower=0> n_strains;
  int<lower=0> n_sera;
  int<lower=0> n_dim;
  matrix[n_strains,n_sera] observed_distances;
  vector[n_dim] first_strain_coords; // Fix the (arbitrary) coordinates of one strain
}
parameters {
  real<lower=0> sigma;
  vector[n_dim] strain_coords[n_strains-1]; // Array of vectors, one for each strain, and each containing n_dim coords
  vector[n_dim] serum_coords[n_sera];
}
model {
  real estimated_distance; 

  for (serum in 1:n_sera){
    // Consider the fixed strain and free serum coords
      estimated_distance = distance(first_strain_coords, serum_coords[serum]);
      observed_distances[1, serum] ~ normal(estimated_distance, sigma);
  // Loop through free strains
  	for (strain in 2:n_strains){
      estimated_distance = distance(strain_coords[strain-1], serum_coords[serum]);
  		observed_distances[strain, serum] ~ normal(estimated_distance, sigma);
  	}
  }
 }


 