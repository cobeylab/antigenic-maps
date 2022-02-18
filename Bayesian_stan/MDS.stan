data {
  int<lower=3> n_antigens;
  int<lower=0> n_sera;
  int<lower=0> n_dim;
  //real sigma;
  real observed_distances[n_antigens,n_sera]; // Array of distances
}
parameters {
  // Note: ag1 coords are all 0 (ag1 falls at the origin)
  //       ag1 coord 1 is free. Ag2 coords 2, 3, etc. are 0 (ag2 is on the x-axis).
  real<lower=0> sigma; // sd of normal density around observed distances
  real<lower=0> ag2_c1; // x coordinate of antigen 2
  vector[n_dim] antigen_coords[n_antigens-2]; // Array of vectors, one for each strain, and each containing n_dim coords
  vector[n_dim] serum_coords[n_sera];
}
model {
  real estimated_distance; 
  vector[n_dim] ag1_coords;
  vector[n_dim] ag2_coords;
  // set all ag1 coords to 0
  ag1_coords = rep_vector(0.0, n_dim);
  // ag2 c1 is a free par, and set all other coords to 0
  ag2_coords = rep_vector(0.0, n_dim);
  ag2_coords[1] = ag2_c1;

  // prior on sigma
  sigma ~ normal(1.0, 1.0);



  for (serum in 1:n_sera){
  // ag1 
      estimated_distance = distance(ag1_coords, serum_coords[serum]);
      observed_distances[1, serum] ~ normal(estimated_distance, sigma);
  // ag2
      estimated_distance = distance(ag2_coords, serum_coords[serum]);
      observed_distances[2, serum] ~ normal(estimated_distance, sigma);
  	for (strain in 3:n_antigens){
      estimated_distance = distance(antigen_coords[strain-2], serum_coords[serum]);
  		observed_distances[strain, serum] ~ normal(estimated_distance, sigma);
  	}
  }
 }

