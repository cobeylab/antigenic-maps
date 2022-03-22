data {
  int<lower=3> n_antigens;
  int<lower=0> n_sera;
  int<lower=0> n_dim;
  real observed_titers[n_antigens,n_sera]; // Array of observed titers
  real smax[n_antigens,n_sera]; // Used to calculate titer in the model
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
  real map_distances[n_antigens, n_sera]; 
  real estimated_titers[n_antigens, n_sera];
  {
    vector[n_dim] ag1_coords;
    vector[n_dim] ag2_coords;
    // set all ag1 coords to 0
    ag1_coords = rep_vector(0.0, n_dim);
    // ag2 c1 is a free par, and set all other coords to 0
    ag2_coords = rep_vector(0.0, n_dim);
    ag2_coords[1] = ag2_c1;

    for (serum in 1:n_sera){
      // ag1 
      map_distances[1,serum] = distance(ag1_coords, serum_coords[serum]);
      estimated_titers[1,serum] = smax[1,serum]-map_distances[1,serum];
      // ag2
      map_distances[2,serum] = distance(ag2_coords, serum_coords[serum]);
      estimated_titers[2,serum] = smax[2,serum] - map_distances[2,serum];
      for (strain in 3:n_antigens){
        map_distances[strain,serum] = distance(antigen_coords[strain-2], serum_coords[serum]);
        estimated_titers[strain,serum] = smax[strain,serum] - map_distances[strain,serum];
      }
    }
  }
}

model {
  // prior on sigma
  sigma ~ normal(1.0, 1.0);

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
