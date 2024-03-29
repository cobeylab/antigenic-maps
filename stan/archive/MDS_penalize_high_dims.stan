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

transformed parameters{
  real map_distances[n_antigens, n_sera]; 
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
      // ag2
      map_distances[2,serum] = distance(ag2_coords, serum_coords[serum]);
      for (strain in 3:n_antigens){
        map_distances[strain,serum] = distance(antigen_coords[strain-2], serum_coords[serum]);
      }
    }
  }
}

model {
  // prior on sigma
  sigma ~ normal(1.0, 1.0);
  
  // priors on coordinate values
  // impose a prior of 0, with variance decreasing in higher dimensions (a higher penalty for non-0 params in higher dims)
  for(strain in 1:(n_antigens-2)){
    for(dim in 1:n_dim){
      antigen_coords[strain][dim] ~ normal(0.0, 50.0/dim);
    }
  }
  for(serum in 1:n_sera){
    for(dim in 1:n_dim){
      serum_coords[serum][dim] ~ normal(0.0, 50.0/dim);
    }
  }

  for (serum in 1:n_sera){
    for (strain in 1:n_antigens){
      observed_distances[strain, serum] ~ normal(map_distances[strain,serum], sigma);
    }
  }
 }
