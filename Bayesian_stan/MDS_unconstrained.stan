data {
  int<lower=1> n_antigens;
  int<lower=1> n_sera;
  int<lower=1> n_dim;
  //real sigma;
  real observed_distances[n_antigens,n_sera]; // Array of distances
}
parameters {
  real<lower=0> sigma; // sd of normal density around observed distances
  vector[n_dim] antigen_coords[n_antigens]; // Array of vectors, one for each strain, and each containing n_dim coords
  vector[n_dim] serum_coords[n_sera];
}

transformed parameters{
  real map_distances[n_antigens, n_sera]; 
  {
    for (serum in 1:n_sera){
      for (strain in 1:n_antigens){
        map_distances[strain,serum] = distance(antigen_coords[strain], serum_coords[serum]);
      }
    }
  }
}

model {
  // assume c1 of Ag1 is 0 to help make the map identifiable
  antigen_coords[1][1] ~ normal(0.0, 0.1); // Fig ag1 coordinates 1 and 2 to the origin

  if(n_dim > 1){
    antigen_coords[1][2] ~ normal(0.0, 0.1);
    antigen_coords[2][2] ~ normal(0.0, 0.1); // If 2D+, fix ag2 to the x axis
    }

  // Prior on sigma --
  // Assume sigma is relatively low. Else sigma grows to infinity to make the observation model permissible.
  sigma ~ normal(1.0, 1.0);

  for (serum in 1:n_sera){
    for (strain in 1:n_antigens){
      observed_distances[strain, serum] ~ normal(map_distances[strain,serum], sigma);
    }
  }
 }
