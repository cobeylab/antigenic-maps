data {
  int<lower=0> n_points;
  int<lower=0> n_dim;
  matrix[n_points, n_points] observed_distances;
}
parameters {
  real<lower=0> sigma;
  vector[n_dim] coords[n_points]; // Array of n-dim coordinate vectors
}

model {
  matrix[n_points, n_points] estimated_distances;

  sigma~normal(1.0,1.0);

  for (point1 in 1:(n_points)){
  	for (point2 in (point1):n_points){
  		estimated_distances[point1, point2] = distance(coords[point1], coords[point2]);
  		observed_distances[point1, point2] ~ normal(estimated_distances[point1, point2], sigma);
  	}
  }
 }

