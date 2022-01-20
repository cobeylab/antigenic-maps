data {
  int<lower=0> n_points;
  int<lower=0> n_dim;
  matrix[n_points, n_points] observed_distances;
}
parameters {
  real<lower=0> sigma;
  matrix[n_dim,n_points] coords;
}
model {
  matrix[n_points, n_points] estimated_distances;
  real l2_summand;


  for (point1 in 1:(n_points)){
  	for (point2 in (point1):n_points){
  	  	l2_summand = 0;
  		for(ii in 1:n_dim){
  			l2_summand += (coords[ii,point1]-coords[ii,point2])^2;
  		}
  		estimated_distances[point1, point2] = sqrt(l2_summand);
  		observed_distances[point1, point2] ~ normal(estimated_distances[point1, point2], sigma);
  	}
  }
 }