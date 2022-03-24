## Least_squares

This directory implements the approach of [Smith et al. 2004]() to infer antigenic maps using a frequentist approach.

The core function, `fit_MDS_least_squares` is defined in `../code.R`.
This directory includes notebooks that demonstrate the method's use.

The approach uses conjugate gradient optimization, with multipe randomly drawn initial states, to minimize the error function:

<img src="https://latex.codecogs.com/svg.image?\sum_{kj}(d_{kj}-\delta_{kj})^2" title="\sum_{kj}(d_{kj}-\delta_{kj})^2" />

Where d<sub>kj</sub> is the observed titer distance between antiserum k and antigen j, and <img src="https://latex.codecogs.com/svg.image?\delta_{kj}" title="\delta_{kj}" /> is the map distance.

The titer distance is defined as <img src="https://latex.codecogs.com/svg.image?d_{kj}&space;=&space;s^{max}_j&space;-&space;s_{kj}" title="d_{kj} = s^{max}_j - s_{kj}" />, where <img src="https://latex.codecogs.com/svg.image?s^{max}_j" title="s^{max}_j" /> is the maximum log2 titer between serum k and any strain, and <img src="https://latex.codecogs.com/svg.image?s^{max}_j" title="s_{kj}" /> is the observed log2 titer between serum k and strain j.

The map distance is defined as <img src="https://latex.codecogs.com/svg.image?\delta_{kj}&space;=&space;||&space;x_k&space;-&space;y_j&space;||_2" title="\delta_{kj} = || x_k - y_j ||_2" />, where x<sub>k</sub> are the inferred coordinates of serum k, y<sub>j</sub> are the inferred map coordinates of strain j, and ||.||<sub>2</sub> is the Euclidean norm. 
