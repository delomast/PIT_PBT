// attempting to implement the PIT tag estimator in Stan

functions{
  
}
// The input data is a vector 'y' of length 'N'.
data {
  int <lower = 1> N; // number of groups
  int <lower = 1> obs[N]; // number of detections for each group
  real <lower = 1e-20, upper = 1> tag[N]; // tag rates
}

transformed data{
  int div;
  div = 1;
}

// The parameters accepted by the model
parameters {
  // vector<lower = 0>[N] abundance;
  real log_abundance[N];
}

transformed parameters{
  real abund[N];
  for (g in 1:N){
    abund[g] = exp(log_abundance[g]);
  }
  
}
// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  // prior
  for (g in 1:N){
    log_abundance[g] ~ normal(0, 30);
  }
  // likelihood
  // int abund[N];
  for (g in 1:N){
     
    obs[g] ~ binomial(abund[g]/div, tag[g]);
  }
}

generated quantities{
  
}
