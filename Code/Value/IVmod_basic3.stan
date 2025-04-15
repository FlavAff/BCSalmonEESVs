//
// This Stan program defines a simple model, with two
// vectors of values, logged and modeled as normally distributed
// with means 'mu' and standard deviations 'sigma'.
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  vector[N] catch_;
  vector[N] harvest;
  vector[N] value;
  
  int<lower=0> N_spp; // Number of species
  array[N] int<lower=1,upper=N_spp> spp;
}

transformed data{
  vector[N] log_catch;
  vector[N] log_harvest;
  vector[N] log_value;
  
  log_catch = log(catch_);
  log_harvest = log(harvest);
  log_value = log(value);
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[N_spp] b0_h;
  vector[N_spp] b0_v;
  vector[N_spp] b1_h;
  vector[N_spp] b1_v;
  real<lower=0> sigma_v;
  real<lower=0> sigma_h;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  b0_h ~ normal(0,1);
  b0_v ~ normal(0,1);
  b1_h ~ normal(0,1);
  b1_v ~ normal(0,1);
  sigma_v ~ exponential(1);
  sigma_h ~ exponential(1);
  
  vector[N] mu_h;
  mu_h = b0_h[spp] + b1_h[spp] .* log_catch;
  vector[N] mu_v;
  mu_v = b0_v[spp] + b1_v[spp] .* log_harvest;
  
  log_harvest ~ normal(mu_h, sigma_h);
  log_value ~ normal(mu_v, sigma_v);
}

generated quantities {
  vector[N] harvest_pred;
  vector[N] value_pred;
  
  vector[N] mu_harvest;
  vector[N] mu_value;
  mu_harvest = b0_h[spp] + b1_h[spp] .* log_catch;
  mu_value = b0_v[spp] + b1_v[spp] .* log_harvest;

  for (i in 1:N) {
    harvest_pred[i] = normal_rng(mu_harvest[i], sigma_h);
    value_pred[i] = normal_rng(mu_value[i], sigma_v);
  }
}
