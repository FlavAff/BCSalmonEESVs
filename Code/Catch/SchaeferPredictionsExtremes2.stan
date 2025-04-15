//
  // This Stan program defines a model for the catch~effort relationship of commercial salmon fleets

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N_tot;
  int<lower=0> N_spp;
  int<lower=0> N_area;
  
  vector[N_tot] log_effortF;
  vector[N_tot] log_abundanceF;

  int<lower=0> N_a_mis; // missing abundance values
  
  array [N_tot] int<lower=1, upper=N_spp> spp;
  array [N_tot] int<lower=1, upper=N_area> area;
}


// The parameters accepted by the model. Our model
// accepts two parameters for regression and one for sd plus hyperparams
parameters {
  matrix[N_spp,N_area] log_q;
  vector<lower=0>[N_spp] sigma_obs;

  // Hyperparameters
  real mu_C;
  real<lower=0> sigma_C;
  
  // hurdle parameter
  matrix<lower=0, upper=1>[N_spp, N_area] theta;
  
  //Values to impute
  vector<lower=0>[N_a_mis] log_abundance_imputed;
}


generated quantities {
  vector[N_tot] preds;  // Simulated data based on the posterior distribution
  vector[N_tot] preds_real;  // Simulated data based on the posterior distribution
  vector[N_tot] mu;
  vector[N_tot] mu_real;
  
  for (n in 1:N_tot) {
     mu[n] = log_q[spp[n], area[n]] + log_effortF[n] + log_abundanceF[n];
  }
  mu_real = exp(mu);
       
  for (n in 1:N_tot) {
    if (bernoulli_rng(theta[spp[n], area[n]]) == 1) {
      // Simulating an excess zero
      preds[n] = 0;
    } else {
      // Simulating a non-zero value from the normal distribution
       preds[n] = cauchy_rng(mu[n], sigma_obs[spp[n]]);
    }
  }
  
  mu_real = exp(mu);
  preds_real = exp(preds);
}
