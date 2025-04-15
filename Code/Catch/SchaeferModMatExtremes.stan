//
  // This Stan program defines a model for the catch~effort relationship of commercial salmon fleets

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N_tot;
  int<lower=0> N_spp;
  int<lower=0> N_area;
  
  vector[N_tot] log_effort;
  vector[N_tot] log_catch;
  
  int<lower=0> N_a_obs; // observed abundance values
  int<lower=0> N_a_mis; // missing abundance values
  array[N_a_obs] int<lower=1, upper=N_tot> ii_a_obs;  // Position of observed values in the column
  array [N_a_mis] int<lower=1, upper=N_tot> ii_a_mis; // Position of the missing values in the column
  vector[N_a_obs] log_abundance_obs;
  
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
  
  // Degrees of freedom for Student's t-distribution
  real<lower=1> df;  // df must be greater than 1 for the distribution to be valid
}

transformed parameters{
  vector[N_tot] log_abundance;      
  log_abundance[ii_a_obs] = log_abundance_obs;        // assign observations to the positions with observations
  log_abundance[ii_a_mis] = log_abundance_imputed;    // assign parameters (y missing) to the positions without observations
}

// The model to be estimated. We model the output
// to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  for (a in 1:N_spp) {
    for (g in 1:N_area) {
      log_q[a, g] ~ normal(mu_C,sigma_C);
      theta[a, g] ~ beta(3, 3);
    }
  }
  
  sigma_obs ~ exponential(1);
  
  mu_C ~ normal(0,1);
  sigma_C ~ exponential(1);
  
  log_abundance_imputed ~ normal(12,1.5);
  
  // Prior for df (degrees of freedom)
  df ~ gamma(2, 0.1);  // A weakly informative prior for df
  
  vector[N_tot] mu;
  for (n in 1:N_tot) {
    mu[n] = log_q[spp[n], area[n]] + log_effort[n] + log_abundance[n];
  }
  
  for (n in 1:N_tot){
    if (log_catch[n] == 0)
      target += log(theta[spp[n], area[n]]);
      else
        target += log1m(theta[spp[n], area[n]]) + student_t_lpdf(log_catch[n] | df, mu[n], sigma_obs[spp[n]]);
  }
}

generated quantities {
  vector[N_tot] preds;  // Simulated data based on the posterior distribution
  vector[N_tot] mu;
  vector[N_tot] log_lik; // Log-likelihood for each observation
  for (n in 1:N_tot) {
    mu[n] = log_q[spp[n], area[n]] + log_effort[n] + log_abundance[n];
  }
  
  for (n in 1:N_tot) {
    if (bernoulli_rng(theta[spp[n], area[n]]) == 1) {
      // Simulating an excess zero
      preds[n] = 0;
      log_lik[n] = log(theta[spp[n], area[n]]);
    } else {
      // Simulating a non-zero value from the normal distribution
      preds[n] = student_t_rng(df, mu[n], sigma_obs[spp[n]]);
      log_lik[n] = log1m(theta[spp[n], area[n]]) + student_t_lpdf(log_catch[n] | df, mu[n], sigma_obs[spp[n]]);
    }
  }
}
