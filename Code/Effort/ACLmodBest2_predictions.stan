//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N_tot;

  int<lower=0> N_gear; // Number of gear types
  int<lower=0> N_area; // Number of gear types
  array[N_tot] int<lower=1,upper=N_gear> gear;     // gear data for all obs
  array[N_tot] int<lower=1,upper=N_area> area;     // area data for all obs
  
  int<lower=0> N_fmis; // missing fleet values
  int<lower=0> N_emis; // missing license values
  int<lower=0> N_lmis; // missing license values

  vector[N_tot] new_fleet;
  vector[N_tot] new_license;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[N_area] beta0_F;
  vector[N_area] beta1_F;
  matrix[N_area,N_gear] beta0_E;
  matrix[N_area,N_gear] beta1_E;
  real<lower=0> sigma_fleet;
  vector<lower=0>[N_gear] sigma_effort;
  
  //Hyper priors
  real mu_0F; // for fleet intercepts
  real mu_1F; // for fleet slopes
  real mu_0E; // for effort intercepts
  real mu_1E; // for effort slopes
  real<lower=0> sigma_0E;
  real<lower=0> sigma_1E;
  real<lower=0> sigma_0F;
  real<lower=0> sigma_1F;
  real<lower=0> sigma_E;
  
  vector<lower=0>[N_fmis] fleet_imputed;
  vector<lower=0>[N_lmis] license_imputed;
  vector<lower=0>[N_emis] effort_imputed;
}

// Generate some predictions
generated quantities {
  vector[N_tot] fleet_pred;
  vector[N_tot] effort_pred;
  
  vector[N_tot] mu_fleet;
  mu_fleet = beta0_F[area] + beta1_F[area] .* new_license;
  vector[N_tot] mu_effort;
  for (n in 1:N_tot) {
    mu_effort[n] = beta0_E[area[n], gear[n]] + beta1_E[area[n], gear[n]] .* new_fleet[n];
  }

  for (i in 1:N_tot) {
    fleet_pred[i] = normal_rng(mu_fleet[i], sigma_fleet);
    effort_pred[i] = normal_rng(mu_effort[i], sigma_effort[gear[i]]);
  }
}

