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
  
  int<lower=0> N_fobs; // observed fleet values
  int<lower=0> N_lobs; // observed license values
  int<lower=0> N_fmis; // missing fleet values
  int<lower=0> N_lmis; // missing license values
  int<lower=0> N_eobs; // missing fleet values
  int<lower=0> N_emis; // missing license values
  array[N_fobs] int<lower=1, upper=N_tot> ii_fobs;  // Position of observed values in the column
  array [N_fmis] int<lower=1, upper=N_tot> ii_fmis; // Position of the missing values in the column
  array[N_lobs] int<lower=1, upper=N_tot> ii_lobs;  // Position of observed values in the column
  array [N_lmis] int<lower=1, upper=N_tot> ii_lmis; // Position of the missing values in the column
  array[N_eobs] int<lower=1, upper=N_tot> ii_eobs;  // Position of observed values in the column
  array [N_emis] int<lower=1, upper=N_tot> ii_emis; // Position of the missing values in the column
  
  vector[N_fobs] fleet_obs;
  vector[N_lobs] license_obs;
  vector[N_eobs] effort_obs;
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

transformed parameters {
  vector<lower=0>[N_tot] fleet;      // create the dataset to fit the likelihood
  fleet[ii_fobs] = fleet_obs;        // assign observations to the positions with observations
  fleet[ii_fmis] = fleet_imputed;    // assign parameters (y missing) to the positions without observations
  
  vector<lower=0>[N_tot] licenses;      // create the dataset to fit the likelihood
  licenses[ii_lobs] = license_obs;        // assign observations to the positions with observations
  licenses[ii_lmis] = license_imputed;    // assign parameters (y missing) to the positions without observations
  
  vector<lower=0>[N_tot] effort;      // create the dataset to fit the likelihood
  effort[ii_eobs] = effort_obs;        // assign observations to the positions with observations
  effort[ii_emis] = effort_imputed;    // assign parameters (y missing) to the positions without observations
}

// The model to be estimated. We model the output
// 'fleet' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {

  beta0_F ~ normal(mu_0F,sigma_0F);
  beta1_F ~ normal(mu_1F,sigma_1F);
  for (a in 1:N_area) {
    for (g in 1:N_gear) {
      beta0_E[a, g] ~ normal(mu_0E,sigma_0E);
      beta1_E[a, g] ~ normal(mu_1E,sigma_1E);
    }
  }
  sigma_fleet ~ exponential(1);
  sigma_effort ~ exponential(sigma_E);
  
  // hyperpriors
  mu_0F ~ normal(2,5);
  mu_1F ~ normal(2,5);
  sigma_0F ~ exponential(1);
  sigma_1F ~ exponential(1);
  mu_0E ~ normal(10,5);
  mu_1E ~ normal(10,5);
  sigma_0E ~ exponential(1);
  sigma_1E ~ exponential(1);
  sigma_E ~ exponential(1);
    
  vector[N_tot] mu_fleet;
  mu_fleet = beta0_F[area] + beta1_F[area] .* licenses;
  vector[N_tot] mu_effort;
  for (n in 1:N_tot) {
    mu_effort[n] = beta0_E[area[n], gear[n]] + beta1_E[area[n], gear[n]] .* fleet[n];
  }

  licenses ~ lognormal(1,.5);
  fleet ~ normal(mu_fleet, sigma_fleet);
  effort ~ normal(mu_effort, sigma_effort[gear]);
}

generated quantities {
  vector[N_tot] fleet_pred;
  vector[N_tot] effort_pred;
  
  vector[N_tot] mu_fleet;
  mu_fleet = beta0_F[area] + beta1_F[area] .* licenses;
  vector[N_tot] mu_effort;
  for (n in 1:N_tot) {
    mu_effort[n] = beta0_E[area[n], gear[n]] + beta1_E[area[n], gear[n]] .* fleet[n];
  }

  for (i in 1:N_tot) {
    fleet_pred[i] = normal_rng(mu_fleet[i], sigma_fleet);
    effort_pred[i] = normal_rng(mu_effort[i], sigma_effort[gear[i]]);
  }
}
