//
// This Stan program defines a linearised Ricker recruitment model.
// This version uses a standardized response variable (log_R) and
// standardized predictors.
//

data {
  int<lower=0> N;       // Total number of observations

  // Data for H (hatchery)
  int<lower=0> N_miss_H;
  int<lower=0> N_obs_H;
  array[N_obs_H] int<lower=1, upper=N> ii_obs_H;
  array[N_miss_H] int<lower=1, upper=N> ii_mis_H;
  vector[N_obs_H] H_new_obs;

  // Data for SST (temperature)
  int<lower=0> N_miss_T;
  int<lower=0> N_obs_T;
  array[N_obs_T] int<lower=1, upper=N> ii_obs_T;
  array[N_miss_T] int<lower=1, upper=N> ii_mis_T;
  vector[N_obs_T] SST_new_obs;

  // Data for log_P (log parental spawners)
  int<lower=0> N_miss_P;
  int<lower=0> N_obs_P;
  array[N_obs_P] int<lower=1, upper=N> ii_obs_P;
  array[N_miss_P] int<lower=1, upper=N> ii_mis_P;
  vector[N_obs_P] log_P_new_obs;
  
  // N_miss_R is needed to declare log_R_imputed_std below
  int<lower=0> N_miss_R;

  // Species / region data
  int<lower=0> spp_n; // number of species/regions
  array[N] int<lower=1, upper=spp_n> spp;

  // Means and SDs for standardization
  real mean_log_P;
  real sd_log_P;
  real mean_SST;
  real sd_SST;
  real mean_H;
  real sd_H;
  real mean_log_R; 
  real sd_log_R;   
}

transformed data {
  // Standardize observed data once at the start
  vector[N_obs_P] log_P_new_obs_std = (log_P_new_obs - mean_log_P) / sd_log_P;
  vector[N_obs_T] SST_new_obs_std = (SST_new_obs - mean_SST) / sd_SST;
  vector[N_obs_H] H_new_obs_std = (H_new_obs - mean_H) / sd_H;
}

parameters {
  // Species-specific parameters
  vector[spp_n] alpha_raw;
  vector[spp_n] beta_raw;
  vector[spp_n] theta;
  vector[spp_n] delta;
  real<lower=0> sigma_std; // Error term is on the standardized scale

  // Hyperparameters
  real mu_alpha;
  real<lower=0> sigma_alpha;
  real mu_beta;
  real<lower=0> sigma_beta;

  // Imputed data (all on standardized scale)
  vector[N_miss_P] log_P_imputed_std;
  vector[N_miss_R] log_R_imputed_std;
  vector[N_miss_T] SST_imputed_std;
  vector[N_miss_H] H_imputed_std;
}

transformed parameters {
  // Hierarchical params
  vector[spp_n] alpha;
  vector<lower=0>[spp_n] beta;
  alpha = mu_alpha + sigma_alpha * alpha_raw;
  beta = exp(mu_beta + sigma_beta * beta_raw);

  // --- Assemble full standardized vectors ---

  vector[N] log_P_std;
  log_P_std[ii_obs_P] = log_P_new_obs_std;
  log_P_std[ii_mis_P] = log_P_imputed_std;

  vector[N] SST_std;
  SST_std[ii_obs_T] = SST_new_obs_std;
  SST_std[ii_mis_T] = SST_imputed_std;

  vector[N] H_std;
  H_std[ii_obs_H] = H_new_obs_std;
  H_std[ii_mis_H] = H_imputed_std;
}

model {
}

generated quantities {
  vector[N] log_mu_unstd; // Mean estimate on original log-recruits scale
  vector[N] log_mu; 
  vector[N] mu; 
  vector[N] log_preds;    // Predictions on original log-recruits scale
  vector[N] preds;        // Predictions on original recruits scale (exp)

  { // Local block for efficiency
    // Un-standardize log_P to use in the non-linear part of the Ricker model
    vector[N] log_P = (log_P_std * sd_log_P) + mean_log_P;
    
    // Calculate the mean of log-recruits on the unstandardized scale
    log_mu_unstd = log_P + alpha[spp] - beta[spp] .* exp(log_P) + theta[spp] .* SST_std + delta[spp] .* H_std;
    
    // Save the means
    log_mu = log_mu_unstd;
    mu = exp(log_mu);
    
    // Standardize the mean prediction to use with the standardized error term
    vector[N] log_mu_std = (log_mu_unstd - mean_log_R) / sd_log_R;

    for (i in 1:N) {
      // Generate a random draw on the standardized scale
      real pred_std = normal_rng(log_mu_std[i], sigma_std);
      
      // Un-standardize the prediction to the original log-recruits scale
      log_preds[i] = (pred_std * sd_log_R) + mean_log_R;
    }
  }
  
  // Exponentiate to get predictions on the natural scale of recruits
  preds = exp(log_preds);
}
