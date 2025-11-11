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
  vector[N_obs_H] H_obs;

  // Data for SST (temperature)
  int<lower=0> N_miss_T;
  int<lower=0> N_obs_T;
  array[N_obs_T] int<lower=1, upper=N> ii_obs_T;
  array[N_miss_T] int<lower=1, upper=N> ii_mis_T;
  vector[N_obs_T] SST_obs;

  // Data for log_R (log recruits)
  int<lower=0> N_miss_R;
  int<lower=0> N_obs_R;
  array[N_obs_R] int<lower=1, upper=N> ii_obs_R;
  array[N_miss_R] int<lower=1, upper=N> ii_mis_R;
  vector[N_obs_R] log_R_obs;

  // Data for log_P (log parental spawners)
  int<lower=0> N_miss_P;
  int<lower=0> N_obs_P;
  array[N_obs_P] int<lower=1, upper=N> ii_obs_P;
  array[N_miss_P] int<lower=1, upper=N> ii_mis_P;
  vector[N_obs_P] log_P_obs;

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
  vector[N_obs_P] log_P_obs_std = (log_P_obs - mean_log_P) / sd_log_P;
  vector[N_obs_T] SST_obs_std = (SST_obs - mean_SST) / sd_SST;
  vector[N_obs_H] H_obs_std = (H_obs - mean_H) / sd_H;
  vector[N_obs_R] log_R_obs_std = (log_R_obs - mean_log_R) / sd_log_R;
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
  vector[N] log_R_std;
  log_R_std[ii_obs_R] = log_R_obs_std;
  log_R_std[ii_mis_R] = log_R_imputed_std;

  vector[N] log_P_std;
  log_P_std[ii_obs_P] = log_P_obs_std;
  log_P_std[ii_mis_P] = log_P_imputed_std;

  vector[N] SST_std;
  SST_std[ii_obs_T] = SST_obs_std;
  SST_std[ii_mis_T] = SST_imputed_std;

  vector[N] H_std;
  H_std[ii_obs_H] = H_obs_std;
  H_std[ii_mis_H] = H_imputed_std;
}

model {
  // Priors
  theta ~ normal(0, 1);
  delta ~ normal(0, 1);
  sigma_std ~ exponential(1); // Prior on standardized sigma

  // Hyperpriors for hierarchical parameters
  // Note: mu_alpha prior is centered on 0, as the mean of the
  // response is now 0.
  mu_alpha ~ normal(0, 1);
  sigma_alpha ~ exponential(1);
  mu_beta ~ normal(0, 1);
  sigma_beta ~ exponential(1);

  // Priors for species-specific parameters
  alpha_raw ~ std_normal();
  beta_raw ~ std_normal();

  // Priors for all imputed values (now standardized)
  log_R_imputed_std ~ std_normal();
  log_P_imputed_std ~ std_normal();
  SST_imputed_std ~ std_normal();
  H_imputed_std ~ std_normal();

  // Model likelihood
  vector[N] log_mu_std;
  { // Local block for efficiency
    vector[N] log_P = (log_P_std * sd_log_P) + mean_log_P; // un-standardize P
    vector[N] log_mu_unstd = log_P + alpha[spp] - beta[spp] .* exp(log_P) + theta[spp] .* SST_std + delta[spp] .* H_std;
    // Standardize the mean prediction
    log_mu_std = (log_mu_unstd - mean_log_R) / sd_log_R;
  }

  log_R_std ~ normal(log_mu_std, sigma_std);
}

generated quantities {
  vector[N] log_mu;     // Mean estimate on original scale
  vector[N] preds;      // Predictions on original scale
  vector[N] log_mu_std; // Mean estimate on standardized scale
  
  { // Local block for efficiency
    vector[N] log_P = (log_P_std * sd_log_P) + mean_log_P;
    vector[N] log_mu_unstd = log_P + alpha[spp] - beta[spp] .* exp(log_P) + theta[spp] .* SST_std + delta[spp] .* H_std;
    log_mu_std = (log_mu_unstd - mean_log_R) / sd_log_R;
  }

  // Un-standardize the mean and the predictions
  log_mu = (log_mu_std * sd_log_R) + mean_log_R;
  
  for (i in 1:N) {
    // Generate prediction on standardized scale, then transform back
    real pred_std = normal_rng(log_mu_std[i], sigma_std);
    preds[i] = (pred_std * sd_log_R) + mean_log_R;
  }
}
