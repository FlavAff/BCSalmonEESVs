//
// This Stan program generates predictions from a fitted standardized
// Ricker recruitment model. It takes new predictor data on the original
// scale, standardizes it using the means/SDs from the original data,
// and produces predictions on the original, un-standardized scale.
//

data {
  // --- New Data for Predictions ---
  int<lower=0> N;                              // Number of new observations for prediction
  vector[N] log_P_new;                         // New parental spawner data (log scale)
  vector[N] SST_new;                           // New sea surface temperature data
  vector[N] H_new;                             // New hatchery data
  int<lower=0> spp_n;                          // Number of species/regions
  array[N] int<lower=1, upper=spp_n> spp;      // Species/region for each new observation
  int<lower=0> N_miss_H;        // Number of missing values
  int<lower=0> N_miss_T;        // Number of missing values
  int<lower=0> N_miss_R;        // Number of missing values
  int<lower=0> N_miss_P;        // Number of missing values

  // --- Means and SDs from Original Model ---
  // These are required to standardize the new data and un-standardize the output
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
  // Standardize the new predictor data using the original means and SDs
  vector[N] log_P_new_std = (log_P_new - mean_log_P) / sd_log_P;
  vector[N] SST_new_std = (SST_new - mean_SST) / sd_SST;
  vector[N] H_new_std = (H_new - mean_H) / sd_H;
}

parameters {
  // This block must declare all the parameters from the original fitted model
  // so that the posterior draws can be accessed.

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

  // Imputed data parameters provided to run but unused
  vector[N_miss_P] log_P_imputed_std;
  vector[N_miss_R] log_R_imputed_std;
  vector[N_miss_T] SST_imputed_std;
  vector[N_miss_H] H_imputed_std;
}

transformed parameters {
  // Reconstruct the hierarchical parameters from their raw components,
  // same as in the original model.
  vector[spp_n] alpha;
  vector<lower=0>[spp_n] beta;
  alpha = mu_alpha + sigma_alpha * alpha_raw;
  beta = exp(mu_beta + sigma_beta * beta_raw);
}

model {
  // The model block is intentionally left empty.
  // We are not fitting the model, only using the parameters
  // from a previous fit to generate new quantities.
}

generated quantities {
  vector[N] log_mu;       // Mean estimate on original log-R scale
  vector[N] mu;       // Mean estimate on original R scale
  vector[N] log_R_preds;  // Predictions on original log-R scale
  vector[N] R_preds;      // Predictions on original R scale (exponentiated)

  { // Local block for efficiency
    vector[N] log_mu_std;   // Mean estimate on the standardized scale

    // 1. Calculate the mean prediction on the original log-R scale.
    //    This follows the exact formula from the original model's likelihood,
    //    using un-standardized log_P and standardized covariates.
    vector[N] log_mu_unstd = log_P_new + alpha[spp] - beta[spp] .* exp(log_P_new) + theta[spp] .* SST_new_std + delta[spp] .* H_new_std;

    // 2. Standardize the mean prediction to match the model's error scale.
    log_mu_std = (log_mu_unstd - mean_log_R) / sd_log_R;

    // Save the mean on the original log-scale for inspection
    log_mu = log_mu_unstd;
    mu = exp(log_mu);

    // 3. Generate predictions for each new data point
    for (i in 1:N) {
      // Generate a random draw on the standardized scale
      real pred_std = normal_rng(log_mu_std[i], sigma_std);

      // Un-standardize the prediction to get back to the original log_R scale
      log_R_preds[i] = (pred_std * sd_log_R) + mean_log_R;

      // Exponentiate to get the final prediction in terms of recruits (R)
      R_preds[i] = exp(log_R_preds[i]);
    }
  }
}
