//
// This Stan program defines a linearised Ricker recruitment model, with a
// vector of values "S" spawners producing "R" recruits
//

data {
  int<lower=0> N;          // Total number of observations

  int<lower=0> N_obs_P;        // Number of observed values
  array[N_obs_P] int<lower=1, upper=N> ii_obs_P;  // Position of observed values in the column
  int<lower=0> N_obs_H;        // Number of observed values
  array[N_obs_H] int<lower=1, upper=N> ii_obs_H;  // Position of observed values in the column
  int<lower=0> N_obs_T;        // Number of observed values
  array[N_obs_T] int<lower=1, upper=N> ii_obs_T;  // Position of observed values in the column

  vector[N_obs_H] H_new_obs;            // Observed spawner values
  vector[N_obs_T] SST_new_obs;            // Observed spawner values
  vector[N_obs_P] log_P_new_obs;            // Observed parental values

  int<lower=0> N_miss_H;        // Number of missing values
  int<lower=0> N_miss_T;        // Number of missing values
  int<lower=0> N_miss_R;        // Number of missing values
  int<lower=0> N_miss_P;        // Number of missing values
  
  array[N_miss_P] int<lower=1, upper=N> ii_mis_P; // Position of missing values in log_P_new
  array[N_miss_T] int<lower=1, upper=N> ii_mis_T; // Position of missing values in SST_new
  array[N_miss_H] int<lower=1, upper=N> ii_mis_H; // Position of missing values in H_new
  
  int<lower=0> spp_n;     // number of regions
  array[N] int<lower=1, upper=spp_n> spp;  // array with region of origin data
}

parameters {
  vector<lower=0>[spp_n] alpha;          // productivity parameter
  vector<lower=0>[spp_n] beta;           // capacity parameter
  vector[spp_n] theta;          // temperature parameter
  vector[spp_n] delta;          // hatcheries parameter
  real<lower=0> log_sigma;     // Recruitment deviation 
  
  //hyper parameters
  real<lower=0> alpha_mu;
  real<lower=0> alpha_sigma;
  real<lower=0> beta_mu;
  
  vector[N_miss_P] log_P_imputed;
  vector[N_miss_R] log_R_imputed;
  vector[N_miss_T] SST_imputed;
  vector[N_miss_H] H_imputed;
}

transformed parameters {
  vector[N] log_P;
  log_P[ii_obs_P] = log_P_new_obs;
  log_P[ii_mis_P] = log_P_imputed;
  
  vector[N] SST;
  SST[ii_obs_T] = SST_new_obs;
  SST[ii_mis_T] = SST_imputed;
  
  vector[N] H;
  H[ii_obs_H] = H_new_obs;
  H[ii_mis_H] = H_imputed;
}

generated quantities {
  vector[N] log_preds;
  vector[N] preds;
  vector[N] mu;      // mean estimated recruits per spawner
  vector[N] log_mu;  // log of the estimate for normal fit
  log_mu = log_P + alpha[spp] - beta[spp] .* exp(log_P) + theta[spp] .* SST + delta[spp] .* H;    // Ricker model ++
  mu = exp(log_mu);

  for (i in 1:N) {
    log_preds[i] = normal_rng(log_mu[i], log_sigma);
    preds[i] = exp(normal_rng(log_mu[i], log_sigma));
  }
}
