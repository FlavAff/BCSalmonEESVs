//
// This Stan program defines a linearised Ricker recruitment model, with a
// vector of values "S" spawners producing "R" recruits
//

data {
  int<lower=0> N;          // Total number of observations

  vector[N] H_new;            // Observed spawner values
  vector[N] SST_new;            // Observed spawner values
  vector[N] log_P_new;            // Observed parental values

  int<lower=0> N_miss_H;        // Number of missing values
  int<lower=0> N_miss_T;        // Number of missing values
  int<lower=0> N_miss_R;        // Number of missing values
  int<lower=0> N_miss_P;        // Number of missing values
  
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

generated quantities {
  vector[N] log_preds;
  vector[N] preds;
  vector[N] mu;      // mean estimated recruits per spawner
  vector[N] log_mu;  // log of the estimate for normal fit
  log_mu = log_P_new + alpha[spp] - beta[spp] .* exp(log_P_new) + theta[spp] .* SST_new + delta[spp] .* H_new;    // Ricker model ++
  mu = exp(log_mu);

  for (i in 1:N) {
    log_preds[i] = normal_rng(log_mu[i], log_sigma);
    preds[i] = exp(log_preds[i]);
  }
}
