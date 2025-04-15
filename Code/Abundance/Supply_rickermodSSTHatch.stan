//
// This Stan program defines a linearised Ricker recruitment model, with a
// vector of values "S" spawners producing "R" recruits
//

data {
  int<lower=0> N;          // Total number of observations
 
  int<lower=0> N_miss_H;        // Number of missing values
  int<lower=0> N_obs_H;        // Number of observed values
  array[N_obs_H] int<lower=1, upper=N> ii_obs_H;  // Position of observed values in the column
  array [N_miss_H] int<lower=1, upper=N> ii_mis_H; // Position of the missing values in the column
  vector[N_obs_H] H_obs;            // Observed spawner values
 
  int<lower=0> N_miss_T;        // Number of missing values
  int<lower=0> N_obs_T;        // Number of observed values
  array[N_obs_T] int<lower=1, upper=N> ii_obs_T;  // Position of observed values in the column
  array [N_miss_T] int<lower=1, upper=N> ii_mis_T; // Position of the missing values in the column
  vector[N_obs_T] SST_obs;            // Observed spawner values
  
  int<lower=0> N_miss_R;        // Number of missing values
  int<lower=0> N_obs_R;        // Number of observed values
  array[N_obs_R] int<lower=1, upper=N> ii_obs_R;  // Position of observed values in the column
  array [N_miss_R] int<lower=1, upper=N> ii_mis_R; // Position of the missing values in the column
  vector[N_obs_R] log_R_obs;            // Observed spawner values
  
  int<lower=0> N_miss_P;        // Number of missing values
  int<lower=0> N_obs_P;        // Number of observed values
  array[N_obs_P] int<lower=1, upper=N> ii_obs_P;  // Position of observed values in the column
  array [N_miss_P] int<lower=1, upper=N> ii_mis_P; // Position of the missing values in the column
  vector[N_obs_P] log_P_obs;            // Observed parental values

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
  vector[N] log_R;      // create the dataset to fit the likelihood
  log_R[ii_obs_R] = log_R_obs;        // assign observations to the positions with observations
  log_R[ii_mis_R] = log_R_imputed;    // assign parameters (y missing) to the positions without observations
  
  vector[N] log_P;
  log_P[ii_obs_P] = log_P_obs;
  log_P[ii_mis_P] = log_P_imputed;
  
  vector[N] SST;
  SST[ii_obs_T] = SST_obs;
  SST[ii_mis_T] = SST_imputed;
  
  vector[N] H;
  H[ii_obs_H] = H_obs;
  H[ii_mis_H] = H_imputed;
}

model {
  // Priors
  alpha ~ gamma(alpha_mu,alpha_sigma); // needs to stay positive
  beta ~ exponential(beta_mu);  //beta itself is a small number, log beta can be 
  theta ~ normal(0,1);
  delta ~ normal(0,1);
  log_sigma ~ exponential(1);
  
  //Hyper parameters
  alpha_mu ~ gamma(5,5);
  alpha_sigma ~ exponential(1);
  beta_mu ~ exponential(2);
  
  log_P_imputed ~ normal(10,5);
  SST_imputed ~ normal(6,1);
  H_imputed ~ normal(12,4);
  
  vector[N] log_mu;  // log of the estimate for normal fit
  log_mu = log_P + alpha[spp] - beta[spp] .* exp(log_P) + theta[spp] .* SST + delta[spp] .* H;    // Ricker model ++
  
  // Likelihood for log of observed recruits with mean estimated from Ricker
  log_R ~ normal(log_mu, log_sigma);
}

generated quantities {
  vector[N] preds;
  //vector[N] mu;      // mean estimated recruits per spawner
  vector[N] mu_log;  // log of the estimate for normal fit
  mu_log = log_P + alpha[spp] - beta[spp] .* exp(log_P) + theta[spp] .* SST + delta[spp] .* H;    // Ricker model ++
  //mu = exp(mu);

  for (i in 1:N) {
    preds[i] = normal_rng(mu_log[i], log_sigma);
  }
}
