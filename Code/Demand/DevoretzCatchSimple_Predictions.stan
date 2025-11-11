//
// This Stan program defines a demand model for wild salmon in BC as function
// of price and of competitive price (farmed salmon)
//


// The response data is a vector 'Q' for quantity of length 'N'
// the explanatory data is a vector 'P' of length N
data {
  int<lower=0> N;
  vector[N] newQ;
  vector[N] newPf;
  array[N] int<lower=1, upper=5> spp;  // array with species data
  int<lower=0> spp_n;     // number of species
 
  int<lower=0> N_miss_P;        // Number of missing values
  int<lower=0> N_miss_Pf;        // Number of missing values
 }

// We will fit the model in log space so we log the data
transformed data{
  vector[N] log_newQ;
  vector[N] log_newPf;

  log_newQ = log(newQ);
  log_newPf = log(newPf);
  
  // *** CENTERING CALCULATIONS ***
  // Calculate means of log predictors
  real mean_log_Q = mean(log_newQ);
  real mean_log_Pf = mean(log_newPf);
  
  // Create centered log predictors
  vector[N] log_Q_c = log_newQ - mean_log_Q;
  vector[N] log_Pf_c = log_newPf - mean_log_Pf;
}

// The parameters accepted by the model.
parameters {
  vector[spp_n] B0; // intercept
  vector[spp_n] B1; // coefficient on harvest
  vector[spp_n] B2; // coefficient on farmed salmon price
  real<lower=0> sigma; // variance in observations

  //mean variance of parameters (hierarchical effect) - hyperparameters
  real mu_spp1;
  real<lower=0> sigma_spp1;
  real mu_spp2;
  real<lower=0> sigma_spp2;

  //missing data parameters
  vector<lower=0>[N_miss_P] P_imputed;
  vector<lower=0>[N_miss_Pf] Pf_imputed;
}


// Generate some predictions
generated quantities {
  vector[N] preds;
  vector[N] real_preds;
  vector[N] mu;      // logged mean
  vector[N] real_mu;      // mean
  mu = B0[spp] 
       + B1[spp] .* log_Q_c 
       + B2[spp] .* log_Pf_c;    // Linear model mean logged
  real_mu = exp(mu);
  
  for (i in 1:N) {
    preds[i] = (normal_rng(mu[i], sigma));
    real_preds[i] = exp(normal_rng(mu[i], sigma));
  }
}
