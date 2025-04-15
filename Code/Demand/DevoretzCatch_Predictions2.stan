//
// This Stan program defines a demand model for wild salmon in BC as function
// of price and of competitive price (farmed salmon)
//


// The response data is a vector 'Q' for quantity of length 'N'
// the explanatory data is a vector 'P' of length N
data {
  int<lower=0> N;
  vector[N] newQ;
  vector[N] newD;
  
  int<lower=0> N_miss_Pf;        // Number of missing values
  int<lower=0> N_obs_Pf;        // Number of observed values
  array[N_obs_Pf] int<lower=1, upper=N> ii_obs_Pf;  // Position of observed values in the column
  array [N_miss_Pf] int<lower=1, upper=N> ii_mis_Pf; // Position of the missing values in the column
  vector[N_obs_Pf] Pf_obs;            // Observed values
  
  array[N] int<lower=1, upper=5> spp;  // array with species data
  int<lower=0> spp_n;     // number of species
 
  int<lower=0> N_miss_P;        // Number of missing values
 }

// We will fit the model in log space so we log the data
transformed data{
  vector[N] log_newQ;
  vector[N] log_newD;

  log_newQ = log(newQ);
  log_newD = log(newD);
}

// The parameters accepted by the model.
parameters {
  vector[spp_n] B0; // intercept
  vector[spp_n] B1; // coefficient on harvest
  vector[spp_n] B2; // coefficient on farmed salmon price
  vector[spp_n] B3; // coefficient on disposable income
  real<lower=0> sigma; // variance in observations

  //mean variance of parameters (hierarchical effect) - hyperparameters
  real mu_spp1;
  real<lower=0> sigma_spp1;
  real mu_spp2;
  real<lower=0> sigma_spp2;
  real mu_spp3;
  real<lower=0> sigma_spp3;
  
  //missing data parameters
  vector<lower=0>[N_miss_P] P_imputed;
  vector<lower=0>[N_miss_Pf] Pf_imputed;
}

// Transformed parameters to merge data with missing values and log for fitting
transformed parameters {
  vector[N] newPf;      
  vector[N] log_newPf;
  newPf[ii_obs_Pf] = Pf_obs;        
  newPf[ii_mis_Pf] = Pf_imputed;    
  log_newPf = log(newPf);
}

// Generate some predictions
generated quantities {
  vector[N] preds;
  vector[N] real_preds;
  vector[N] mu;      // logged mean
  vector[N] real_mu;      // mean
  mu = B0[spp] + B1[spp] .* log_newQ + B2[spp] .* log_newPf + B3[spp] .* log_newD;    // Linear model mean logged
  real_mu = exp(mu);
  
  for (i in 1:N) {
    preds[i] = (normal_rng(mu[i], sigma));
    real_preds[i] = exp(normal_rng(mu[i], sigma));
  }
}
