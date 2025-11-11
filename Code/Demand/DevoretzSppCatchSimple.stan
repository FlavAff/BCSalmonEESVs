//
// This Stan program defines a demand model for wild salmon in BC as function
// of price and of competitive price (farmed salmon)
//


// The response data is a vector 'Q' for quantity of length 'N'
// the explanatory data is a vector 'P' of length N
data {
  int<lower=0> N;
  
  int<lower=0> N_miss_P;        // Number of missing values
  int<lower=0> N_obs_P;        // Number of observed values
  array[N_obs_P] int<lower=1, upper=N> ii_obs_P;  // Position of observed values in the column
  array [N_miss_P] int<lower=1, upper=N> ii_mis_P; // Position of the missing values in the column
  vector[N_obs_P] P_obs;            // Observed values
  
  int<lower=0> N_miss_Pf;        // Number of missing values
  int<lower=0> N_obs_Pf;        // Number of observed values
  array[N_obs_Pf] int<lower=1, upper=N> ii_obs_Pf;  // Position of observed values in the column
  array [N_miss_Pf] int<lower=1, upper=N> ii_mis_Pf; // Position of the missing values in the column
  vector[N_obs_Pf] Pf_obs;            // Observed values
  
  vector[N] Q;
  array[N] int<lower=1, upper=5> spp;  // array with species data
  int<lower=0> spp_n;     // number of species
}

// We will fit the model in log space so we log the data that is not missing
transformed data{
  vector[N] log_Q;
  log_Q = log(Q);
  
  // Calculate means of log predictors
  real mean_log_Q = mean(log_Q);
  
  // Create centered log predictors
  vector[N] log_Q_c = log_Q - mean_log_Q;
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

// Transformed parameters to merge data with missing values and log for fitting
transformed parameters {
  vector[N] P;      // create the dataset
  vector[N] log_P;      // create the log dataset to fit the likelihood
  P[ii_obs_P] = P_obs;        // assign observations to the positions with observations
  P[ii_mis_P] = P_imputed;    // assign parameters (y missing) to the positions without observations
  log_P = log(P);
  
  vector[N] Pf;      
  vector[N] log_Pf;
  Pf[ii_obs_Pf] = Pf_obs;        
  Pf[ii_mis_Pf] = Pf_imputed;    
  log_Pf = log(Pf);
  
  real mean_log_Pf = mean(log_Pf);
  vector[N] log_Pf_c = log_Pf - mean_log_Pf;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  //imputed data priors (non-negative)
  P_imputed ~ lognormal(2,2);
  Pf_imputed ~ normal(6,2);

  //Parameter fit priors
  B0 ~ normal(0,1);
  B1 ~ normal(mu_spp1,sigma_spp1);
  B2 ~ normal(mu_spp2,sigma_spp2);
  sigma ~ exponential(1);
  
  // Hyper priors
  mu_spp1 ~ normal(0,1);
  sigma_spp1 ~ exponential(.5);
  mu_spp2 ~ normal(0,1);
  sigma_spp2 ~ exponential(.5);
  
  // Mean term now uses the centered variables (log_Q_c, log_D_c, log_Pf_c)
  vector[N] mu;       // log mean
  mu = B0[spp] 
       + B1[spp] .* log_Q_c 
       + B2[spp] .* log_Pf_c; // Centered interaction 1

  //Log likelihood without autoregression
  log_P ~ normal(mu, sigma);
}

//Generate the average
generated quantities {
  vector[N] real_mu;       // log mean
  vector[N] preds;
  
  // The expected value (real_mu) calculation must use the centered variables
  real_mu = exp(B0[spp] 
                + B1[spp] .* log_Q_c 
                + B2[spp] .* log_Pf_c);  

  for (i in 1:N) {
    preds[i] = (normal_rng(log(real_mu[i]), sigma));
  }
 }
