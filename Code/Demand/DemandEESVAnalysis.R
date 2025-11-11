setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(cmdstanr)
library(ggplot2)
library(rethinking)
library(tidyverse)
library(tidybayes)
library(shinystan)
library(rstan)
library(posterior)
library(loo)

source("DemandModels/DemandData.R")

#Do the modelling
Demand_mod <- cmdstan_model(stan_file = "DemandModels/DevoretzSppCatchSimple.stan", pedantic=T) #not as good .88/.89 but better as data less related
datalist <- list(N = nrow(dat),
                 N_years = length(unique(dat$year)),
                 D = dat$disposable.income,
                 Q = dat$catch,
                 
                 #Missing demand values
                 N_miss_P = sum(is.na(dat$price.per.kilo.wholesale)),
                 N_obs_P = sum(!is.na(dat$price.per.kilo.wholesale)),
                 ii_obs_P = which(!is.na(dat$price.per.kilo.wholesale)),
                 ii_mis_P = which(is.na(dat$price.per.kilo.wholesale)),
                 P_obs = dat$price.per.kilo.wholesale[!is.na(dat$price.per.kilo.wholesale)],
                 #Missing alternative price values
                 N_miss_Pf = sum(is.na(dat$price.per.kilo.farmed.wholesale)),
                 N_obs_Pf = sum(!is.na(dat$price.per.kilo.farmed.wholesale)),
                 ii_obs_Pf = which(!is.na(dat$price.per.kilo.farmed.wholesale)),
                 ii_mis_Pf = which(is.na(dat$price.per.kilo.farmed.wholesale)),
                 Pf_obs = dat$price.per.kilo.farmed.wholesale[!is.na(dat$price.per.kilo.farmed.wholesale)],
                 #Missing harvest values
                 N_miss_Q = sum(is.na(dat$harvest)),
                 N_obs_Q = sum(!is.na(dat$harvest)),
                 ii_obs_Q = which(!is.na(dat$harvest)),
                 ii_mis_Q = which(is.na(dat$harvest)),
                 Q_obs = dat$harvest[!is.na(dat$harvest)],
                 
                 spp_n = length(unique(dat$species)),
                 spp = dat$spp.n)
Demand <- Demand_mod$sample(datalist, parallel_chains = 4)


#Posterior predictive checks
norm_count_draws <- Demand$draws(variables = "preds")
#convert the drawn data into a matrix of 22 observations of 1000 draws
norm_draws_mat <- posterior::as_draws_matrix(norm_count_draws)
#now to plotting
bayesplot::ppc_dens_overlay(y = log(dat$price.per.kilo.wholesale[!is.na(dat$price.per.kilo.wholesale)]),
                            yrep = head((norm_draws_mat[,which(!is.na(dat$price.per.kilo.wholesale))]), 500))


#Correct Bayes R squared from Gelman: https://sites.stat.columbia.edu/gelman/research/unpublished/bayes_R2.pdf
bayes_R2 <- function(fit) {
  # Extract posterior draws of the linear predictor or predictions
  y_pred <- posterior::as_draws_matrix(fit$draws("preds"))
  # Calculate variance of the predictions for each posterior draw
  var_fit <- apply(y_pred, 1, var)
  # Extract posterior draws of the residual standard deviation (sigma) and square it to get residual variance
  sigma <- posterior::as_draws_matrix(fit$draws("sigma"))  # Replace "sigma" with the name of your residual standard deviation parameter
  var_res <- sigma^2
  # Compute Bayesian R-squared
  r2 <- var_fit / (var_fit + var_res)
  return(r2)}
## Compute Bayesian R2
rsq_bayes <- bayes_R2(Demand)
hist(rsq_bayes)
print(c(mean(rsq_bayes), median(rsq_bayes), sd(rsq_bayes)))

# Get unique species values
unique_species <- unique(dat$species)
# Loop over each species
plots <- map(unique_species, function(species_of_interest) {
  # Subset observed data, excluding NAs in price.per.kilo.wholesale
  observed_data_species <- log(dat$price.per.kilo.wholesale[dat$species == species_of_interest & !is.na(dat$price.per.kilo.wholesale)])
  # Subset posterior predictive draws, excluding NAs in price.per.kilo.wholesale
  species_indices <- which(dat$species == species_of_interest & !is.na(dat$price.per.kilo.wholesale))
  yrep_species <- norm_draws_mat[, species_indices]
  # Create plot
  plot <- bayesplot::ppc_dens_overlay(y = observed_data_species, yrep = yrep_species[1:500, ]) +
    ggtitle(paste("Posterior Predictive Check for", species_of_interest))
  # Save the plot
  file_name <- paste0("../Results/Demand/Diagnostics/Posterior_Predictive_Check_", species_of_interest, ".png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
  # Return the plot for further use (if needed)
  return(plot)
})


write_csv(Demand$summary(), "../Results/Demand/ModelParams.csv")


#Model predictions
Nsamp_pred <- cmdstan_model(stan_file = "DemandModels/DevoretzCatchSimple_Predictions.stan", pedantic=T)
reps <- 30
newdata <- expand_grid(species=c("Chinook","Coho","Chum","Pink","Sockeye"),
                       price.farmed= c(6.655788), #seq(from=4,to=11,length.out=reps), #
                       catch= seq(from=1,to=1000,length.out=reps), #c(250),#
                       disposable.income= c(226.243), #seq(from=123,to=360,length.out=reps), #
)
newdata$spp.n <- c(rep(5,reps),rep(4,reps),rep(3,reps),rep(2,reps),rep(1,reps))
newdata$price <- c(rep(12.399,reps),rep(28.6,reps),rep(5.483,reps),rep(4.123,reps),rep(22.909,reps))

data <-  list(newQ = newdata$catch,
              newP = newdata$price,
              newPf = newdata$price.farmed,
              newD = newdata$disposable.income,
              N = nrow(newdata),
              spp_n = length(unique(newdata$species)),
              spp = newdata$spp.n,
              N_miss_Pf = sum(is.na(newdata$price.farmed)),
              N_miss_P = sum(is.na(newdata$catch))
)
DemandPreds <- Nsamp_pred$generate_quantities(data, fitted_params = Demand)

png(file = "../Results/Demand/DemandCatch.png", width = 3600, height = 2400, res = 300)
DemandPreds |> 
  gather_rvars(real_mu[i]) |> 
  #gather_rvars(mu[i]) |> 
  bind_cols(newdata) |> 
  ggplot(aes(x = (catch), dist = .value)) + 
  #ggplot(aes(x = (disposable.income), dist = .value)) + 
  #ggplot(aes(x = (price.farmed), dist = .value)) + 
  stat_lineribbon() +
  #geom_point(data = dat, aes(x=price.per.kilo.farmed.wholesale, y=catch), colour = "salmon", alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(.~species, scales = "free") + 
  scale_fill_brewer(palette = "Greens", direction = -1) + 
  theme(legend.position = "none") + 
  #labs(x = "Disposable income (1000$)", y = "Predicted price ($/kg)")
  #labs(x = "Farmed salmon price ($/kg)", y = "Predicted price ($/kg)")
  labs(x = "Catch", y = "Predicted price ($/kg)")
dev.off()


Nsamp_pred <- cmdstan_model(stan_file = "DemandModels/DevoretzCatchSimple_Preds.stan", pedantic=T)
data <-  list(newQ = dat$catch,
              #Missing alternative price values
              N_miss_Pf = sum(is.na(dat$price.per.kilo.farmed.wholesale)),
              N_obs_Pf = sum(!is.na(dat$price.per.kilo.farmed.wholesale)),
              ii_obs_Pf = which(!is.na(dat$price.per.kilo.farmed.wholesale)),
              ii_mis_Pf = which(is.na(dat$price.per.kilo.farmed.wholesale)),
              Pf_obs = dat$price.per.kilo.farmed.wholesale[!is.na(dat$price.per.kilo.farmed.wholesale)],
              newD = dat$disposable.income,
              N = nrow(dat),
              spp_n = length(unique(dat$species)),
              spp = dat$spp.n,
              N_miss_P = sum(is.na(dat$price.per.kilo.wholesale))
)
DemandPreds <- Nsamp_pred$generate_quantities(data, fitted_params = Demand)
png(file = "../Results/Demand/DemandTimePoints.png", width = 3400, height = 2800, res = 300)
DemandPreds |> 
  gather_rvars(real_preds[i]) |> 
  #gather_rvars(mu[i]) |> 
  bind_cols(dat) |> 
  ggplot(aes(x = (year), dist = .value)) + 
  stat_lineribbon() +
  geom_point(data = dat, aes(x=year, y=price.per.kilo.wholesale), colour = "darkorange", alpha = 0.7, inherit.aes = FALSE, size = 2) +
  facet_wrap(.~species, scales = "free", ncol = 2) + 
  scale_fill_brewer(palette = "Greens", direction = -1, name = "Credible\ninterval") + 
  theme_minimal() +
  labs(x = "Year", y = "Predicted price ($/kg)") +
  theme(
    legend.position = c(0.8, 0.1),  # Place legend inside the plot at bottom right
    legend.justification = c(.5, .3),  # Anchor legend to the bottom right
    text = element_text(size = 25),  # Increase font size for all text
    strip.text = element_text(face = "bold", size = 20)  # Bold and increase facet title size
  )
dev.off()

