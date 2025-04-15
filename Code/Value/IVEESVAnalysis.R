setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(cmdstanr)
library(ggplot2)
library(tidyverse)
library(tidybayes)
library(shinystan)
library(rstan)
library(posterior)
library(loo)

source("Value/IVData.R")

#THe modelling
IV_mod <- cmdstan_model(stan_file = "Value/IVmod_basic3.stan", pedantic=T) #adds spp effect on both (.96/.93)
datalist <- list(N = nrow(dat),
                 catch_ = dat$catch,
                 harvest = dat$harvest,
                 value = dat$landed.value,
                 
                 N_spp = length(unique(dat$species)),
                 spp = dat$spp.n)
IV <- IV_mod$sample(datalist, parallel_chains = 4)

#Posterior predictive checks HARVEST
norm_count_draws <- IV$draws(variables = "harvest_pred")
#convert the drawn data into a matrix of 22 observations of 1000 draws
norm_draws_mat <- posterior::as_draws_matrix(norm_count_draws)
#now to plotting
bayesplot::ppc_dens_overlay(y = log(dat$harvest),
                            yrep = head((norm_draws_mat), 500))


#Correct Bayes R squared from Gelman: https://sites.stat.columbia.edu/gelman/research/unpublished/bayes_R2.pdf
bayes_R2 <- function(fit) {
  # Extract posterior draws of the linear predictor or predictions
  y_pred <- posterior::as_draws_matrix(fit$draws("harvest_pred"))
  # Calculate variance of the predictions for each posterior draw
  var_fit <- apply(y_pred, 1, var)
  # Extract posterior draws of the residual standard deviation (sigma) and square it to get residual variance
  sigma <- posterior::as_draws_matrix(fit$draws("sigma_h"))  # Replace "sigma" with the name of your residual standard deviation parameter
  var_res <- sigma^2
  # Compute Bayesian R-squared
  r2 <- var_fit / (var_fit + var_res)
  return(r2)}
## Compute Bayesian R2
rsq_bayes <- bayes_R2(IV)
hist(rsq_bayes)
print(c(mean(rsq_bayes), median(rsq_bayes), sd(rsq_bayes)))


# Get unique species values
unique_species <- unique(dat$species)
# Loop over each species
plots <- map(unique_species, function(species_of_interest) {
  # Subset observed data
  observed_data_species <- log(dat$harvest[dat$species == species_of_interest])
  # Subset posterior predictive draws
  species_indices <- which(dat$species == species_of_interest)
  yrep_species <- norm_draws_mat[, species_indices]
  # Create plot
  plot <- bayesplot::ppc_dens_overlay(y = observed_data_species, yrep = yrep_species[1:500, ]) +
    ggtitle(paste("Posterior Predictive Check for", species_of_interest))
  # Save the plot
  file_name <- paste0("../Results/Value/Diagnostics/PP_Check_Harvest_", species_of_interest, ".png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
  # Return the plot for further use (if needed)
  return(plot)
})


#Posterior predictive checks VALUE
norm_count_draws <- IV$draws(variables = "value_pred")
#convert the drawn data into a matrix of 22 observations of 1000 draws
norm_draws_mat <- posterior::as_draws_matrix(norm_count_draws)
#now to plotting
bayesplot::ppc_dens_overlay(y = log(dat$landed.value),
                            yrep = head((norm_draws_mat), 500))


#Correct Bayes R squared from Gelman: https://sites.stat.columbia.edu/gelman/research/unpublished/bayes_R2.pdf
bayes_R2 <- function(fit) {
  # Extract posterior draws of the linear predictor or predictions
  y_pred <- posterior::as_draws_matrix(fit$draws("value_pred"))
  # Calculate variance of the predictions for each posterior draw
  var_fit <- apply(y_pred, 1, var)
  # Extract posterior draws of the residual standard deviation (sigma) and square it to get residual variance
  sigma <- posterior::as_draws_matrix(fit$draws("sigma_v"))  # Replace "sigma" with the name of your residual standard deviation parameter
  var_res <- sigma^2
  # Compute Bayesian R-squared
  r2 <- var_fit / (var_fit + var_res)
  return(r2)}
## Compute Bayesian R2
rsq_bayes <- bayes_R2(IV)
hist(rsq_bayes)
print(c(mean(rsq_bayes), median(rsq_bayes), sd(rsq_bayes)))


# Get unique species values
unique_species <- unique(dat$species)
# Loop over each species
plots <- map(unique_species, function(species_of_interest) {
  # Subset observed data
  observed_data_species <- log(dat$landed.value[dat$species == species_of_interest])
  # Subset posterior predictive draws
  species_indices <- which(dat$species == species_of_interest)
  yrep_species <- norm_draws_mat[, species_indices]
  # Create plot
  plot <- bayesplot::ppc_dens_overlay(y = observed_data_species, yrep = yrep_species[1:500, ]) +
    ggtitle(paste("Posterior Predictive Check for", species_of_interest))
  # Save the plot
  file_name <- paste0("../Results/Value/Diagnostics/PP_Check_Value_", species_of_interest, ".png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
  # Return the plot for further use (if needed)
  return(plot)
})


IV |> 
  posterior::summarise_draws() |> 
  flextable::flextable()

write_csv(IV$summary(), "../Results/Value/ModelParams.csv")

################# PREDICTIONS
#Model predictions
Nsamp_pred <- cmdstan_model(stan_file = "Value/IVmod3_predictions.stan", pedantic=T) #non hierarchical hurdle model with matrix beta
reps <- 50
newdata <- expand_grid(species=c("Pink","Chum","Coho","Chinook","Sockeye"),
                       catch= seq(from=0.01,to=30,length.out=reps))
newdata$spp.n <- c(rep(5,reps),rep(4,reps),rep(3,reps),rep(2,reps),rep(1,reps))
newdata$harvest <- rep(seq(from=0.1,to=28,length.out=reps),5)
data <-  list(new_catch = (newdata$catch),
              new_harvest = (newdata$harvest),
              N = nrow(newdata),
              N_spp = length(unique(newdata$species)),
              spp = newdata$spp.n)
ModPreds <- Nsamp_pred$generate_quantities(data, fitted_params = IV)

png(file = "../Results/Value/HarvestCatch.png", width = 3600, height = 2400, res = 300)
ModPreds |> 
  gather_rvars(mu_harvest[i]) |> 
  #gather_rvars(harvest_pred[i]) |> 
  bind_cols(newdata) |> 
  ggplot(aes(x = (catch), dist = exp(.value))) + 
  stat_lineribbon() +
  #geom_point(data = dat, aes(x=catch, y=harvest), colour = "salmon", alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(.~species, scales = "free") + 
  scale_fill_brewer(palette = "Blues", direction = -1) + 
  theme(legend.position = "none") + 
  labs(x = "Catch", y = "Predicted harvest (tonnes)")
dev.off()

png(file = "../Results/Value/ValueHarvest.png", width = 3600, height = 2400, res = 300)
ModPreds |> 
  gather_rvars(mu_value[i]) |> 
  #gather_rvars(value_pred[i]) |> 
  bind_cols(newdata) |> 
  ggplot(aes(x = (harvest), dist = exp(.value))) + 
  stat_lineribbon() +
  #geom_point(data = dat, aes(x=harvest, y=landed.value), colour = "salmon", alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(.~species, scales = "free") + 
  scale_fill_brewer(palette = "Blues", direction = -1) + 
  theme(legend.position = "none") + 
  labs(x = "Harvest (tonnes)", y = "Predicted value (M$)")
dev.off()


#With real data for TS
newdata2 <- dat
data <-  list(new_catch = (newdata2$catch),
              new_harvest = (newdata2$harvest),
              N = nrow(newdata2),
              N_spp = length(unique(newdata2$species)),
              spp = newdata2$spp.n)
ModPreds <- Nsamp_pred$generate_quantities(data, fitted_params = IV)

png(file = "../Results/Value/HarvestYearPoints.png", width = 3000, height = 1800, res = 300)
ModPreds |> 
  gather_rvars(harvest_pred[i]) |> 
  #gather_rvars(mu_harvest[i]) |> 
  bind_cols(newdata2) |> 
  ggplot(aes(x = (year), dist = exp(.value))) + 
  stat_lineribbon() +
  geom_point(data = dat, aes(x=year, y=harvest), colour = "red", alpha = 0.5, inherit.aes = FALSE, size = 2) +
  facet_wrap(.~species, scales = "free") + 
  scale_fill_brewer(palette = "Blues", direction = -1, name = "Credible\ninterval") + 
  theme_minimal() +
  labs(x = "Year", y = "Predicted harvest (tonnes)") +
  theme(
    legend.position = c(0.85, 0.15),  # Place legend inside the plot at bottom right
    legend.justification = c(.5, .3),  # Anchor legend to the bottom right
    text = element_text(size = 18),  # Increase font size for all text
    strip.text = element_text(face = "bold", size = 14)  # Bold and increase facet title size
  )
dev.off()

png(file = "../Results/Value/ValueYearPoints.png", width = 3000, height = 1800, res = 300)
ModPreds |> 
  gather_rvars(value_pred[i]) |> 
  #gather_rvars(mu_value[i]) |> 
  bind_cols(newdata2) |> 
  ggplot(aes(x = (year), dist = exp(.value))) + 
  stat_lineribbon() +
  geom_point(data = dat, aes(x=year, y=landed.value), colour = "red", alpha = 0.5, inherit.aes = FALSE, size = 2) +
  facet_wrap(.~species, scales = "free") + 
  scale_fill_brewer(palette = "Blues", direction = -1, name = "Credible\ninterval") + 
  theme_minimal() + 
  labs(x = "Year", y = "Predicted value (M$)") +
  theme(
    legend.position = c(0.85, 0.15),  # Place legend inside the plot at bottom right
    legend.justification = c(.5, .3),  # Anchor legend to the bottom right
    text = element_text(size = 18),  # Increase font size for all text
    strip.text = element_text(face = "bold", size = 14)  # Bold and increase facet title size
  )
dev.off()
