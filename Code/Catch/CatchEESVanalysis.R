setwd("~/Documents/Code/")
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

dat <- read_csv("../Data/Salmon/AreaEESVs2.csv")
Dmatrix <- as.matrix(read_csv("../Data/BC_management_areas/AreaMatrix.csv"))

#convert 0s to 1s so when logged in stan it doesn't explode
dat <- dat |> mutate(catch = ifelse(catch == 0, 1, catch)) |> 
  mutate(effort = ifelse(effort == 0, 1, effort)) |> 
  mutate(abundance = ifelse(abundance == 0, NA, abundance)) |> 
  #use this one if using dat as 0s are NAs
  #mutate(abundance = ifelse(abundance == 0, NA, abundance)) |> 
  mutate(log_catch = log(catch), log_effort = log(effort), log_abundance = log(abundance))
dat <- dat |> filter(log_catch != "NA")

#Start modelling
catch_mod <- cmdstan_model(stan_file = "Catch/SchaeferModMatExtremes.stan", pedantic=T) #seems to be even bestester .73/.46
datalist <- list(N_tot = nrow(dat),
                 log_effort = dat$log_effort,
                 log_catch = dat$log_catch,
                 Catch = dat$catch,
                 effort = dat$effort,

                 #Missing abundance values
                 N_a_mis = sum(is.na(dat$log_abundance)),
                 N_a_obs = sum(!is.na(dat$log_abundance)),
                 ii_a_obs = which(!is.na(dat$log_abundance)),
                 ii_a_mis = which(is.na(dat$log_abundance)),
                 log_abundance_obs = dat$log_abundance[!is.na(dat$log_abundance)],
                 #Missing catch values
                 # N_c_mis = sum(is.na(dat$log_catch)),
                 # N_c_obs = sum(!is.na(dat$log_catch)),
                 # ii_c_obs = which(!is.na(dat$log_catch)),
                 # ii_c_mis = which(is.na(dat$log_catch)),
                 # log_catch_obs = dat$log_catch[!is.na(dat$log_catch)],
                 
                 N_spp = length(unique(dat$species)),
                 spp = dat$species.n,
                 N_area = length(unique(dat$area)),
                 area = dat$area.n,
                 N_gear = length(unique(dat$gear)),
                 gear = dat$gear.n,
                 D = Dmatrix
                 )
Catch <- catch_mod$sample(datalist, parallel_chains = 4)

#Posterior predictive checks
norm_count_draws <- Catch$draws(variables = "preds")
#convert the drawn data into a matrix of 22 observations of 1000 draws
norm_draws_mat <- posterior::as_draws_matrix(norm_count_draws)
#now to plotting
bayesplot::ppc_dens_overlay(y = log(dat$catch),
                            yrep = head((norm_draws_mat), 500)) #+ xlim(-100,100)


# Get unique species values
unique_species <- unique(dat$species)

# Loop over each species
plots <- map(unique_species, function(species_of_interest) {
  # Subset observed data
  observed_data_species <- log(dat$catch[dat$species == species_of_interest])
  
  # Subset posterior predictive draws
  species_indices <- which(dat$species == species_of_interest)
  yrep_species <- norm_draws_mat[, species_indices]
  
  # Create plot
  plot <- bayesplot::ppc_dens_overlay(y = observed_data_species, yrep = yrep_species[1:500, ]) +
    ggtitle(paste("Posterior Predictive Check for", species_of_interest)) + xlim(-100,100)
  
  # Save the plot
  file_name <- paste0("../Results/Catch/Diagnostics/PP_Check_", species_of_interest, ".png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
  
  # Return the plot for further use (if needed)
  return(plot)
})

# Get unique area values
unique_area <- unique(dat$area)
# Loop over each area
plots <- map(unique_area, function(area_of_interest) {
  # Subset observed data
  observed_data_species <- log(dat$catch[dat$area == area_of_interest])
  
  # Subset posterior predictive draws
  area_indices <- which(dat$area == area_of_interest)
  yrep_area <- norm_draws_mat[, area_indices]
  
  # Create plot
  plot <- bayesplot::ppc_dens_overlay(y = observed_data_species, yrep = yrep_area[1:500, ]) +
    ggtitle(paste("Posterior Predictive Check for Area", area_of_interest)) + xlim(-100,100)
  
  # Save the plot
  file_name <- paste0("../Results/Catch/Diagnostics/PP_Check_Area_", area_of_interest, ".png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
  
  # Return the plot for further use (if needed)
  return(plot)
})


#Bayes R squared from Gelman: https://sites.stat.columbia.edu/gelman/research/unpublished/bayes_R2.pdf
bayes_R2 <- function(fit) {
  # Extract posterior draws of the linear predictor or predictions
  y_pred <- posterior::as_draws_matrix(fit$draws("preds"))
  # Calculate variance of the predictions for each posterior draw
  var_fit <- apply(y_pred, 1, var)
  # Extract posterior draws of the residual standard deviation (sigma) and square it to get residual variance
  sigma <- posterior::as_draws_matrix(fit$draws("sigma_obs"))  # Replace "sigma" with the name of your residual standard deviation parameter
  var_res <- sigma^2
  # Compute Bayesian R-squared
  r2 <- var_fit / (var_fit + var_res)
  return(r2)}
## Compute Bayesian R2
rsq_bayes <- bayes_R2(Catch)
hist(rsq_bayes)
print(c(mean(rsq_bayes), median(rsq_bayes), sd(rsq_bayes)))

Catch |> 
  posterior::summarise_draws() |> 
  flextable::flextable()

write_csv(Catch$summary(), "../Results/Catch/ModelParams.csv")


#Model predictions
Nsamp_pred <- cmdstan_model(stan_file = "Catch/SchaeferPredictionsExtremes2.stan", pedantic=T) 

reps <- 50
newdata <- expand_grid(species=c("Pink","Chum","Coho","Chinook","Sockeye"),
                       area = c("A","B","C","D","E","F","G","H"),
                       effort = c(8.2553), #seq(from=0.05,to=60,length.out=reps),  #
                       abundance = seq(from=1000,to=3000000,length.out=reps) #c(200000), #
)
newdata$spp.n <- c(rep(5,reps*8),rep(4,reps*8),rep(3,reps*8),rep(2,reps*8),rep(1,reps*8))
newdata$area.n <- rep(c(rep(1,reps),rep(2,reps),rep(3,reps),rep(4,reps),rep(5,reps),rep(6,reps),rep(7,reps),rep(8,reps)),5)


data <-  list(log_abundanceF = log(newdata$abundance),
              log_effortF = log(newdata$effort),
              N_tot = nrow(newdata),
              N_spp = length(unique(newdata$species)),
              spp = newdata$spp.n,
              N_area = length(unique(newdata$area)),
              area = newdata$area.n,
              N_a_mis=sum(is.na(newdata$abundance)))
ModPreds <- Nsamp_pred$generate_quantities(data, fitted_params = Catch)

#Plotting effort/abundance instead requires changing the newdata list and the hashed out lines here
#We end up using mu and not preds for the student t allows for such wide bounds that exponentiating the CIs goes to infinity.
png(file = "../Results/Catch/CatchAbundance.png", width = 7200, height = 3600, res = 300)
ModPreds |> 
  gather_rvars(mu_real[i]) |> 
  #gather_rvars(mu[i]) |> 
  bind_cols(newdata) |> 
  ggplot(aes(x = (abundance), dist = .value)) + 
  #ggplot(aes(x = (effort), dist = .value)) + 
  stat_lineribbon() +
  #geom_point(data = dat, aes(x=abundance, y=catch), colour = "salmon", alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(species~area, scales = "free") + 
  scale_fill_brewer(palette = "Reds", direction = -1) + 
  theme(legend.position = "none") + #xlim(1,60) +
  labs(x = "Abundance", y = "Predicted catch")
  #labs(x = "Effort (boat days)", y = "Predicted catch")
dev.off()

#Plot sums per species
png(file = "../Results/Catch/CatchAbundanceSpp.png", width = 3600, height = 2400, res = 300)
ModPreds |> 
  gather_rvars(mu_real[i]) |> 
  #gather_rvars(preds_real[i]) |> 
  bind_cols(newdata) |> 
  group_by(species,abundance) |>
  #group_by(species,effort) |>
  summarise(.value = rvar_sum(.value, na.rm = TRUE), .groups = 'drop') |>
  ggplot(aes(x = (abundance), dist = .value)) + 
  #ggplot(aes(x = (effort), dist = .value)) + 
  stat_lineribbon() +
  #geom_point(data = dat, aes(x=abundance, y=catch), colour = "green", alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(.~species, scales = "free") + 
  scale_fill_brewer(palette = "Reds", direction = -1) + 
  theme(legend.position = "none") + #xlim(1,60) +
  labs(x = "Abundance", y = "Predicted catch")
  #labs(x = "Effort (boat days)", y = "Predicted catch")
dev.off()

#############################
#Use actual data for time series 
#Nsamp_pred <- cmdstan_model(stan_file = "Catch/CatchPredictions6.stan", pedantic=T) #non hierarchical hurdle model with matrix beta
Nsamp_pred <- cmdstan_model(stan_file = "Catch/SchaeferPredictionsExtremes3.stan", pedantic=T) 
data <-  list(N_a_mis = sum(is.na(dat$log_abundance)),
              N_a_obs = sum(!is.na(dat$log_abundance)),
              ii_a_obs = which(!is.na(dat$log_abundance)),
              ii_a_mis = which(is.na(dat$log_abundance)),
              log_abundance_obs = dat$log_abundance[!is.na(dat$log_abundance)],
              log_effortF = dat$log_effort,
              N_tot = nrow(dat),
              N_spp = length(unique(dat$species)),
              spp = dat$species.n,
              N_area = length(unique(dat$area)),
              area = dat$area.n)
ModPreds <- Nsamp_pred$generate_quantities(data, fitted_params = Catch)
#png(file = "../Results/Catch/CatchTime.png", width = 7200, height = 3600, res = 300)
png(file = "../Results/Catch/SchaeferTime.png", width = 7200, height = 3600, res = 300)
ModPreds |> 
  gather_rvars(mu_real[i]) |> 
  bind_cols(dat) |> 
  ggplot(aes(x = year, dist = .value)) + 
  stat_lineribbon() +
  #geom_point(data = dat, aes(x=year, y=catch), colour = "darkgreen", alpha = 0.7, inherit.aes = FALSE, size = 2) +
  facet_wrap(species~area, scales = "free") + 
  scale_fill_brewer(palette = "Reds", direction = -1, name = "Credible\ninterval") + 
  theme_minimal() + 
  labs(x = "Year", y = "Predicted catch") + 
  theme(
    legend.position = c(0.85, 0.05),  # Place legend inside the plot at bottom right
    legend.justification = c(.5, .3),  # Anchor legend to the bottom right
    text = element_text(size = 18),  # Increase font size for all text
    strip.text = element_text(face = "bold", size = 14)  # Bold and increase facet title size
  )
dev.off()


#sum the posteriors to get graphs per species
datsum <- dat |> group_by(year, species) |> summarise(catch = sum(catch, na.rm = TRUE), .groups = 'drop')
#png(file = "../Results/Catch/CatchTimeSpp.png", width = 3600, height = 2400, res = 300)
png(file = "../Results/Catch/SchaeferTimeSppPoints.png", width = 3600, height = 1800, res = 300)
ModPreds |>
  gather_rvars(mu_real[i]) |> 
  bind_cols(dat) |>
  group_by(year, species) |>
  summarise(.value = rvar_sum(.value, na.rm = TRUE), .groups = 'drop') |>
  ggplot(aes(x = year, dist = (.value))) + 
  stat_lineribbon() +
  geom_point(data = datsum, aes(x=year, y=catch), colour = "darkgreen", alpha = 0.7, inherit.aes = FALSE, size = 2) +
  facet_wrap(.~species, scales = "free") + 
  scale_fill_brewer(palette = "Reds", direction = -1, name = "Credible\ninterval") + 
  theme_minimal() +
  labs(x = "Year", y = "Predicted catch") +
  theme(
    legend.position = c(0.85, 0.15),  # Place legend inside the plot at bottom right
    legend.justification = c(.5, .3),  # Anchor legend to the bottom right
    text = element_text(size = 18),  # Increase font size for all text
    strip.text = element_text(face = "bold", size = 14)  # Bold and increase facet title size
  )
dev.off()
