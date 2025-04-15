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

dat <- read_csv("../Data/Salmon/Ecological Supply/SupplyEESV.csv") |> 
  filter(!(species == "Chum" & region == "Fraser") & !(species == "Coho" & region == "VIMI"))


#THe modelling
Supply_mod <- cmdstan_model(stan_file = "Abundance/Supply_rickermodSSTHatch.stan", pedantic=T) #adds hatcheries (.91)
datalist <- list(N = nrow(dat),
                 
                 N_miss_T = sum(is.na(dat$parental_sst)),
                 N_obs_T = sum(!is.na(dat$parental_sst)),
                 ii_obs_T = which(!is.na(dat$parental_sst)),
                 ii_mis_T = which(is.na(dat$parental_sst)),
                 SST_obs = log(dat$parental_sst[!is.na(dat$parental_sst)]),
                 
                 N_miss_H = sum(is.na(dat$parental_release)),
                 N_obs_H = sum(!is.na(dat$parental_release)),
                 ii_obs_H = which(!is.na(dat$parental_release)),
                 ii_mis_H = which(is.na(dat$parental_release)),
                 H_obs = log(dat$parental_release[!is.na(dat$parental_release)]),
                 
                 N_miss_R = sum(is.na(dat$abundance)),
                 N_obs_R = sum(!is.na(dat$abundance)),
                 ii_obs_R = which(!is.na(dat$abundance)),
                 ii_mis_R = which(is.na(dat$abundance)),
                 log_R_obs = log(dat$abundance[!is.na(dat$abundance)]),
                 
                 N_miss_P = sum(is.na(dat$parental_abundance)),
                 N_obs_P = sum(!is.na(dat$parental_abundance)),
                 ii_obs_P = which(!is.na(dat$parental_abundance)),
                 ii_mis_P = which(is.na(dat$parental_abundance)),
                 log_P_obs = log(dat$parental_abundance[!is.na(dat$parental_abundance)]),
                 
                 spp_n = length(unique(dat$species)),
                 spp = dat$spp.n,
                 reg_n = length(unique(dat$region)),
                 reg = dat$reg.n)
Supply <- Supply_mod$sample(datalist, parallel_chains = 4,iter_warmup = 500,iter_sampling = 8000)

shinystan::launch_shinystan(Supply)

#Posterior predictive checks
norm_count_draws <- Supply$draws(variables = "preds")
#convert the drawn data into a matrix of 22 observations of 1000 draws
norm_draws_mat <- posterior::as_draws_matrix(norm_count_draws)
#now to plotting
bayesplot::ppc_dens_overlay(y = log(dat$abundance[!is.na(dat$abundance)]),
                            yrep = head((norm_draws_mat[,which(!is.na(dat$abundance))]), 500))


#Correct Bayes R squared from Gelman: https://sites.stat.columbia.edu/gelman/research/unpublished/bayes_R2.pdf
bayes_R2 <- function(fit) {
  # Extract posterior draws of the linear predictor or predictions
  y_pred <- posterior::as_draws_matrix(fit$draws("preds"))
  # Calculate variance of the predictions for each posterior draw
  var_fit <- apply(y_pred, 1, var)
  # Extract posterior draws of the residual standard deviation (sigma) and square it to get residual variance
  sigma <- posterior::as_draws_matrix(fit$draws("log_sigma"))  # Replace "sigma" with the name of your residual standard deviation parameter
  var_res <- sigma^2
  # Compute Bayesian R-squared
  r2 <- var_fit / (var_fit + var_res)
  return(r2)}
## Compute Bayesian R2
rsq_bayes <- bayes_R2(Supply)
hist(rsq_bayes)
print(c(mean(rsq_bayes), median(rsq_bayes), sd(rsq_bayes)))

# Get unique species values
unique_species <- unique(dat$species)
# Loop over each species
plots <- map(unique_species, function(species_of_interest) {
  # Subset observed data, excluding NAs in price.per.kilo.wholesale
  observed_data_species <- log(dat$abundance[dat$species == species_of_interest & !is.na(dat$abundance)])
  # Subset posterior predictive draws, excluding NAs in price.per.kilo.wholesale
  species_indices <- which(dat$species == species_of_interest & !is.na(dat$abundance))
  yrep_species <- norm_draws_mat[, species_indices]
  # Create plot
  plot <- bayesplot::ppc_dens_overlay(y = observed_data_species, yrep = yrep_species[1:500, ]) +
    ggtitle(paste("Posterior Predictive Check for", species_of_interest))
  # Save the plot
  file_name <- paste0("../Results/Supply/Diagnostics/MatrixParams/PP_Check_", species_of_interest, ".png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
  # Return the plot for further use (if needed)
  return(plot)
})

# Get unique species values
unique_region <- unique(dat$region)
# Loop over each species
plots <- map(unique_region, function(region_of_interest) {
  # Subset observed data, excluding NAs in price.per.kilo.wholesale
  observed_data_species <- log(dat$abundance[dat$region == region_of_interest & !is.na(dat$abundance)])
  # Subset posterior predictive draws, excluding NAs in price.per.kilo.wholesale
  species_indices <- which(dat$region == region_of_interest & !is.na(dat$abundance))
  yrep_species <- norm_draws_mat[, species_indices]
  # Create plot
  plot <- bayesplot::ppc_dens_overlay(y = observed_data_species, yrep = yrep_species[1:500, ]) +
    ggtitle(paste("Posterior Predictive Check for", region_of_interest))
  # Save the plot
  file_name <- paste0("../Results/Supply/Diagnostics/MatrixParams/PP_Check_", region_of_interest, ".png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
  # Return the plot for further use (if needed)
  return(plot)
})


#extract posterior predictions into a df
draws_df <- tidybayes::spread_rvars(Supply, preds[i], mu_log[i])|>
  median_qi(preds,mu_log) |> mutate(abundance = dat$abundance[i],
                                    parental_abundance = dat$parental_abundance[i],
                                    parental_release = dat$parental_release[i],
                                    parental_sst = dat$parental_sst[i],
                                    region = dat$region[i],
                                    species = dat$species[i])
# posterior predictive check
ggplot() + 
  geom_density(data = draws_df, aes((preds),fill = 'Posterior Predictive'), alpha = 0.5) + 
  geom_density(aes(log(dat$abundance), fill = 'Observed'), alpha = 0.5)
# easy plot
draws_df |> 
  ggplot(aes(x = log(parental_abundance), y = (preds))) + 
  geom_ribbon(aes(ymin=(preds.lower), ymax=(preds.upper)), fill = "salmon", alpha = 0.5) +
  geom_line(alpha = 1, col = "red") + 
  geom_point(aes(x = log(parental_abundance), y = log(abundance)),
             data = dat, 
             inherit.aes = FALSE, alpha = 0.7) + 
  facet_wrap(.~species, scales = "free") 


write_csv(Supply$summary(), "../Results/Supply/ModelParamsNew.csv")



## Predictions ##
Nsamp_pred <- cmdstan_model(stan_file = "Abundance/Supply_rickermodSSTHatchnoNA_predictions.stan", pedantic=T)

reps <- 50
newdata <- expand_grid(species=c("Pink","Chum","Coho","Chinook","Sockeye"),
                       #Remove region if not using Mat
                       #region = c("Yukon","HG","Nass","Skeena","CC","Fraser","Columbia","VIMI"),
                       parental_abundance = c(222186), #seq(from=5,to=815850,length.out=reps),  #
                       parental_sst = c(5.689), #seq(from=4.724,to=7.296,length.out=reps), #
                       parental_release = seq(from=163007,to=10763799,length.out=reps) #c(0) #
)
newdata$spp.n <- c(rep(5,reps),rep(4,reps),rep(3,reps),rep(2,reps),rep(1,reps))
#If using regions:
#newdata$spp.n <- c(rep(5,reps*8),rep(4,reps*8),rep(3,reps*8),rep(2,reps*8),rep(1,reps*8))
#newdata$reg.n <- rep(c(rep(1,reps),rep(2,reps),rep(3,reps),rep(4,reps),rep(5,reps),rep(6,reps),rep(7,reps),rep(8,reps)),5)

data <-  list(N = nrow(newdata),
              N_miss_R = 0,
              
              N_miss_T = sum(is.na(newdata$parental_sst)),
              SST_new = (newdata$parental_sst[!is.na(newdata$parental_sst)]),
              
              N_miss_H = sum(is.na(newdata$parental_release)),
              H_new = (newdata$parental_release[!is.na(newdata$parental_release)]),
              
              N_miss_P = sum(is.na(newdata$parental_abundance)),
              log_P_new = log(newdata$parental_abundance[!is.na(newdata$parental_abundance)]),
              
              #reg_n = length(unique(newdata$region)),
              #reg = newdata$reg.n,
              spp_n = length(unique(newdata$species)),
              spp = newdata$spp.n)
ModPreds <- Nsamp_pred$generate_quantities(data, fitted_params = Supply)

#Plotting SST/parental spawners/hatcheries instead requires changing the newdata list and the hashed out lines here
png(file = "../Results/Supply/SupplyReleaseLog.png", width = 3600, height = 2400, res = 300)
ModPreds |> 
  #gather_rvars(log_mu[i]) |> 
  gather_rvars(mu[i]) |> 
  #gather_rvars(preds[i]) |> 
  bind_cols(newdata) |> 
  #ggplot(aes(x = log(parental_abundance), dist = .value)) + 
  #ggplot(aes(x = (parental_sst), dist = .value)) + 
  ggplot(aes(x = (parental_release), dist = .value)) + 
  stat_lineribbon() +
  #geom_point(data = dat, aes(x=log(abundance), y=log(parental_abundance)), colour = "orange", alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(.~species, scales = "free") + 
  scale_fill_brewer(palette = "Purples", direction = -1) + 
  theme(legend.position = "none") + 
  #labs(x = "Parental abundance", y = "Predicted abundance")
  #labs(x = "Sea surface temperature", y = "Predicted abundance")
  labs(x = "Hatcheries", y = "Predicted abundance")
dev.off()




#Use actual data for time series
Nsamp_pred <- cmdstan_model(stan_file = "Abundance/Supply_rickermodSSTHatch_predictions.stan", pedantic=T)
data <-  list(N = nrow(dat),
              N_miss_R = sum(is.na(dat$abundance)),
              
              N_miss_T = sum(is.na(dat$parental_sst)),
              N_obs_T = sum(!is.na(dat$parental_sst)),
              ii_obs_T = which(!is.na(dat$parental_sst)),
              ii_mis_T = which(is.na(dat$parental_sst)),
              SST_new_obs = log(dat$parental_sst[!is.na(dat$parental_sst)]),
              
              N_miss_H = sum(is.na(dat$parental_release)),
              N_obs_H = sum(!is.na(dat$parental_release)),
              ii_obs_H = which(!is.na(dat$parental_release)),
              ii_mis_H = which(is.na(dat$parental_release)),
              H_new_obs = log(dat$parental_release[!is.na(dat$parental_release)]),
              
              N_miss_P = sum(is.na(dat$parental_abundance)),
              N_obs_P = sum(!is.na(dat$parental_abundance)),
              ii_obs_P = which(!is.na(dat$parental_abundance)),
              ii_mis_P = which(is.na(dat$parental_abundance)),
              log_P_new_obs = log(dat$parental_abundance[!is.na(dat$parental_abundance)]),
              
              spp_n = length(unique(dat$species)),
              spp = dat$spp.n)
ModPreds <- Nsamp_pred$generate_quantities(data, fitted_params = Supply)
png(file = "../Results/Supply/AbundanceTime.png", width = 4200, height = 3600, res = 300)
ModPreds |> 
  #gather_rvars(mu[i]) |> 
  gather_rvars(preds[i]) |> 
  bind_cols(dat) |> 
  ggplot(aes(x = year, dist = .value)) + 
  stat_lineribbon() +
  #geom_point(data = dat, aes(x=year, y=abundance), colour = "orange", alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(region~species, scales = "free") + 
  scale_fill_brewer(palette = "Purples", direction = -1) + 
  theme(legend.position = "none") + 
  labs(x = "Year", y = "Predicted abundance")
dev.off()


#sum the posteriors to get graphs per species
datsum <- dat |> group_by(year, species) |> summarise(abundance = sum(abundance, na.rm = TRUE), .groups = 'drop')
png(file = "../Results/Supply/AbundanceTimeSppPoints.png", width = 3600, height = 2400, res = 300)
ModPreds |> 
  gather_rvars(preds[i]) |> 
  #expgather_rvars(mu[i]) |> 
  bind_cols(dat) |>
  group_by(year, species) |>
  summarise(.value = rvar_sum(.value, na.rm = TRUE), .groups = 'drop') |>
  ggplot(aes(x = year, dist = (.value))) + 
  stat_lineribbon() +
  geom_point(data = datsum, aes(x=year, y=abundance), colour = "orange", alpha = 0.5, inherit.aes = FALSE, size = 2) +
  facet_wrap(.~species, scales = "free") + 
  scale_fill_brewer(palette = "Purples", direction = -1, name = "Credible\ninterval") + 
  theme_minimal() + 
  labs(x = "Year", y = "Predicted abundance") +
  theme(
    legend.position = c(0.85, 0.15),  # Place legend inside the plot at bottom right
    legend.justification = c(.5, .3),  # Anchor legend to the bottom right
    text = element_text(size = 18),  # Increase font size for all text
    strip.text = element_text(face = "bold", size = 14)  # Bold and increase facet title size
  )
dev.off()
