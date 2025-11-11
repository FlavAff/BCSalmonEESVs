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

source("Effort/EffortData.R")
dat$fleet <- dat$fleet/100
dat$licenses <- dat$licenses/100

#We are running two models here (effort & fleet)
AC_mod <- cmdstan_model(stan_file = "Effort/EffortMod.stan", pedantic=T) #ACLmod with licenses given a spread (.90/.65)
datalist <- list(N_tot=nrow(dat),
            N_fobs=sum(!is.na(dat$fleet)),
            N_lobs=sum(!is.na(dat$licenses)),
            N_eobs=sum(!is.na(dat$effort)),
            N_fmis=sum(is.na(dat$fleet)),
            N_lmis=sum(is.na(dat$licenses)),
            N_emis=sum(is.na(dat$effort)),
            
            ii_lobs = which(!is.na(dat$licenses)),
            ii_lmis = which(is.na(dat$licenses)),
            ii_fobs = which(!is.na(dat$fleet)),
            ii_fmis = which(is.na(dat$fleet)),
            ii_eobs = which(!is.na(dat$effort)),
            ii_emis = which(is.na(dat$effort)),
            
            fleet_obs=dat$fleet[!is.na(dat$fleet)],
            license_obs=dat$licenses[!is.na(dat$licenses)],
            effort_obs=dat$effort[!is.na(dat$effort)],
            
            N_gear=length(unique(dat$gear)),
            gear = dat$gear.n,
            N_area=length(unique(dat$area)),
            area = dat$area.n)
AC <- AC_mod$sample(datalist, parallel_chains = 4)


#Posterior predictive checks EFFORT
norm_count_draws <- AC$draws(variables = "effort_pred")
#convert the drawn data into a matrix of 22 observations of 1000 draws
norm_draws_mat <- posterior::as_draws_matrix(norm_count_draws)
#now to plotting
bayesplot::ppc_dens_overlay(y = (dat$effort[!is.na(dat$effort)]),
                            yrep = (head((norm_draws_mat[,which(!is.na(dat$effort))]), 500)))


#Correct Bayes R squared from Gelman: https://sites.stat.columbia.edu/gelman/research/unpublished/bayes_R2.pdf
bayes_R2 <- function(fit) {
  # Extract posterior draws of the linear predictor or predictions
  y_pred <- posterior::as_draws_matrix(fit$draws("effort_pred"))
  # Calculate variance of the predictions for each posterior draw
  var_fit <- apply(y_pred, 1, var)
  # Extract posterior draws of the residual standard deviation (sigma) and square it to get residual variance
  sigma <- posterior::as_draws_matrix(fit$draws("sigma_effort"))  # Replace "sigma" with the name of your residual standard deviation parameter
  var_res <- sigma^2
  # Compute Bayesian R-squared
  r2 <- var_fit / (var_fit + var_res)
  return(r2)}
## Compute Bayesian R2
rsq_bayes <- bayes_R2(AC)
hist(rsq_bayes)
print(c(mean(rsq_bayes), median(rsq_bayes), sd(rsq_bayes)))

# Get unique species values
unique_area <- unique(dat$area)
# Loop over each species
plots <- map(unique_area, function(area_of_interest) {
  # Subset observed data, excluding NAs in price.per.kilo.wholesale
  observed_data_area <- (dat$effort[dat$area == area_of_interest & !is.na(dat$effort)])
  # Subset posterior predictive draws, excluding NAs in price.per.kilo.wholesale
  area_indices <- which(dat$area == area_of_interest & !is.na(dat$effort))
  yrep_area <- norm_draws_mat[, area_indices]
  # Create plot
  plot <- bayesplot::ppc_dens_overlay(y = observed_data_area, yrep = yrep_area[1:500, ]) +
    ggtitle(paste("Posterior Predictive Check for", area_of_interest))
  # Save the plot
  file_name <- paste0("../Results/Effort/Diagnostics/PP_Check_Effort_", area_of_interest, ".png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
  # Return the plot for further use (if needed)
  return(plot)
})

# Get unique species values
unique_gear <- unique(dat$gear)
# Loop over each species
plots <- map(unique_gear, function(gear_of_interest) {
  # Subset observed data, excluding NAs in price.per.kilo.wholesale
  observed_data_gear <- (dat$effort[dat$gear == gear_of_interest & !is.na(dat$effort)])
  # Subset posterior predictive draws, excluding NAs in price.per.kilo.wholesale
  gear_indices <- which(dat$gear == gear_of_interest & !is.na(dat$effort))
  yrep_gear <- norm_draws_mat[, gear_indices]
  # Create plot
  plot <- bayesplot::ppc_dens_overlay(y = observed_data_gear, yrep = yrep_gear[1:500, ]) +
    ggtitle(paste("Posterior Predictive Check for", gear_of_interest))
  # Save the plot
  file_name <- paste0("../Results/Effort/Diagnostics/PP_Check_Effort_", gear_of_interest, ".png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
  # Return the plot for further use (if needed)
  return(plot)
})


####################################################
#Posterior predictive checks FLEET
norm_count_draws <- AC$draws(variables = "fleet_pred")
#convert the drawn data into a matrix of 22 observations of 1000 draws
norm_draws_mat <- posterior::as_draws_matrix(norm_count_draws)
#now to plotting
bayesplot::ppc_dens_overlay(y = (dat$fleet[!is.na(dat$fleet)]),
                            yrep = head((norm_draws_mat[,which(!is.na(dat$fleet))]), 500))

# Loop over each species
plots <- map(unique_area, function(area_of_interest) {
  # Subset observed data, excluding NAs in price.per.kilo.wholesale
  observed_data_area <- (dat$fleet[dat$area == area_of_interest & !is.na(dat$fleet)])
  # Subset posterior predictive draws, excluding NAs in price.per.kilo.wholesale
  area_indices <- which(dat$area == area_of_interest & !is.na(dat$fleet))
  yrep_area <- norm_draws_mat[, area_indices]
  # Create plot
  plot <- bayesplot::ppc_dens_overlay(y = observed_data_area, yrep = yrep_area[1:500, ]) +
    ggtitle(paste("Posterior Predictive Check for", area_of_interest))
  # Save the plot
  file_name <- paste0("../Results/Effort/Diagnostics/PP_Check_Fleet_", area_of_interest, ".png")
  ggsave(file_name, plot, width = 8, height = 6, dpi = 300)
  # Return the plot for further use (if needed)
  return(plot)
})

#Correct Bayes R squared from Gelman: https://sites.stat.columbia.edu/gelman/research/unpublished/bayes_R2.pdf
bayes_R2 <- function(fit) {
  # Extract posterior draws of the linear predictor or predictions
  y_pred <- posterior::as_draws_matrix(fit$draws("fleet_pred"))
  # Calculate variance of the predictions for each posterior draw
  var_fit <- apply(y_pred, 1, var)
  # Extract posterior draws of the residual standard deviation (sigma) and square it to get residual variance
  sigma <- posterior::as_draws_matrix(fit$draws("sigma_fleet"))  # Replace "sigma" with the name of your residual standard deviation parameter
  var_res <- sigma^2
  # Compute Bayesian R-squared
  r2 <- var_fit / (var_fit + var_res)
  return(r2)}
## Compute Bayesian R2
rsq_bayes <- bayes_R2(AC)
hist(rsq_bayes)
print(c(mean(rsq_bayes), median(rsq_bayes), sd(rsq_bayes)))

#Check parameters
AC |> 
  posterior::summarise_draws() |> 
  flextable::flextable()
shinystan::launch_shinystan(AC)

write_csv(AC$summary(), "../Results/Effort/ModelParams.csv")


##############################################################
#Model predictions
Nsamp_pred <- cmdstan_model(stan_file = "Effort/EffortMod_predictions.stan", pedantic=T)
reps <- 30
newdata <- expand_grid(area=c("A","B","C","D","E","F","G","H"),
                       fleet= seq(from=0,to=10,length.out=reps)
)
newdata <- newdata %>%
  mutate(
    gear = case_when(
      area %in% c("A", "B") ~ "Seine",
      area %in% c("C", "D", "E") ~ "Gillnet",
      area %in% c("F", "G", "H") ~ "Troll"
    ),
    area.n = match(area, c("A", "B", "C", "D", "E", "F", "G", "H")),
    gear.n = case_when(
      gear == "Seine" ~ 1,
      gear == "Gillnet" ~ 3,
      gear == "Troll" ~ 2
    )
  )
newdata$licenses <- rep(seq(from=0,to=9,length.out=reps),8)

data <-  list(new_fleet = newdata$fleet,
              new_license = newdata$licenses,
              N_tot = nrow(newdata),
              N_gear = length(unique(newdata$gear)),
              gear = newdata$gear.n,
              N_area = length(unique(newdata$area)),
              area = newdata$area.n,
              N_fmis = sum(is.na(newdata$fleet)),
              N_emis = sum(is.na(newdata$fleet)),
              N_lmis = sum(is.na(newdata$fleet)),
              license_mean = mean(log(dat$licenses[!is.na(dat$licenses)])),
              fleet_mean = mean(log(dat$fleet[!is.na(dat$fleet)]))
)
ACPreds <- Nsamp_pred$generate_quantities(data, fitted_params = AC)

png(file = "../Results/Effort/EffortFleet.png", width = 2400, height = 1800, res = 300)
ACPreds |> 
  gather_rvars(effort_pred[i]) |> 
  #gather_rvars(mu_effort[i]) |> 
  bind_cols(newdata) |> 
  ggplot(aes(x = (fleet - mean(fleet)), dist = .value)) + 
  stat_lineribbon() +
  #geom_point(data = dat, aes(x=fleet, y=effort), colour = "black", alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(area~gear, scales = "free") + 
  scale_fill_brewer(palette = "Oranges", direction = -1) + 
  theme(legend.position = "none") + 
  labs(x = "Average fleet size", y = "Predicted effort (boat days)")
dev.off()



png(file = "../Results/Effort/FleetLicense.png", width = 2400, height = 1800, res = 300)
ACPreds |> 
  gather_rvars(fleet_pred[i]) |> 
  #gather_rvars(mu_fleet[i]) |> 
  bind_cols(newdata) |> 
  ggplot(aes(x = (licenses - mean(licenses)), dist = .value)) + 
  stat_lineribbon() +
  #geom_point(data = dat, aes(x=licenses, y=fleet), colour = "black", alpha = 0.5, inherit.aes = FALSE) +
  facet_wrap(.~area, scales = "free") + 
  scale_fill_brewer(palette = "Oranges", direction = -1) + 
  theme(legend.position = "none") + 
  labs(x = "Average licenses", y = "Predicted fleet size")
dev.off()


#Plot time series with missing data
Nsamp_pred <- cmdstan_model(stan_file = "Effort/EffortMod_preds.stan", pedantic=T)
data <-  list(N_tot = nrow(dat),
              N_gear = length(unique(dat$gear)),
              gear = dat$gear.n,
              N_area = length(unique(dat$area)),
              area = dat$area.n,
              
              N_emis = sum(is.na(dat$effort)),
              ii_emis = which(is.na(dat$effort)),
              
              N_lmis = sum(is.na(dat$licenses)),
              N_lobs = sum(!is.na(dat$licenses)),
              ii_lobs = which(!is.na(dat$licenses)),
              ii_lmis = which(is.na(dat$licenses)),
              L_new_obs = (dat$licenses[!is.na(dat$licenses)]),
              
              N_fmis = sum(is.na(dat$fleet)),
              N_fobs = sum(!is.na(dat$fleet)),
              ii_fobs = which(!is.na(dat$fleet)),
              ii_fmis = which(is.na(dat$fleet)),
              F_new_obs = (dat$fleet[!is.na(dat$fleet)]),
              
              license_mean = mean((dat$licenses[!is.na(dat$licenses)])),
              fleet_mean = mean((dat$fleet[!is.na(dat$fleet)])))
ACPreds <- Nsamp_pred$generate_quantities(data, fitted_params = AC)


png(file = "../Results/Effort/EffortYearPoints.png", width = 3400, height = 3000, res = 300)
ACPreds |> 
  gather_rvars(effort_pred[i]) |> 
  #gather_rvars(mu_effort[i]) |> 
  bind_cols(dat) |> 
  ggplot(aes(x = (year), dist = .value)) + 
  stat_lineribbon() +
  geom_point(data = dat, aes(x=year, y=effort), colour = "blue", alpha = 0.7, inherit.aes = FALSE, size = 2) +
  facet_wrap(~ interaction(area, gear, sep = " - "), scales = "free", ncol = 2) + 
  scale_fill_brewer(palette = "Oranges", direction = -1, name = "Credible\ninterval") + 
  theme_minimal() + 
  labs(x = "Year", y = "Predicted effort (boat days)") +
  theme(
    legend.position = "none",  # Place legend inside the plot at bottom right
    legend.justification = c(.5, .3),  # Anchor legend to the bottom right
    text = element_text(size = 20),  # Increase font size for all text
    strip.text = element_text(face = "bold", size = 16)  # Bold and increase facet title size
  )
dev.off()

png(file = "../Results/Effort/FleetYearPoints.png", width = 3400, height = 1800, res = 300)
ACPreds |> 
  gather_rvars(fleet_pred[i]) |> 
  #gather_rvars(mu_fleet[i]) |> 
  bind_cols(dat) |> 
  ggplot(aes(x = (year), dist = .value)) + 
  stat_lineribbon() +
  geom_point(data = dat, aes(x=year, y=fleet), colour = "blue", alpha = 0.5, inherit.aes = FALSE, size = 2) +
  facet_wrap(~ interaction(area, gear, sep = " - "), scales = "free", ncol = 3) + 
  scale_fill_brewer(palette = "Oranges", direction = -1, name = "Credible\ninterval") + 
  theme_minimal() + xlim(2005,2023) +
  labs(x = "Year", y = "Predicted fleet size") +
  theme(
    legend.position = c(0.85, 0.05),  # Place legend inside the plot at bottom right
    legend.justification = c(.5, .3),  # Anchor legend to the bottom right
    text = element_text(size = 18),  # Increase font size for all text
    strip.text = element_text(face = "bold", size = 14)  # Bold and increase facet title size
  )
dev.off()

#Sum totals per year
datsum <- dat |> group_by(year) |> summarise(effort = sum(effort, na.rm = TRUE), .groups = 'drop')
png(file = "../Results/Effort/TotalEffortYear.png", width = 1200, height = 1200, res = 300)
ACPreds |> 
  gather_rvars(effort_pred[i]) |> 
  #gather_rvars(mu_effort[i]) |> 
  bind_cols(dat) |> 
  group_by(year) |>
  summarise(.value = rvar_sum(.value, na.rm = TRUE), .groups = 'drop') |>
  ggplot(aes(x = (year), dist = .value)) + 
  stat_lineribbon() +
  geom_point(data = datsum, aes(x=year, y=effort), colour = "blue", alpha = 0.5, inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Oranges", direction = -1) + 
  theme(legend.position = "none") + 
  labs(x = "Year", y = "Predicted effort (boat days)")
dev.off()

datsum <- dat |> group_by(year) |> summarise(fleet = sum(fleet, na.rm = TRUE), .groups = 'drop')
png(file = "../Results/Effort/TotalFleetYear.png", width = 1200, height = 1200, res = 300)
ACPreds |> 
  gather_rvars(fleet_pred[i]) |> 
  #gather_rvars(mu_effort[i]) |> 
  bind_cols(dat) |> 
  group_by(year) |>
  summarise(.value = rvar_sum(.value, na.rm = TRUE), .groups = 'drop') |>
  ggplot(aes(x = (year), dist = .value)) + 
  stat_lineribbon() +
  geom_point(data = datsum, aes(x=year, y=fleet), colour = "blue", alpha = 0.5, inherit.aes = FALSE) +
  scale_fill_brewer(palette = "Oranges", direction = -1) + 
  theme(legend.position = "none") + xlim(2005,2023) +
  labs(x = "Year", y = "Predicted fleet size")
dev.off()
