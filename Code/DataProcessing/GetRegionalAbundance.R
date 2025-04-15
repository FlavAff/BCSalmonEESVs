setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)
library(zoo)
library(boot)

salmon <- read_csv("../Data/Salmon/Ecological Supply/COSEWIC-compilation-main/Output/CU_Spawner_Abund_20241115.csv")

Chum <- filter(salmon, Species =="Chum")
Coho <- filter(salmon, Species =="Coho")
Chinook <- filter(salmon, Species =="Chinook")
Pink.E <- filter(salmon, Species =="Pink (Even)")
Pink.O <- filter(salmon, Species =="Pink (Odd)")
Sockeye <- filter(salmon, Species =="Sockeye")

plotregionaloutliers <- function(dat) {
  #fish <- dat
  dat$YearF <- as.factor(dat$Year)
  name <- dat$Species[1]
  p <- dat |> ggplot(aes(x=YearF, y = Spawner.Abundance)) + geom_boxplot(outlier.colour="red") + facet_wrap(~Region)
  png(paste0("../Results/RegionalAbundances/",name,"Abundances.png"),width=30, height=20, units="cm", res=1200)
  print(p)
  dev.off()
}

plotregionaloutliers(Chum)
plotregionaloutliers(Coho)
plotregionaloutliers(Sockeye)
plotregionaloutliers(Chinook)
plotregionaloutliers(Pink.E)
plotregionaloutliers(Pink.O)

getTS <- function(fish) {
  df <- NULL
  for (i in unique(fish$Region)) {
    region <- filter(fish, Region==i)
    TS <- tapply(region$Spawner.Abundance, region$Year, mean, na.rm=T)
    dat <- tibble(Abundance=TS,Region=i,Year=seq(min(region$Year),max(region$Year),1), species = rep(fish$Species[1]))
    df <- rbind(dat,df)
  }
  return(df)
}

ChumTS <- getTS(Chum)
CohoTS <- getTS(Coho)
ChinookTS <- getTS(Chinook)
PinkE.TS <- getTS(Pink.E)
PinkO.TS <- getTS(Pink.O)
SockeyeTS <- getTS(Sockeye)

regionalaverages <- rbind(ChumTS,CohoTS,ChinookTS,PinkE.TS,PinkO.TS,SockeyeTS)
write.csv(regionalaverages, "../Data/Salmon/Ecological Supply/AverageAbundanceRegion.csv")

# write.csv(ChumTS, "../Data/Salmon/Ecological Supply/ChumAverageRegion.csv")
# write.csv(CohoTS, "../Data/Salmon/Ecological Supply/CohoAverageRegion.csv")
# write.csv(ChinookTS, "../Data/Salmon/Ecological Supply/ChinookAverageRegion.csv")
# write.csv(SockeyeTS, "../Data/Salmon/Ecological Supply/SockeyeAverageRegion.csv")
# write.csv(PinkE.TS, "../Data/Salmon/Ecological Supply/PinkEAverageRegion.csv")
# write.csv(PinkO.TS, "../Data/Salmon/Ecological Supply/PinkOAverageRegion.csv")


gettotalTS <- function(fish) {
  df <- NULL
  TS <- tapply(fish$Spawner.Abundance, fish$Year, mean, na.rm=T)
  dat <- data_frame(Abundance=TS,Year=seq(min(fish$Year),max(fish$Year),1), species = rep(fish$Species[1]))
  df <- rbind(dat,df)
  return(df)
}

ChumTS <- gettotalTS(Chum)
CohoTS <- gettotalTS(Coho)
ChinookTS <- gettotalTS(Chinook)
PinkE.TS <- gettotalTS(Pink.E)
PinkO.TS <- gettotalTS(Pink.O)
SockeyeTS <- gettotalTS(Sockeye)

totalaverages <- rbind(ChumTS,CohoTS,ChinookTS,PinkE.TS,PinkO.TS,SockeyeTS)
write.csv(totalaverages, "../Data/Salmon/Ecological Supply/YearlyAverageAbundance.csv")



######### ISSUES OF MISSING DATA IN CUIDs to get a yearly sum - INTERPOLATION options ##############

#Simple attempts with reduced dataset
region <- filter(Chum, CUID == 211)   
#take out the tails to avoid GAM estimating empty data there
region <- filter(region, Year < 2018)
region <- filter(region, Year > 1953)

# Perform linear interpolation on Spawner.Abundance
regionapprox <- region %>%
  arrange(Year) %>%
  mutate(Spawner.Abundance = na.approx(Spawner.Abundance, x = Year, na.rm = FALSE))

# Perform GAM interpolation on Spawner.Abundance
model <- gam(Spawner.Abundance ~ s(Year), data = region, method = "REML")
predictions <- predict(model, newdata = data.frame(Year = seq(min(region$Year), max(region$Year), 1)), se.fit = TRUE)
region2$Spawner.Abundance <- ifelse(is.na(region$Spawner.Abundance), predictions$fit, region$Spawner.Abundance)
region2$Lower <- predictions$fit - 1.96 * predictions$se.fit
region2$Upper <- predictions$fit + 1.96 * predictions$se.fit

#Check the differences between the original dataset, the linear interpolation and GAM
plot(regionapprox$Spawner.Abundance)
plot(region$Spawner.Abundance)
plot(region2$Spawner.Abundance)

#A Bayesian approach using brms
library(tidyr)
library(brms)

#Back to using all Chum data
Chum <- Chum %>%
  mutate(CUID = factor(CUID), Year = as.numeric(Year))

# Fit the model using brms
fit <- brm(
  Spawner.Abundance ~ s(Year) + (1 | CUID),
  data = Chum,
  family = lognormal(),
  chains = 4, cores = 4, iter = 2000, control = list(adapt_delta = 0.95)
)

# Add predictions to the data
Chum$Spawner.Abundance <- ifelse(
  is.na(Chum$Spawner.Abundance),
  posterior_predict(fit, newdata = Chum, re_formula = NA) %>% colMeans(),
  Chum$Spawner.Abundance
)

# Aggregate data
# Extract posterior draws for total abundance
posterior_totals <- Chum %>%
  group_by(Year) %>%
  summarise(Total.Abundance = sum(Spawner.Abundance, na.rm = TRUE))

# Add credible intervals for each year
posterior_samples <- posterior_predict(fit, newdata = Chum)
total_by_year <- apply(posterior_samples, 1, function(row) {
  tapply(row, Chum$Year, sum, na.rm = TRUE)
})

# Calculate credible intervals
# Transpose total_by_year matrix for easier handling
total_by_year_df <- as.data.frame(t(total_by_year))
# Convert to a long format for summarization
total_by_year_long <- total_by_year_df %>%
  pivot_longer(cols = everything(), names_to = "Posterior_Sample", values_to = "Posterior_Abundance")
# Summarize posterior distributions for each year
total_by_year_summary <- total_by_year_long %>%
  group_by(Posterior_Sample) %>%
  summarise(
    Mean = mean(Posterior_Abundance),
    Lower.CI = quantile(Posterior_Abundance, probs = 0.025),
    Upper.CI = quantile(Posterior_Abundance, probs = 0.975)
  )

#Visualise results
ggplot(total_by_year_summary, aes(x = Posterior_Sample, y = Mean)) +
  geom_line() +
  #geom_ribbon(aes(ymin = Lower.CI, ymax = Upper.CI), alpha = 0.2, fill = "blue") +
  theme_minimal() +
  labs(
    title = "Total Abundance of Chum Salmon Across BC",
    x = "Year",
    y = "Total Abundance"
  )
