setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)
library(zoo)
library(boot)
library(tidyr)
library(brms)

salmon <- read_csv("../Data/Salmon/Ecological Supply/COSEWIC-compilation-main/Output/CU_Spawner_Abund_20241115.csv")

Chum <- filter(salmon, Species =="Chum")
Coho <- filter(salmon, Species =="Coho")
Chinook <- filter(salmon, Species =="Chinook")
Pink.E <- filter(salmon, Species =="Pink (Even)") |> filter(Year %% 2 == 0)
Pink.O <- filter(salmon, Species =="Pink (Odd)") |> filter(Year %% 2 != 0)
Sockeye <- filter(salmon, Species =="Sockeye")

getfullTS <- function(fish,scale) {
  # Fit the model using brms
  fit <- brm(
    Spawner.Abundance ~ s(Year) + (1 | CUID),
    data = fish,
    family = lognormal(),
    chains = 4, cores = 4, iter = 2000, control = list(adapt_delta = 0.95)
  )
  # Add predictions to the data
  fish$Spawner.Abundance <- ifelse(
    is.na(fish$Spawner.Abundance),
    posterior_predict(fit, newdata = fish, re_formula = NA) %>% colMeans(),
    fish$Spawner.Abundance
  )
  if (scale == "total") {
  # Aggregate data
  posterior_totals <- fish %>%
      group_by(Year) %>%
      summarise(Total.Abundance = sum(Spawner.Abundance, na.rm = TRUE))
  }
  if (scale == "region") {
    # Aggregate data
  posterior_totals <- fish %>%
    group_by(Year, Region) %>%
    summarise(Total.Abundance = sum(Spawner.Abundance, na.rm = TRUE), .groups = 'drop')
  }
  return(posterior_totals)
}


ChumTS <- getfullTS(Chum,"total")
CohoTS <- getfullTS(Coho,"total")
ChinookTS <- getfullTS(Chinook,"total")
PinkE.TS <- getfullTS(Pink.E,"total")
PinkO.TS <- getfullTS(Pink.O,"total")
SockeyeTS <- getfullTS(Sockeye,"total")

totals <- rbind(ChumTS,CohoTS,ChinookTS,PinkE.TS,PinkO.TS,SockeyeTS)
write.csv(totals, "../Data/Salmon/Ecological Supply/YearlyTotalAbundance.csv")
#This dataset is then cleaned manually to remove any out of sample imputations




ChumTS <- getfullTS(Chum,"region")
CohoTS <- getfullTS(Coho,"region")
ChinookTS <- getfullTS(Chinook,"region")
PinkE.TS <- getfullTS(Pink.E,"region")
PinkO.TS <- getfullTS(Pink.O,"region")
SockeyeTS <- getfullTS(Sockeye,"region")

ChumTS$species <- "Chum"
CohoTS$species <- "Coho"
ChinookTS$species <- "Chinook"
PinkE.TS$species <- "Pink"
PinkO.TS$species <- "Pink"
SockeyeTS$species <- "Sockeye"

totals <- rbind(ChumTS,CohoTS,ChinookTS,PinkE.TS,PinkO.TS,SockeyeTS)
write.csv(totals, "../Data/Salmon/Ecological Supply/YearlyTotalRegionalAbundance.csv")
#This dataset is then cleaned manually to remove any out of sample imputations