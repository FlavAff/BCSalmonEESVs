setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
#rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)
library(slider)

source("DataProcessing/GetRegionalGenLengths.R")
gens$genlength <- round(gens$genlength, digits=0)
source("DataProcessing/GetRegionalRelease.R") 
SST <- read.csv("../Data/dataset-satellite-sea-surface-temperature/Spring/yearlyMeanVarSpringSST.csv")
# salmon <- read_csv("../Data/Salmon/Ecological Supply/YearlyTotalRegionalAbundance.csv") |> rename_with(str_to_lower) |> 
#   rename(abundance = 'total.abundance')
salmon <- read_csv("../Data/Salmon/Ecological Supply/AverageAbundanceRegion.csv") |> rename_with(str_to_lower)
#salmon <- salmon[,c(2:5)]

finaldat <- left_join(salmon, gens,
                       by = c("region","species")) |> mutate(parental_abundance = NA)
finaldat <- left_join(finaldat, hatchery,
                      by = c("region","species","year")) |> mutate(parental_release = NA)
finaldat <- left_join(finaldat, SST,
                      by = c("year")) |> mutate(parental_sst = NA)

datout <- NULL
for (g in unique(finaldat$region)) {
  finaldat.r <- filter(finaldat, region == g)
  
  for (i in 1:nrow(finaldat.r)) {
    genlength <- finaldat.r$genlength[i]
    yearsago <- i - genlength
    
    if (!is.na(genlength) && yearsago >= 1) {  # Make sure yearsago is not negative
      mum <- finaldat.r$abundance[yearsago]
      mum2 <- finaldat.r$release[yearsago]
      if (finaldat.r$species[i] == "Sockeye") {
        mum3 <- finaldat.r$temperature[(yearsago+2)]
      }
      if (finaldat.r$species[i] != "Sockeye") {
        mum3 <- finaldat.r$temperature[(yearsago+1)]
      }
      
      
      if (!is.na(mum)) {
        finaldat.r$parental_abundance[i] <- mum
      }
      if (!is.na(mum2)) {
        finaldat.r$parental_release[i] <- mum2
      }
      if (!is.na(mum3)) {
        finaldat.r$parental_sst[i] <- mum3
      }
    }
  }
  datout <- rbind(finaldat.r,datout)
}

rm(finaldat,finaldat.r,gens,salmon,g,genlength,i,mum,mum2,mum3,yearsago,SST)
salmon <- datout |> filter(year>1995)

