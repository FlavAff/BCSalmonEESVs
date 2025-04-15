setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)

catchdat <- read.csv("../Data/Salmon/Catch & Effort/RegionalCatch.csv")

datout <- NULL
for (fish in unique(catchdat$species)) {
  spp <- filter(catchdat, species == fish)
  for (i in c(1996:2023)) {
    year <- filter(spp, year == i)
    out <- data_frame(catch = tapply(year$catch, year$region, mean, na.rm=T),
                      region = sort(unique(year$region)))
    out$species <- fish
    out$year <- i
    datout <- bind_rows(datout,out)
  }
}

datout2 <- NULL
for (fish in unique(datout$species)) {
  spp <- filter(datout, species == fish)
  for (i in c(1996:2023)) {
    year <- filter(spp, year == i)
    for (l in c(1:nrow(year))) {
      year$catch_prop[l] <- as.numeric(year$catch[l]/sum(year$catch, na.rm = T))
    }
    datout2 <- bind_rows(datout2,year)
  }
}

value <- read_csv("../Data/Salmon/Landings/DetailedLandings_1997-2021.csv") |> rename_with(str_to_lower) |> 
  filter(species.group == "Salmon")

valuecatch <- left_join(datout2, value, 
                        by = c("species", "year"))
valuecatch$landed <- valuecatch$catch_prop*valuecatch$`landed value`
valuecatch$wholesale <- valuecatch$catch_prop*valuecatch$`wholesale value`

finalvalue <- valuecatch[,c(1,2,3,4,11,12)]
rm(list=setdiff(ls(), c("finalvalue")))
finalvalue[finalvalue == "NaN"] <- NA
