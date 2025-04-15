setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
#rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)

hatcheries <- read_csv("../Data/Salmon/Ecological Supply/hatchery_releases_2024.csv")

hatchery <- NULL
for (fish in unique(hatcheries$species)) {
  spp <- filter(hatcheries, species == fish)
  for (i in c(1956:2023)) {
    year <- filter(spp, year == i)
    out <- data_frame(release = tapply(year$total_release, year$region, sum, na.rm=T),
                      region = sort(unique(year$region)))
    out$species <- fish
    out$year <- i
    hatchery <- bind_rows(hatchery,out)
  }
}

rm(hatcheries,fish,i,spp,out,year)
