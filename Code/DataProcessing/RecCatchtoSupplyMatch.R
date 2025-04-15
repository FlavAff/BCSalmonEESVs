setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)

source("GetRegionalProportions.R")
Allprops <- Allprops %>% 
  rename(year = Year,
         species = Species,
         mgmt_area = MGMT_AREA
  )
catch <- read.csv("../Data/Salmon/Catch & Effort/Recreational/YearlyRecCatch.csv") |> 
  rename_with(str_to_lower)
matchdat <- read.csv("../Data/Salmon/Catch & Effort/DispersalGrounds.csv") |> 
  rename_with(str_to_lower)

areas <- sort(unique(matchdat$mgmt_area))
multiareas <- sort(c(4,5,6,18,19,29,101,102,104,105,106,107,108,109,110,111,130))
singleareas <- areas[! areas %in% multiareas]
matchdatuni <- filter(matchdat, mgmt_area %in% singleareas)

catch$region <- NA
for (i in 1:nrow(matchdatuni)) {
  for (j in 1:nrow(catch)) {
    if (matchdatuni$mgmt_area[i] == catch$area[j]) {
      catch$region[j] <- matchdatuni$region[i]
    }
  }
}
catchuni <- filter(catch, region != "NA")


catch$mgmt_area <- catch$area
catchmulti <- left_join(Allprops, catch, 
                        by = c("year", "mgmt_area", "species")) |> 
  mutate(catch_prop = catch*reg_prop) |> filter(year>2011)


catch1 <- data_frame(year = catchuni$year, mgmt_area = catchuni$area,
                     catch = catchuni$catch, region = catchuni$region,
                     species = catchuni$species)
catch2 <- data_frame(year = catchmulti$year, mgmt_area = catchmulti$mgmt_area,
                     catch = catchmulti$catch_prop, region = catchmulti$region.x,
                     species = catchmulti$species)
catchfinal <- rbind(catch1,catch2)
rm(list=setdiff(ls(), c("catchfinal")))
write.csv(catchfinal, "../Data/Salmon/Catch & Effort/RegionalRecCatch.csv")




datout <- NULL
for (fish in unique(catchfinal$species)) {
  spp <- filter(catchfinal, species == fish)
  for (i in c(2012:2022)) {
    year <- filter(spp, year == i)
    out <- data_frame(catch = tapply(year$catch, year$region, sum, na.rm=T),
                      region = sort(unique(year$region)))
    out$species <- fish
    out$year <- i
    datout <- bind_rows(datout,out)
  }
}
write.csv(datout, "../Data/Salmon/Catch & Effort/Recreational/YearlyRegionalRecCatch.csv")
