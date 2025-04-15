setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)
source("DataProcessing/GetRegionalValue.R")
source("DataProcessing/GetRegionalEffort.R")
source("DataProcessing/GetRegionalLicenses.R") 
source("DataProcessing/GetTminusGenlength.R")

matchdat <- read.csv("../Data/Salmon/Catch & Effort/DispersalGrounds.csv") |> rename_with(str_to_lower)
pens <- read.csv("../Data/Salmon/Aquaculture/FarmsperRegion.csv") |> rename_with(str_to_lower) |> select(c(year,region,active.farms))
reccatch <- read.csv("../Data/Salmon/Catch & Effort/Recreational/YearlyRegionalRecCatch.csv") |> mutate(rec.catch=catch) |> 
  select(-c(X,catch))


addedfarms <- left_join(finalvalue, pens, 
                        by = c("region","year"))
addedAbund <- left_join(addedfarms, salmon, 
                      by = c("year", "species","region"))
addedeffort <- left_join(addedAbund, effort, 
                         by = c("year","region"))
addedrec <- left_join(addedeffort, reccatch,
                      by = c("year", "species","region"))

rm(list=setdiff(ls(), c("addedrec")))
write.csv(addedrec, "../Data/Salmon/RegionalEESVsNew2.csv")
