setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
#rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)

effort <- read_csv("../Data/Salmon/Catch & Effort/All_1996-2023.csv") |> rename_with(str_to_lower)
matchdat <- read.csv("../Data/Salmon/Catch & Effort/DispersalGrounds.csv") |> 
  rename_with(str_to_lower)

match <- left_join(effort, matchdat, 
                   by = c("mgmt_area"))

effort <- match[,c(4,6,7,15)]


datout <- NULL
for (i in c(1996:2023)) {
    year <- filter(effort, year == i)
    days <- data_frame(boat_days = tapply(year$boat_days, year$region, sum, na.rm=T),
                      region = sort(unique(year$region)))
    vessels <- data_frame(vessel_count = tapply(year$vessel_count, year$region, sum, na.rm=T))
    out <- cbind(days,vessels)
    out$year <- i
    datout <- bind_rows(datout,out)
}

for (i in 1:nrow(datout)) {
  if (datout$vessel_count[i] == 0) {
    datout$vessel_count[i] <- NA
  }
}

effort <- datout
rm(days,datout,match,matchdat,out,vessels,year,i)
