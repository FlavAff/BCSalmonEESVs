setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
#rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)

licenses <- read_csv("../Data/Salmon/Licenses/TotalLicenses_1998-2024.csv") |> rename_with(str_to_lower)
areamatch <- read_csv("../Data/Salmon/Licenses/AreaToRegionMatch.csv") |> rename_with(str_to_lower)

license.per.region <- left_join(licenses, areamatch, 
                                by = c("area"))


licenses <- NULL
for (reg in unique(license.per.region$region)) {
  regio <- filter(license.per.region, region == reg)
  out <- data_frame(licenses = tapply(regio$count, regio$year, sum, na.rm=T),
                    region = sort(unique(regio$region)),
                    year = seq(1998,2024,1))
  licenses <- bind_rows(licenses,out)
}

rm(license.per.region,areamatch,reg,out,regio)
