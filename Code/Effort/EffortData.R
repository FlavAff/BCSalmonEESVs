#setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
#rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)
library(ggplot2)

dat <- read_csv("../Data/Salmon/Catch & Effort/ACEESV.csv")
dat$licenses[213:219] <- c(94, 153, 545, 351, 197, 80, 59) #note that some areas are missing the data
licenses <- read_csv("../Data/Salmon/Licenses/TotalLicenses_1998-2024.csv") |> rename_with(str_to_lower) |> rename(gear = type)

#Some values were removed in the summing due to missing data
# Define all possible areas
all_areas <- c("A", "B", "C", "D", "E", "F", "G", "H")
# Identify missing areas and add them
complete_data <- dat %>%
  complete(year, area = all_areas, fill = list(Value = NA, Gear = NA)) %>%
  mutate(gear = case_when(
    is.na(gear) & area == "B" ~ "Seine",
    is.na(gear) & area == "E" ~ "Gillnet",
    is.na(gear) & area == "H" ~ "Troll",
    TRUE ~ gear),
    gear.n = case_when(
      gear == "Seine" ~ 1,
      gear == "Gillnet" ~ 3,
      gear == "Troll" ~ 2,
      TRUE ~ gear.n),
    area.n = case_when(
      area == "B" ~ 2,
      area == "E" ~ 5,
      area == "H" ~ 8,
      TRUE ~ area.n)
  )
complete_data <- left_join(complete_data, licenses, by = c("year", "gear", "area"))
dat <- complete_data[, !colnames(complete_data) %in% c("licenses","...5")]
dat <- dat |> rename(licenses = count)
#rm(list=setdiff(ls(), "dat"))