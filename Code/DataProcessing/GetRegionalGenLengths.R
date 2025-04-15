setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
#rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)

generations <- read_csv("../Data/Salmon/Ecological Supply/GenLengths.csv") |> rename_with(str_to_lower) |>
  filter(species != "NA")
generations <- generations[c(1:400),]

gens <- NULL
for (reg in unique(generations$region)) {
  regio <- filter(generations, region == reg)
  out <- data_frame(genlength = tapply(regio$gen_length, regio$species, mean, na.rm=T),
                    region = sort(unique(regio$region)),
                    species = sort(unique(regio$species)))
  gens <- bind_rows(gens,out)
}

rm(generations,out,regio,reg)
