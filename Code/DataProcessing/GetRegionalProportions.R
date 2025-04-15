setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)
library(tidyr)

salmon <- read.csv("~/Documents/McGill/PhD/Chapter 2/Data/Salmon/Ecological Supply/YearlyTotalRegionalAbundance.csv")
salmon <- salmon[,c(2:5)]

# Pivot wider to align data correctly
data_wide <- salmon %>%
  pivot_wider(
    names_from = Region,
    values_from = Abundance,
    values_fill = list(Abundance = NA) # Fill missing combinations with NA
  )

chum <- data_wide |> filter(species == "Chum")
coho <- data_wide |> filter(species == "Coho")
chinook <- data_wide |> filter(species == "Chinook")
pinkE <- data_wide |> filter(species == "Pink (Even)")
pinkO <- data_wide |> filter(species == "Pink (Odd)")
sockeye <- data_wide |> filter(species == "Sockeye")

areas <- tribble(
  ~area, ~regions,
  "4", c("Nass", "Skeena"),
  "5", c("Skeena", "CC"),
  "104", c("Nass", "Skeena", "CC"),

  "101", c("Nass", "Skeena", "CC", "HG"),
  "102", c("Nass", "CC", "HG","Skeena"),
  "105", c("Nass", "Skeena", "CC"),
  "106", c("Skeena", "CC"),
  
  "111", c("CC", "VIMI","Fraser"),
  "11", c("VIMI","Fraser"),
  "12", c("VIMI","Fraser"),
  "13", c("VIMI","Fraser"),
  "14", c("VIMI","Fraser"),
  "15", c("VIMI","Fraser"),
  "16", c("VIMI","Fraser"),
  "17", c("VIMI","Fraser"),
  "18", c("Fraser", "VIMI"),
  "19", c("Fraser", "VIMI"),
  "29", c("Fraser", "VIMI"),
  "20", c("Fraser", "VIMI"),
  "121", c("Fraser", "VIMI"),
  "123", c("Fraser", "VIMI"),
  "124", c("Fraser", "VIMI"),
  "125", c("Fraser", "VIMI"),
  "126", c("Fraser", "VIMI"),
  "127", c("Fraser", "VIMI"),
  
  "130", c("Nass", "Skeena", "CC", "HG", "VIMI","Fraser")
)

getRegionProps <- function(area, regions, abund_table) {
  prop_table <- map(regions,
                    ~ abund_table |>
                      rowwise() |>
                      transmute("{.x}_prop" := .data[[.x]] / 
                                  sum(c_across(regions),
                                      na.rm = TRUE))
  ) |> list_cbind()
  
  bind_cols(abund_table, prop_table) %>%
    mutate(area = !!area, .after = Year) |>
    pivot_longer(cols = ends_with("_prop"),
                 names_to = "region",
                 values_to = "reg_prop") |>
    mutate(region = region |> str_remove("_prop"))
}


getprops <- function(fish) {
  props <- pmap(areas, getRegionProps, abund_table=fish) %>%
    list_rbind()
  props$MGMT_AREA <- as.numeric(props$area)
  propsonly <- props[c("Year","MGMT_AREA","region","reg_prop","species")]
  return(propsonly)
}

Chumprops <- getprops(chum)
Cohoprops <- getprops(coho)
Chinookprops <- getprops(chinook)
PinkEprops <- getprops(pinkE)
PinkOprops <- getprops(pinkO)
Sockeyeprops <- getprops(sockeye)

Pinkprops <- PinkEprops[,c(1,2,3,5)]
Pinkprops$reg_prop <- ifelse(is.na(PinkEprops$reg_prop), PinkOprops$reg_prop, PinkEprops$reg_prop)

Allprops <- rbind(Chumprops,Cohoprops,Chinookprops,Pinkprops,Sockeyeprops)
rm(list=setdiff(ls(), c("Allprops")))
