setwd("~/Documents/McGill/PhD/Chapter 2/Code/")

catchdat <- read_csv("../Data/Salmon/AreaEESVs.csv")
catchdat <- catchdat %>%
  group_by(year, species) %>%
  summarise(catch = sum(catch, na.rm = TRUE)) %>%
  ungroup()
dat <- read_csv("../Data/Salmon/Landings/DetailedLandings_1997-2023.csv") |> rename_with(str_to_lower) |> filter(source == "Wild commercial")|> 
  filter(species.group == "Salmon") |> filter(harvest != "NA") |> rename(landed.value = 'landed value') |> rename(wholesale.value = 'wholesale value')
dat <- left_join(dat, catchdat, by = c("year","species"))
dat$catch <- dat$catch / 100000 #make catch smaller by 100,000 to scale closer to units of harvest (tonnes) and value (million $)
dat <- dat %>%
  mutate(spp.n = case_when(
    species == "Sockeye" ~ 1,
    species == "Chinook" ~ 2,
    species == "Coho"    ~ 3,
    species == "Chum"    ~ 4,
    species == "Pink"    ~ 5,
    TRUE ~ NA_real_  # Assign NA for any species not listed
  ))
