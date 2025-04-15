catch <- read_csv("../Data/Salmon/Catch & Effort/All_1996-2023.csv") |> rename_with(str_to_lower) 
catch <- catch %>%
  group_by(year) %>%
  summarise(Coho = sum(coho_kept, na.rm = T), Chum = sum(chum_kept, na.rm = T), 
            Sockeye = sum(sockeye_kept, na.rm = T), Pink = sum(pink_kept, na.rm = T), Chinook = sum(chinook_kept, na.rm = T))
catch <- tibble(year = rep(seq(1996,2023,1),5),
                species = c(rep("Coho",28),rep("Chum",28),rep("Pink",28),rep("Sockeye",28),rep("Chinook",28)),
                catch = c(catch$Coho,catch$Chum,catch$Pink,catch$Sockeye,catch$Chinook))


prices.wild <- read.csv("../Data/Salmon/Demand/Prices_1997-2023.csv") |> rename_with(str_to_lower) |> filter(source == "Wild commercial")
prices.farm <- read.csv("../Data/Salmon/Demand/Prices_1997-2023.csv") |> rename_with(str_to_lower) |> filter(source == "Farmed") |>
  filter(species == "Farmed")
prices.farm <- prices.farm |> rename(price.per.kilo.farmed.wholesale = price.per.kilo.wholesale) |> 
  rename(price.per.kilo.farmed.landed = price.per.kilo.landed)
prices.farm <- prices.farm[,c(1,8,9)]
disposable.income <- read.csv("../Data/Salmon/Demand/36100112.csv") |> rename_with(str_to_lower) |> 
  filter(estimates == "Equals: household disposable income") |> filter(year > 1995) |> 
  filter(seasonal.adjustment == "Unadjusted")
disposable.income <- tibble(year = seq(from=1996,to=2023,by=1),
                            disposable.income = tapply(disposable.income$value, disposable.income$year, mean))

dat1 <- left_join(prices.wild, prices.farm, by = "year")
dat2 <- left_join(disposable.income, catch, by = "year")
dat <- left_join(dat2, dat1, by = c("year","species"))

spp_mapping <- c("Sockeye" = 1, "Pink" = 2, "Chum" = 3, "Coho" = 4, "Chinook" = 5)
dat$spp.n <- spp_mapping[as.character(dat$species)]
dat$catch <- dat$catch/1000
dat$disposable.income <- dat$disposable.income/1000
rm(list=setdiff(ls(), c("dat","finaldat")))
dat <- dat %>% arrange(species, year)
