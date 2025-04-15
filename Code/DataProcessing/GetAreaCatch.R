setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)

catch <- read_csv("../Data/Salmon/Catch & Effort/All_1996-2023.csv")


sumso <- list()
sumpi <- list()
sumchi <- list()
sumco <- list()
sumchu <- list()
for (i in unique(catch$year)) {
  cal <- filter(catch, year == i)
  socksum <- list(tapply(cal$SOCKEYE_KEPT, cal$MGMT_AREA, sum, na.rm=T))
  pinksum <- list(tapply(cal$PINK_KEPT, cal$MGMT_AREA, sum, na.rm=T))
  chinsum <- list(tapply(cal$CHINOOK_KEPT, cal$MGMT_AREA, sum, na.rm=T))
  cohosum <- list(tapply(cal$COHO_KEPT, cal$MGMT_AREA, sum, na.rm=T))
  chumsum <- list(tapply(cal$CHUM_KEPT, cal$MGMT_AREA, sum, na.rm=T))
  
  sumso <- c(sumso,socksum)
  sumpi <- c(sumpi,pinksum)
  sumchi <- c(sumchi,chinsum)
  sumco <- c(sumco,cohosum)
  sumchu <- c(sumchu,chumsum)
}


listtodf <- function(dat) {
  Catch <- NULL
  for (i in 1:28) {
    Catch1 <- c(dat[[i]])
    Catch2 <- as.data.frame(Catch1)
    Catch2$MGMT_AREA <- rownames(Catch2)
    Catch2$YEAR <- i + 1995
    Catch <- rbind(Catch, Catch2)
  }
  return(Catch)
}
sockeyecatch <- listtodf(sumso)
sockeyecatch$Species <- rep("Sockeye",nrow(sockeyecatch))
pinkcatch <- listtodf(sumpi)
pinkcatch$Species <- rep("Pink",nrow(pinkcatch))
cohocatch <- listtodf(sumco)
cohocatch$Species <- rep("Coho",nrow(cohocatch))
chumcatch <- listtodf(sumchu)
chumcatch$Species <- rep("Chum",nrow(chumcatch))
chinookcatch <- listtodf(sumchi)
chinookcatch$Species <- rep("Chinook",nrow(chinookcatch))

sumcatch <- rbind(sockeyecatch,pinkcatch,cohocatch,
                  chumcatch,chinookcatch)

write.csv(sumcatch, "../Data/Salmon/Catch & Effort/Catchperarea.csv")
