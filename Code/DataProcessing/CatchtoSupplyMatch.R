setwd("~/Documents/McGill/PhD/Chapter 2/Code/")
rm(list = ls(all=TRUE)); #Remove all the objects in the memory

library(tidyverse)

source("GetRegionalProportions.R")
Allprops <- Allprops %>% 
                rename(YEAR = Year,
                       Region = region
  )
catch <- read.csv("../Data/Salmon/Catch & Effort/Catchperarea.csv")
catch <- catch[,c(2:5)]
matchdat <- read.csv("../Data/Salmon/Catch & Effort/DispersalGrounds.csv")

areas <- sort(unique(matchdat$MGMT_AREA))
multiareas <- sort(c(4,5,11,12,13,14,15,16,17,18,19,20,29,101,102,104,105,106,111,121,123,124,125,126,127,130))
singleareas <- areas[! areas %in% multiareas]
matchdatuni <- filter(matchdat, MGMT_AREA %in% singleareas)

catch$Region <- NA
for (i in 1:nrow(matchdatuni)) {
  for (j in 1:nrow(catch)) {
    if (matchdatuni$MGMT_AREA[i] == catch$MGMT_AREA[j]) {
      catch$Region[j] <- matchdatuni$Region[i]
    }
  }
}
catchuni <- filter(catch, Region != "NA")

catchmulti <- left_join(Allprops, catch, 
                              by = c("YEAR", "MGMT_AREA", "species")) |> 
  mutate(catch_prop = Catch1*reg_prop) |> filter(YEAR > 1995)


catch1 <- data_frame(year = catchuni$YEAR, mgmt_area = catchuni$MGMT_AREA,
                     catch = catchuni$Catch1, region = catchuni$Region,
                     species = catchuni$species)
catch2 <- data_frame(year = catchmulti$YEAR, mgmt_area = catchmulti$MGMT_AREA,
                     catch = catchmulti$catch_prop, region = catchmulti$Region,
                     species = catchmulti$species)
catchfinal <- rbind(catch1,catch2)
rm(list=setdiff(ls(), c("catchfinal")))
write.csv(catchfinal, "../Data/Salmon/Catch & Effort/RegionalCatch.csv")


###############################################################################
############## OLD CODE
###############################################################################

#rm(list=setdiff(ls(), c("salmon","catch")))
#matchdat <- read.csv("../Data/Salmon/Catch & Effort/DispersalGrounds.csv")

#salmon$Nass.SC <- NA
#salmon$Skeena.NC <- NA
#salmon$CC.SN <- NA
#salmon$Skeena.C <- NA
#salmon$CC.S <- NA
#salmon$Fraser.V <- NA
#salmon$VIMI.F <- NA
#salmon$Nass.SCH <- NA
#salmon$Skeena.NCH <- NA
#salmon$CC.SNH <- NA
#salmon$HG.SCN <- NA
#salmon$Nass.CH <- NA
#salmon$CC.NH <- NA
#salmon$HG.NC <- NA
#salmon$Skeena.CV <- NA
#salmon$CC.SV <- NA
#salmon$VIMI.SC <- NA
#salmon$Skeena.NCHV <- NA
#salmon$Nass.SCHV <- NA
#salmon$CC.NSHV <- NA
#salmon$HG.NCSV <- NA
#salmon$VIMI.NCHs <- NA

#for (i in 1:nrow(salmon)) {
#  salmon$Nass.SC[i] <- salmon$Nass[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i])
#  salmon$Skeena.NC[i] <- salmon$Skeena[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i])
#  salmon$CC.SN[i] <- salmon$CC[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i])
#  
#  salmon$Skeena.C[i] <- salmon$Skeena[i]/(salmon$Skeena[i]+salmon$CC[i])
#  salmon$CC.S[i] <- salmon$CC[i]/(salmon$Skeena[i]+salmon$CC[i])
#  
#  salmon$Fraser.V[i] <- salmon$Fraser[i]/(salmon$Fraser[i]+salmon$VIMI[i])
#  salmon$VIMI.F[i] <- salmon$VIMI[i]/(salmon$Fraser[i]+salmon$VIMI[i])
#  
#  salmon$Nass.SCH[i] <- salmon$Nass[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i]+salmon$HG[i])
#  salmon$Skeena.NCH[i] <- salmon$Skeena[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i]+salmon$HG[i])
#  salmon$CC.SNH[i] <- salmon$HG[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i]+salmon$HG[i])
#  salmon$HG.SCN[i] <- salmon$CC[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i]+salmon$HG[i])
#  
#  salmon$Nass.CH[i] <- salmon$Nass[i]/(salmon$Nass[i]+salmon$CC[i]+salmon$HG[i])
#  salmon$CC.NH[i] <- salmon$CC[i]/(salmon$Nass[i]+salmon$CC[i]+salmon$HG[i])
#  salmon$HG.NC[i] <- salmon$HG[i]/(salmon$Nass[i]+salmon$CC[i]+salmon$HG[i])
#  
#  salmon$Skeena.CV[i] <- salmon$Skeena[i]/(salmon$Skeena[i]+salmon$VIMI[i]+salmon$CC[i])
#  salmon$CC.SV[i] <- salmon$CC[i]/(salmon$Skeena[i]+salmon$VIMI[i]+salmon$CC[i])
#  salmon$VIMI.SC[i] <- salmon$VIMI[i]/(salmon$Skeena[i]+salmon$VIMI[i]+salmon$CC[i])
#  
#  salmon$Skeena.NCHV[i] <- salmon$Skeena[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i]+salmon$HG[i]+salmon$VIMI[i])
#  salmon$Nass.SCHV[i] <- salmon$Nass[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i]+salmon$HG[i]+salmon$VIMI[i])
#  salmon$CC.NSHV[i] <- salmon$CC[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i]+salmon$HG[i]+salmon$VIMI[i])
#  salmon$HG.NCSV[i] <- salmon$HG[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i]+salmon$HG[i]+salmon$VIMI[i])
#  salmon$VIMI.NCHs[i] <- salmon$VIMI[i]/(salmon$Skeena[i]+salmon$Nass[i]+salmon$CC[i]+salmon$HG[i]+salmon$VIMI[i])
#}

#sockeyeS <- filter(salmon, Species == "Sockeye")
#sockeyeC <- filter(catch, Species == "Sockeye")

#for (i in unique(sockeyeS$Year)) {
#  yearS <- subset(sockeyeS, sockeyeS$Year==i)
#  yearC <- subset(sockeyeC, sockeyeC$YEAR==i)
#  for (j in multiareas) {
#    area <- subset(yearC, yearC$MGMT_AREA==j)
#    if (j == 18) {
#      skeena18 <- area$Catch1*yearS$Fraser.V
#      nass18 <- area$Catch1*yearS$VIMI.F
#      #print(sockeyeC$Catch1)
#      print(skeena18+nass18)
#    }
#  }
#}


#getcatchperarea <- function(dat1,dat2,area) {
#  areaC <- subset(dat2, dat2$MGMT_AREA==area)
#  df <- NULL
#  
# for (i in c(1996:2018)) {
#    yearC <- subset(areaC, areaC$YEAR==i)
#    yearS <- subset(dat1, dat1$Year==i)
#    
#    if (area %in% c(18,19,29)) {
#      skeena <- yearC$Catch1*yearS$Fraser.V
#      vimi <- yearC$Catch1*yearS$VIMI.F
#      areaSV <- c(skeena,vimi)
#      if (length(areaSV) > 0) {
#        dfSV <- data_frame(Year = rep(i,2), MGMT_AREA = rep(area,2), Catch = c(skeena,vimi), Region = c("Skeena","VIMI"))
#        df <- rbind(df,dfSV)
#      }
#    }
#    if (area %in% c(4,5,104)) {
#      skeena <- yearC$Catch1*yearS$Skeena.NC
#      nass <- yearC$Catch1*yearS$Nass.SC
#      cc <- yearC$Catch1*yearS$CC.SN
#      areaNSC <- c(nass,skeena,cc)
#      if (length(areaNSC) > 0) {
#        dfNSC <- data_frame(Year = rep(i,3), MGMT_AREA = rep(area,3), Catch = c(skeena,nass,cc), Region = c("Skeena","Nass","CC"))
#        df <- rbind(df,dfNSC)
#      }
#    }
#  }
#  return(df)
#}

#area18 <- getcatchperarea(sockeyeS,sockeyeC,5)