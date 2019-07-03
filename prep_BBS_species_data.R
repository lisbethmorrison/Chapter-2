###########################################################
## Title: Prep BBS species data for analysis 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2019
##########################################################

rm(list=ls()) # clear R

#install.packages("dplyr", dep = TRUE)
library(dplyr)

## read in data
BBS_2018 <- read.csv("../Data/BBS_data_2018.csv", header=TRUE)
BBS_2004_2017 <- read.csv("../Data/BBS_data_2004_2017.csv", header=TRUE)

BBS_2018$YEAR <- as.factor(BBS_2018$YEAR)

### first calculate TOT for 2004-2017 data
## this is summed across all distance categories and transect sections for each year and visit
## first change NAs to zeros (for some reason zeros have changed to NAs)
BBS_2004_2017[is.na(BBS_2004_2017)] <- 0

str(BBS_2004_2017)
BBS_2004_2017$ï..CBC_CODE <- as.character(BBS_2004_2017$ï..CBC_CODE)
BBS_2004_2017$ENGLISH_NAME <- as.character(BBS_2004_2017$ENGLISH_NAME)
BBS_2004_2017$YEAR <- as.factor(BBS_2004_2017$YEAR)
BBS_2004_2017$BAND_1 <- as.numeric(BBS_2004_2017$BAND_1)
BBS_2004_2017$BAND_2 <- as.numeric(BBS_2004_2017$BAND_2)
BBS_2004_2017$BAND_3 <- as.numeric(BBS_2004_2017$BAND_3)
BBS_2004_2017$FLYING <- as.numeric(BBS_2004_2017$FLYING)

## calculate sum across all distance categories
BBS_2004_2017$SUM <- rowSums(BBS_2004_2017[,8:11])

## calculate TOT across all transects
BBS_2004_2017_tot <- BBS_2004_2017 %>%
  group_by(ï..CBC_CODE, ENGLISH_NAME, YEAR, OBS_DT, VISIT, GRIDREF) %>% 
  summarise(TOT = sum(SUM))
BBS_2004_2017_tot <- as.data.frame(BBS_2004_2017_tot)
head(BBS_2004_2017_tot)
## change first column to CBC_CODE
names(BBS_2004_2017_tot)[1]<-paste("CBC_CODE")
## remove some columns not needed (tenkm and region_code) from 2018 dataset
BBS_2018 <- BBS_2018[, -c(2:3)]
head(BBS_2018)

## merge together
BBS_final <- rbind(BBS_2004_2017_tot, BBS_2018)
## take max value across both visits
BBS_final <- BBS_final %>% group_by(YEAR, GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(TOT = max(TOT))
BBS_final <- as.data.frame(BBS_final)
head(BBS_final)

## remove mammals from dataset (they have numbers instead of characters as code)
BBS_final <- subset(BBS_final, !grepl("[0-9]", CBC_CODE))
## remove rows with zero TOT counts (this removes sheep and some house martin and rook records)
BBS_final<-BBS_final[!(BBS_final$TOT==0),]

## only include species with at least 8 years of data (not continuous)
BBS_final2 <- BBS_final %>% group_by(GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(no_years = n())
## remove site-species combos with less than 8 years of data
BBS_final2 <- BBS_final2[BBS_final2$no_years>=8,]
length(unique(BBS_final2$ENGLISH_NAME)) ## 124 species
length(unique(BBS_final2$GRIDREF)) ## 200 sites still

## merge species and sites into BBS_final
BBS_final3 <- merge(BBS_final2, BBS_final, by=c("ENGLISH_NAME","CBC_CODE", "GRIDREF"), all=FALSE)
length(unique(BBS_final3$ENGLISH_NAME)) ## 124 species
length(unique(BBS_final3$GRIDREF)) ## 200 sites still

## save BBS_final
write.csv(BBS_final3, file="../Data/BBS_2004_2018.csv", row.names=FALSE) ## 200 sites and 209 species from 2004-2018
spp_list <- BBS_final3[,c(1,2)]
spp_list <- unique(spp_list)
## save species list
write.csv(spp_list, file="../Data/BBS_species_list.csv", row.names=FALSE)
