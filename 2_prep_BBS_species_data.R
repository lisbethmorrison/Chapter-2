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
BBS_2018 <- read.csv("../Data/BBS_abund_data/BBS_data_2018.csv", header=TRUE)
BBS_2004_2017 <- read.csv("../Data/BBS_abund_data/BBS_data_2004_2017.csv", header=TRUE)
seed_spp_list <- read.csv("../Data/Trait data/seed_dispersing_spp_list_JAT.csv", header=TRUE)
invert_spp_list <- read.csv("../Data/Trait data/insectivores_spp_list_JAT.csv", header=TRUE)

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

## change some common names to merge with trait species list
BBS_final3$ENGLISH_NAME[ BBS_final3$ENGLISH_NAME == "Common Crossbill" ] <- "Crossbill"
BBS_final3$ENGLISH_NAME[ BBS_final3$ENGLISH_NAME == "Feral Pigeon" ] <- "Rock Dove"
BBS_final3$ENGLISH_NAME[ BBS_final3$ENGLISH_NAME == "Cetti's Warbler" ] <- "Cetti Warbler" ## apostrophe's don't merge

## merge seed dispersing species list (so we only have BBS data for seed dispersers)
seed_spp_list <- seed_spp_list[,-c(2,3,5)]
seed_spp_list <- seed_spp_list[seed_spp_list$yes_no==1,] ## 53 species which are seed dispersers
BBS_final4 <- merge(BBS_final3, seed_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_final4$ENGLISH_NAME)) ## 43 species which have BBS data

## merge invertivore species list (so we only have BBS data for pest control species)
invert_spp_list <- invert_spp_list[,-c(2,4)]
invert_spp_list <- invert_spp_list[invert_spp_list$yes_no==1,] ## 108 species which are pest controllers
invert_spp_list$species <- as.character(invert_spp_list$species)
invert_spp_list$species[ invert_spp_list$species == "Cetti’s Warbler" ] <- "Cetti Warbler" ## won't merge with apostrophe
BBS_final5 <- merge(BBS_final3, invert_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_final5$ENGLISH_NAME)) ## 74 species which have BBS data

## save seed dispersal files
write.csv(BBS_final4, file="../Data/Analysis_data/Seed dispersal/BBS_2004_2018_seed.csv", row.names=FALSE) ## 200 sites and 43 species from 2004-2018
spp_list2 <- BBS_final4[,c(1,2)] 
spp_list2 <- unique(spp_list2)
## save species list (43 species)
write.csv(spp_list2, file="../Data/Analysis_data/Seed dispersal/BBS_seed_species_list.csv", row.names=FALSE) ## 43 species

## save invertivore files
write.csv(BBS_final5, file="../Data/Analysis_data/Pest control/BBS_2004_2018_invert.csv", row.names=FALSE) ## 200 sites and 43 species from 2004-2018
spp_list3 <- BBS_final5[,c(1,2)] 
spp_list3 <- unique(spp_list3)
## save species list (73 species)
write.csv(spp_list3, file="../Data/Analysis_data/Pest control/BBS_invert_species_list.csv", row.names=FALSE) ## 73 species
