###########################################################
## Title: Prep BBS species data for analysis 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2019
##########################################################

rm(list=ls()) # clear R

#install.packages("dplyr", dep = TRUE)
library(dplyr)
library(tidyr)
library(ggplot2)

## read in data
BBS_2018 <- read.csv("../Data/BBS_abund_data/BBS_data_2018.csv", header=TRUE)
BBS_2004_2017 <- read.csv("../Data/BBS_abund_data/BBS_data_2004_2017.csv", header=TRUE)
## detecability data
seed_detect <- read.csv("../Data/BBS_abund_data/seed_dispersers_with_det_prob.csv", header=TRUE)
invert_detect <- read.csv("../Data/BBS_abund_data/insect_pred_with_det_prob.csv", header=TRUE)
#BBS_final <- read.csv("../Data/BBS_final.csv", header=TRUE)
#seed_list <- read.csv("../Data/site_species_seed_dispersers.csv", header=TRUE)
#invert_list <- read.csv("../Data/site_species_insect_pred.csv", header=TRUE)

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

## calculate sum across all distance categories for each transect
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

## remove mammals from dataset (they have numbers instead of characters as code)
BBS_final <- subset(BBS_final, !grepl("[0-9]", CBC_CODE))
BBS_final <-  BBS_final[!(is.na(BBS_final$CBC_CODE) | BBS_final$CBC_CODE==""), ] ## remove blanks (sheep)

## change some common names to merge with trait & detectability species list
BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Common Crossbill" ] <- "Crossbill"
BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Feral Pigeon" ] <- "Rock Dove"
BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Cetti's Warbler" ] <- "Cetti Warbler" ## apostrophe's don't merge
BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Pied/White Wagtail" ] <- "Pied Wagtail" ## apostrophe's don't merge
## change rock dove CBC code to DV, not FP
BBS_final$CBC_CODE[BBS_final$ENGLISH_NAME=='Rock Dove'] <- 'DV'

## fill in data with zeros - if there is no record for a speciecs in a site-year combination, the species was not recorded = zero
BBS_final$YEAR <- as.numeric(levels(BBS_final$YEAR))[BBS_final$YEAR]
BBS_final2 <- BBS_final %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  complete(ENGLISH_NAME, CBC_CODE, VISIT, YEAR = 2004:2018)
BBS_final2$TOT[is.na(BBS_final2$TOT)] <- 0
## remove OBS_DATE column - not required
BBS_final2 <- BBS_final2[,-6]
## this is the complete dataset filled with zeros where a species wasn't seen at a site but the site WAS visited
length(unique(BBS_final2$ENGLISH_NAME)) ## 208 species
length(unique(BBS_final2$GRIDREF)) ## 200 sites

## merge seed dispersing species list (so we only have BBS data for seed dispersers)
seed_spp_list <- seed_spp_list[,-c(2,3,5)]
seed_spp_list <- seed_spp_list[seed_spp_list$yes_no==1,] ## 53 species which are seed dispersers
BBS_seed <- merge(BBS_final2, seed_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_seed$ENGLISH_NAME)) ## 49 species which have BBS data

## merge invertivore species list (so we only have BBS data for pest control species)
invert_spp_list <- invert_spp_list[,-c(2,4)]
invert_spp_list <- invert_spp_list[invert_spp_list$yes_no==1,] ## 100 species which are pest controllers
invert_spp_list$species <- as.character(invert_spp_list$species)
invert_spp_list$species[ invert_spp_list$species == "Cetti’s Warbler" ] <- "Cetti Warbler" ## won't merge with apostrophe
BBS_invert <- merge(BBS_final2, invert_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_invert$ENGLISH_NAME)) ## 89 species which have BBS data

BBS_final3 <- merge(BBS_seed, BBS_invert, all=TRUE) ## merge both together
BBS_final3 <- BBS_final3[,-7]
## save file to send to Dario
write.csv(BBS_final3, file="../Data/BBS_abund_data/BBS_2004_2018_zeros.csv", row.names=FALSE)
length(unique(BBS_final3$ENGLISH_NAME)) ## 108 unique species

#### READ IN ABUNDANCE DATA WITH DETECTABILITY
rm(list=ls()) # clear R
BBS_final_det <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_zeros_with_det_prob.csv", header=TRUE) ## 108 species at 200 sites

## which species are completely missing detectability data?
dat <- BBS_final_det %>% 
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarize_all(.funs = funs('NA' = sum(is.na(.))))
## remove NAs
BBS_final_det_comp <- BBS_final_det[!is.na(BBS_final_det$det.prob), ]
length(unique(BBS_final_det_comp$GRIDREF))
length(unique(BBS_final_det_comp$YEAR))
length(unique(BBS_final_det_comp$ENGLISH_NAME)) ## 96 species - 12 species removed which don't have ANY detectability data
BBS_final_det_comp <- droplevels(BBS_final_det_comp)
sp_pre <- as.data.frame(unique(BBS_final_det$ENGLISH_NAME)) ## 108
sp_post <- as.data.frame(unique(BBS_final_det_comp$ENGLISH_NAME)) ## 96
sp <- merge(sp_pre$`unique(BBS_final_det$ENGLISH_NAME)`, sp_post$`unique(BBS_final_det_comp$ENGLISH_NAME)`)
matched <- intersect(sp_pre$`unique(BBS_final_det$ENGLISH_NAME)`, sp_post$`unique(BBS_final_det_comp$ENGLISH_NAME)`)
all <-  union(sp_pre$`unique(BBS_final_det$ENGLISH_NAME)`, sp_post$`unique(BBS_final_det_comp$ENGLISH_NAME)`)
non.matched <- all[!all %in% matched] ## 12 species with no detectability (9 of these have det data from Johnsone 2014)

##############################
## COMPLETE CASE DETECTABILITY 

## Must be done **BEFORE** I take the max across visits - different detectability values for early vs late depending on habitat
## re-calculate TOT weighted by detection
## TOT / det.prob
BBS_final_det_comp$TOT_det <- BBS_final_det_comp$TOT/BBS_final_det_comp$det.prob
## take max value across both visits (using TOT_det)
BBS_final_det_comp <- BBS_final_det_comp %>% group_by(YEAR, GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(TOT_det = max(TOT_det))
BBS_final_det_comp <- as.data.frame(BBS_final_det_comp)
head(BBS_final_det_comp)
length(unique(BBS_final_det_comp$ENGLISH_NAME)) ## 96 species
length(unique(BBS_final_det_comp$GRIDREF)) ## 200

## only include species with at all 15 years of data
BBS_final_det_comp2 <- BBS_final_det_comp %>% group_by(GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(no_years = n())
## remove site-species combos with less than 12 years of data
BBS_final_det_comp2 <- BBS_final_det_comp2[BBS_final_det_comp2$no_years>=12,]
length(unique(BBS_final_det_comp2$ENGLISH_NAME)) ## still 96 species
length(unique(BBS_final_det_comp2$GRIDREF)) ## 130 sites - 70 removed
## these 70 sites have no habitat data recorded for at least one year - so tot_det can't be calculated and whole site is removed

## merge back to match BBS data
BBS_final_det_comp3 <- merge(BBS_final_det_comp, BBS_final_det_comp2, by=c("ENGLISH_NAME","CBC_CODE", "GRIDREF"), all=FALSE)
length(unique(BBS_final_det_comp3$ENGLISH_NAME)) ## 96 species
length(unique(BBS_final_det_comp3$GRIDREF)) ## 130 sites

## split file back into seed disperser and insect predators
## to use for future analysis
seed_spp_list <- read.csv("../Data/Trait data/seed_dispersing_spp_list.csv", header=TRUE)
invert_spp_list <- read.csv("../Data/Trait data/insectivores_spp_list.csv", header=TRUE)

seed_spp_list <- seed_spp_list[,-c(2,3,5)]
seed_spp_list <- seed_spp_list[seed_spp_list$yes_no==1,] ## 53 species which are seed dispersers
BBS_seed_det_comp <- merge(BBS_final_det_comp3, seed_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_seed_det_comp$ENGLISH_NAME)) ## 45 species which have BBS data
length(unique(BBS_seed_det_comp$GRIDREF)) ## 130 sites

## merge invertivore species list (so we only have BBS data for pest control species)
invert_spp_list <- invert_spp_list[,-c(2,4)]
invert_spp_list <- invert_spp_list[invert_spp_list$yes_no==1,] ## 100 species which are pest controllers
invert_spp_list$species <- as.character(invert_spp_list$species)
invert_spp_list$species[ invert_spp_list$species == "Cetti’s Warbler" ] <- "Cetti Warbler" ## won't merge with apostrophe
BBS_invert_det_comp <- merge(BBS_final_det_comp3, invert_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_invert_det_comp$ENGLISH_NAME)) ## 80 species which have BBS data 
length(unique(BBS_invert_det_comp$GRIDREF)) ## 130 sites

## save seed dispersal files
write.csv(BBS_seed_det_comp, file="../Data/BBS_abund_data/BBS_2004_2018_seed_det_comp.csv", row.names=FALSE) ## 180 sites and 39 species from 2004-2018
## save insectivore files
write.csv(BBS_invert_det_comp, file="../Data/BBS_abund_data/BBS_2004_2018_invert_det_comp.csv", row.names=FALSE) ## 180 sites and 55 species from 2004-2018
##########################################################################

site_pre <- as.data.frame(unique(BBS_final_det$GRIDREF)) ## 200
site_post <- as.data.frame(unique(BBS_final_det_comp3$GRIDREF)) ## 130
sites <- merge(site_pre$`unique(BBS_final_det$GRIDREF)`, site_post$`uniqueBBS_final_det_comp3$GRIDREF)`)
matched <- intersect(site_pre$`unique(BBS_final_det$GRIDREF)`, site_post$`unique(BBS_final_det_comp3$GRIDREF)`)
all <-  union(site_pre$`unique(BBS_final_det$GRIDREF)`, site_post$`unique(BBS_final_det_comp3$GRIDREF)`)
non.matched <- all[!all %in% matched] ## 70 sites with missing habitat data (& therefore detectability) for at least 1 year

################################## 2 FILL IN MISSING DETECTABILITY ################################## 
######## create BBS_final_det with interpolated missing values
## out of the 12 species without det data - 9 of these can be found from Johnstone et al 2014 Bird Study
## read in data
rm(list=ls()) # clear R

BBS_final_det <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_zeros_with_det_prob.csv", header=TRUE) ## 108 species at 200 sites
missing_det <- read.csv("../Data/johnstone_det_prob.csv", header=TRUE)

## now fill in detectability NAs for the 9 species
sum(is.na(BBS_final_det$det.prob)) ## 19366 rows of NAs in det.prob (~7% missing data)
BBS_final_det_interpol <- merge(BBS_final_det, missing_det, by=c("ENGLISH_NAME"), all=TRUE) ## creates two det.prob columns
BBS_final_det_interpol <- BBS_final_det_interpol %>% mutate(det.prob = coalesce(det.prob.x, det.prob.y)) %>%
  select(ENGLISH_NAME, GRIDREF, YEAR, VISIT, CBC_CODE, TOT, det.prob) ## merge two det.prob columns together

### Now fill in site-species combinations with NAs using average non-missing values
## group_by species, fill NAs with average of non-NA values
BBS_final_det_interpol <- BBS_final_det_interpol %>% group_by(ENGLISH_NAME) %>% 
                      mutate(det.prob = ifelse(is.na(det.prob), mean(det.prob, na.rm = TRUE), det.prob))
## remove NaN values (these are for common scoter, crane and honey-buzzard which have NO detectability data anywhere)
BBS_final_det_interpol <- na.omit(BBS_final_det_interpol)

## re-calculate TOT weighted by detection
## TOT / det.prob
BBS_final_det_interpol$TOT_det <- BBS_final_det_interpol$TOT/BBS_final_det_interpol$det.prob
## take max value across both visits (using TOT_det)
BBS_final_det_interpol <- BBS_final_det_interpol %>% group_by(YEAR, GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(TOT_det = max(TOT_det))
BBS_final_det_interpol <- as.data.frame(BBS_final_det_interpol)
head(BBS_final_det_interpol)
length(unique(BBS_final_det_interpol$ENGLISH_NAME)) ## 105 species
length(unique(BBS_final_det_interpol$GRIDREF)) ## 200 sites

## check each species-site combo has 15 years of data
BBS_final_det_interpol2 <- BBS_final_det_interpol %>% group_by(GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(no_years = n())
## YES!

## split file back into seed disperser and insect predators
## to use for future analysis
seed_spp_list <- read.csv("../Data/Trait data/seed_dispersing_spp_list.csv", header=TRUE)
invert_spp_list <- read.csv("../Data/Trait data/insectivores_spp_list.csv", header=TRUE)

seed_spp_list <- seed_spp_list[,-c(2,3,5)]
seed_spp_list <- seed_spp_list[seed_spp_list$yes_no==1,] ## 53 species which are seed dispersers
BBS_seed_det_interpol <- merge(BBS_final_det_interpol, seed_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_seed_det_interpol$ENGLISH_NAME)) ## 48 species which have BBS data
length(unique(BBS_seed_det_interpol$GRIDREF)) ## 200 sites

## merge invertivore species list (so we only have BBS data for pest control species)
invert_spp_list <- invert_spp_list[,-c(2,4)]
invert_spp_list <- invert_spp_list[invert_spp_list$yes_no==1,] ## 100 species which are pest controllers
invert_spp_list$species <- as.character(invert_spp_list$species)
invert_spp_list$species[ invert_spp_list$species == "Cetti’s Warbler" ] <- "Cetti Warbler" ## won't merge with apostrophe
BBS_invert_det_interpol <- merge(BBS_final_det_interpol, invert_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_invert_det_interpol$ENGLISH_NAME)) ## 87 species which have BBS data 
length(unique(BBS_invert_det_interpol$GRIDREF)) ## 200 sites

## save seed dispersal files
write.csv(BBS_seed_det_interpol, file="../Data/BBS_abund_data/BBS_2004_2018_seed_det_interpol.csv", row.names=FALSE) ## 180 sites and 39 species from 2004-2018
## save insectivore files
write.csv(BBS_invert_det_interpol, file="../Data/BBS_abund_data/BBS_2004_2018_invert_det_interpol.csv", row.names=FALSE) ## 180 sites and 55 species from 2004-2018
##########################################################################


