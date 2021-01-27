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

#### ADD IN DETECTABILITY
## Must be done **BEFORE** I take the max across visits - different detectability values for early vs late depending on habitat

### Merge seed and invert detectability (should have 200 sites & 74 species)
bbs_detect <- rbind(seed_detect, invert_detect)
bbs_detect <- unique(bbs_detect)
length(unique(bbs_detect$ENGLISH_NAME)) ## 74 species
length(unique(bbs_detect$GRIDREF)) ## 200 sites

################################## 1 COMPLETE CASES ################################## 
### merge in with BBS data => COMPLETE CASES (NAs are removed)
bbs_detect2 <- na.omit(bbs_detect)
length(unique(bbs_detect2$ENGLISH_NAME)) ## 72 species (chough and rock pipit don't have detectability data)
## some site-species combinations will also be removed (rook, stock dove, house martin, sand martin, swallow and swift)
length(unique(bbs_detect2$GRIDREF)) ## 200 sites
## Remove 2nd column (ENGLISH_NAME2)
bbs_detect2 <- bbs_detect2[,-2]
## BBS_final_detect1 = complete cases
## rename columns in bbs_detect2
colnames(bbs_detect2) <- c("GRIDREF", "ENGLISH_NAME", "CBC_CODE", "YEAR", "VISIT", "det.prob")
bbs_detect2$YEAR <- as.factor(bbs_detect2$YEAR)
## merge two together (72 species, 200 sites)
BBS_final_detect1 <- left_join(bbs_detect2, BBS_final, by=c("ENGLISH_NAME", "CBC_CODE", "GRIDREF", "YEAR", "VISIT"))
BBS_final_detect1 <- na.omit(BBS_final_detect1)
## for some reason there is detection probabilities for site-species-visit combinations where there is no count data
## this file is smaller than BBS_final because 6681 rows (site-species combos) from detectability
## are missing due to missing habitat data

## re-calculate TOT weighted by detection
## TOT / det.prob
BBS_final_detect1$TOT_det <- BBS_final_detect1$TOT/BBS_final_detect1$det.prob

## take max value across both visits (using TOT_det)
BBS_final_detect1 <- BBS_final_detect1 %>% group_by(YEAR, GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(TOT_det = max(TOT_det))
BBS_final_detect1 <- as.data.frame(BBS_final_detect1)
head(BBS_final_detect1)
length(unique(BBS_final_detect1$ENGLISH_NAME)) ## 72 species
length(unique(BBS_final_detect1$GRIDREF)) ## 200

## only include species with at least 12 years of data (non-consecutive)
BBS_final2 <- BBS_final_detect1 %>% group_by(GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(no_years = n())
## remove site-species combos with less than 12 years of data
BBS_final2 <- BBS_final2[BBS_final2$no_years>=12,]
length(unique(BBS_final2$ENGLISH_NAME)) ## 70 species (2 species removed)
## grey wagtail and whinchat removed --> when no detectability, they were only found at one site each
## with detectability --> these sites are missing a few years of detectability, so are removed completely
length(unique(BBS_final2$GRIDREF)) ## 180 sites (20 sites removed)
## most likely because these 20 sites are missing habitat data for some years, therefore detectability
## cannot be calculated and results in the whole site being removed
table(BBS_final2$GRIDREF)

## merge species and sites into BBS_final
BBS_final3 <- merge(BBS_final2, BBS_final_detect1, by=c("ENGLISH_NAME","CBC_CODE", "GRIDREF"), all=FALSE)
length(unique(BBS_final3$ENGLISH_NAME)) ## 70 species
length(unique(BBS_final3$GRIDREF)) ## 180 sites still
## 52,120 rows
##########################################################################

################################## 2 FILL IN MISSING DETECTABILITY ################################## 
######## create BBS_final_detect2 with interpolated missing values
## Chough and Rock pipit - need all values interpolated (only occur at one site after filtering)
## Get this data from Johnstone et al 2014 Bird Study
chough_det <- 0.804
rock_pipit_det <- 0.427

bbs_detect <- bbs_detect[,-2] ## remove duplicated english_name column
colnames(bbs_detect) <- c("GRIDREF", "ENGLISH_NAME", "CBC_CODE", "YEAR", "VISIT", "det.prob") ## rename columns
BBS_final_detect2 <- merge(bbs_detect, BBS_final, by=c("ENGLISH_NAME", "CBC_CODE", "GRIDREF", "YEAR", "VISIT"),
                           all=TRUE) ## merge detectability with BBS_final
## remove NAs in TOT column (for some reason there is detectability where there is no count data)
BBS_final_detect2 <- BBS_final_detect2[!is.na(BBS_final_detect2$TOT), ] ## same length at BBS_final
## now fill in detectability NAs for Chough and Rock pipit
sum(is.na(BBS_final_detect2$det.prob)) ## 60,609 rows of NAs in det.prob (~37% missing data)
BBS_final_detect2$det.prob[BBS_final_detect2$ENGLISH_NAME=="Chough" & is.na(BBS_final_detect2$det.prob)] <- chough_det
BBS_final_detect2$det.prob[BBS_final_detect2$ENGLISH_NAME=="Rock Pipit" & is.na(BBS_final_detect2$det.prob)] <- rock_pipit_det
### Now fill in site-species combinations with NAs using average non-missing values
## group_by species, fill NAs with average of non-NA values
BBS_final_detect2 <- BBS_final_detect2 %>% group_by(ENGLISH_NAME) %>% 
                      mutate(det.prob = ifelse(is.na(det.prob), mean(det.prob, na.rm = TRUE), det.prob))

## re-calculate TOT weighted by detection
## TOT / det.prob
BBS_final_detect2$TOT_det <- BBS_final_detect2$TOT/BBS_final_detect2$det.prob

## take max value across both visits (using TOT_det)
BBS_final_detect2 <- BBS_final_detect2 %>% group_by(YEAR, GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(TOT_det = max(TOT_det))
BBS_final_detect2 <- as.data.frame(BBS_final_detect2)
head(BBS_final_detect2)
## remove NaN (species without detection which get removed later anyway)
BBS_final_detect2 <- BBS_final_detect2[complete.cases(BBS_final_detect2), ]
length(unique(BBS_final_detect2$ENGLISH_NAME)) ## 74 species
length(unique(BBS_final_detect2$GRIDREF)) ## 200 sites

## only include species with at least 12 years of data (non-consecutive)
BBS_final2 <- BBS_final_detect2 %>% group_by(GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(no_years = n())
## remove site-species combos with less than 12 years of data
BBS_final2 <- BBS_final2[BBS_final2$no_years>=12,]
length(unique(BBS_final2$ENGLISH_NAME)) ## 74 species
length(unique(BBS_final2$GRIDREF)) ## 200 sites

## merge species and sites into BBS_final
BBS_final3 <- merge(BBS_final2, BBS_final_detect2, by=c("ENGLISH_NAME","CBC_CODE", "GRIDREF"), all=FALSE)
length(unique(BBS_final3$ENGLISH_NAME)) ## 74 species
length(unique(BBS_final3$GRIDREF)) ## 200 sites
## 61,234 rows
##########################################################################


######################################### FILL IN MISSING ABUNDANCE DATA ###############################################

## remove no.years column
BBS_final4 <- BBS_final3[,-4]
BBS_final4$YEAR <- as.factor(BBS_final4$YEAR)
BBS_final4$ENGLISH_NAME <- as.character(BBS_final4$ENGLISH_NAME)
BBS_final4$CBC_CODE <- as.character(BBS_final4$CBC_CODE)
BBS_final4$GRIDREF <- as.character(BBS_final4$GRIDREF)

## complete data so each species/site combo has 15 years filled with NAs
BBS_final4$YEAR <- as.numeric(levels(BBS_final4$YEAR))[BBS_final4$YEAR]
BBS_final5 <- BBS_final4 %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  complete(ENGLISH_NAME, CBC_CODE, YEAR = 2004:2018)
############ WORKS!! ############ 

## create two new columns: TOT_lag and DIFF
BBS_final5 <- BBS_final5 %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  arrange(YEAR) %>%
  mutate(TOT_lag=lag(TOT_det), DIFF=TOT_det-TOT_lag)

### take random value from DIFF to fill NAs from missing years - run multiple times? e.g. 100 for now
BBS_interpol <- NULL
for(i in 1:100){
  print(i)
  BBS_final6 <- BBS_final5 %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  mutate(DIFF=ifelse(is.na(DIFF) | is.na(TOT_det), sample(na.omit(DIFF)), DIFF))
  results_table_temp <- data.frame(BBS_final6$GRIDREF, BBS_final6$ENGLISH_NAME, BBS_final6$CBC_CODE, BBS_final6$YEAR, BBS_final6$TOT_det,
                                   BBS_final6$TOT_lag, BBS_final6$DIFF,i, stringsAsFactors = FALSE)
  BBS_interpol <- as.data.frame(rbind(BBS_interpol, results_table_temp))
}
colnames(BBS_interpol) <- c("GRIDREF", "ENGLISH_NAME", "CBC_CODE", "YEAR", "TOT_det", "TOT_lag", "DIFF", "i")

## yes/no for missing data
BBS_interpol$missing_data <- ifelse(is.na(BBS_interpol$TOT_det), "yes", "no")
BBS_interpol2 <- BBS_interpol
BBS_interpol2$i <- as.factor(BBS_interpol2$i)

## calculate new TOT using randomly sampled differences
## first descending order to fill in early years 
## run this multiple times
BBS_interpol3 <- BBS_interpol2 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(desc(YEAR)) %>%
  mutate(TOT_det = ifelse(is.na(TOT_det), (lag(TOT_det) - lag(DIFF)), TOT_det))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT_det <- ifelse(BBS_interpol3$TOT_det<0, 0, BBS_interpol3$TOT_det)
## run again
BBS_interpol3 <- BBS_interpol3 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(desc(YEAR)) %>%
  mutate(TOT_det = ifelse(is.na(TOT_det), (lag(TOT_det) - lag(DIFF)), TOT_det))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT_det <- ifelse(BBS_interpol3$TOT_det<0, 0, BBS_interpol3$TOT_det)
## run again
BBS_interpol3 <- BBS_interpol3 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(desc(YEAR)) %>%
  mutate(TOT_det = ifelse(is.na(TOT_det), (lag(TOT_det) - lag(DIFF)), TOT_det))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT_det <- ifelse(BBS_interpol3$TOT_det<0, 0, BBS_interpol3$TOT_det)
## then run in increasing order to fill in late years
## run multiple times
BBS_interpol3 <- BBS_interpol3 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(YEAR) %>%
  mutate(TOT_det = ifelse(is.na(TOT_det), (lag(TOT_det) + DIFF), TOT_det))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT_det <- ifelse(BBS_interpol3$TOT_det<0, 0, BBS_interpol3$TOT_det)
## run again
BBS_interpol3 <- BBS_interpol3 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(YEAR) %>%
  mutate(TOT_det = ifelse(is.na(TOT_det), (lag(TOT_det) + DIFF), TOT_det))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT_det <- ifelse(BBS_interpol3$TOT_det<0, 0, BBS_interpol3$TOT_det)
## run again
BBS_interpol3 <- BBS_interpol3 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(YEAR) %>%
  mutate(TOT_det = ifelse(is.na(TOT_det), (lag(TOT_det) + DIFF), TOT_det))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT_det <- ifelse(BBS_interpol3$TOT_det<0, 0, BBS_interpol3$TOT_det)

sum(is.na(BBS_interpol3$TOT_det)) ## zero NAs 

## save dataframe
write.csv(BBS_interpol3, file="../Data/BBS_abund_data/BBS_2004_2018_interpol_repeats_det2.csv", row.names=FALSE) ## 180 sites and 70 species from 2004-2018

rm(list=ls()) # clear R
seed_spp_list <- read.csv("../Data/Trait data/seed_dispersing_spp_list.csv", header=TRUE)
invert_spp_list <- read.csv("../Data/Trait data/insectivores_spp_list.csv", header=TRUE)
BBS_final <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_interpol_repeats_det2.csv", header=TRUE)

## merge seed dispersing species list (so we only have BBS data for seed dispersers)
seed_spp_list <- seed_spp_list[,-c(2,3,5)]
seed_spp_list <- seed_spp_list[seed_spp_list$yes_no==1,] ## 53 species which are seed dispersers
BBS_seed <- merge(BBS_final, seed_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_seed$ENGLISH_NAME)) ## 39 species which have BBS data
length(unique(BBS_seed$GRIDREF)) ## 200 sites

## merge invertivore species list (so we only have BBS data for pest control species)
invert_spp_list <- invert_spp_list[,-c(2,4)]
invert_spp_list <- invert_spp_list[invert_spp_list$yes_no==1,] ## 100 species which are pest controllers
invert_spp_list$species <- as.character(invert_spp_list$species)
invert_spp_list$species[ invert_spp_list$species == "Cetti’s Warbler" ] <- "Cetti Warbler" ## won't merge with apostrophe
BBS_invert <- merge(BBS_final, invert_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_invert$ENGLISH_NAME)) ## 59 species which have BBS data 
length(unique(BBS_invert$GRIDREF)) ## 200 sites
## 2 species missing all detectability data == chough and rock pipit
## 2 species missing some years of detectability which were filtered out == grey wagtail and whinchat

## save seed dispersal files
write.csv(BBS_seed, file="../Data/BBS_abund_data/BBS_2004_2018_seed_interpol_repeat_det2.csv", row.names=FALSE) ## 180 sites and 39 species from 2004-2018
## save insectivore files
write.csv(BBS_invert, file="../Data/BBS_abund_data/BBS_2004_2018_invert_interpol_repeat_det2.csv", row.names=FALSE) ## 180 sites and 55 species from 2004-2018


