###########################################################
## Title: Calculate mean and stability of seed dispersal function
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: June 2019
##########################################################

rm(list=ls()) # clear R
library(tidyverse)

## add data
BBS_final <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_2004_2018_seed.csv", header=TRUE)
effect <- read.csv("../Data/Trait data/Seed dispersal/seed_effect_traits.csv", header=TRUE) ## 32 spp
response <- read.csv("../Data/Trait data/Seed dispersal/seed_response_traits.csv", header=TRUE) ## 38 spp
all_traits <- read.csv("../Data/Trait data/Seed dispersal/seed_all_traits.csv", header=TRUE) ## 28 spp
## read these files in to get species list

# ## remove unecessary columns from each trait data (only need name and code to merge with BBS data)
# ## only need english name and yes/no
# effect <- effect[,c(1,2)]
# response <- response[,c(1,2)]
# effect_response <- effect_response[,c(1,2)]

## merge each file with BBS_final abundance data
BBS_final_effect <- merge(BBS_final, effect, by=c("ENGLISH_NAME", "CBC_CODE"), all=FALSE)
length(unique(BBS_final_effect$ENGLISH_NAME)) ## 32 species
BBS_final_response <- merge(BBS_final, response, by=c("ENGLISH_NAME", "CBC_CODE"), all=FALSE)
length(unique(BBS_final_response$ENGLISH_NAME)) ## 38 species
BBS_final_all <- merge(BBS_final, all_traits, by=c("ENGLISH_NAME", "CBC_CODE"), all=FALSE)
length(unique(BBS_final_all$ENGLISH_NAME)) ## 28 species

############################################################################################################
############################################################################################################
#### MEAN FUNCTION

# BBS_mean_effect <- BBS_final_effect %>% group_by(GRIDREF) %>% summarise(Mean_abund = mean(TOT)) ## 32 species at 199 sites
# BBS_mean_response <- BBS_final_response %>% group_by(GRIDREF) %>% summarise(Mean_abund = mean(TOT)) ## 38 species at 200 sites
# BBS_mean_effect_response <- BBS_final_effect_response %>% group_by(GRIDREF) %>% summarise(Mean_abund = mean(TOT)) ## 28 species at 199 sites

## wrong way to calculate mean abundance --> should take total across all years for each species, then mean across all species at each site
## this is done below anyway to calculate stability

############################################################################################################
############################################################################################################
#### STABILITY OF FUNCTION
## calculate mean and SD for each species and site

### effect species (=32)
## calc total abundance of each species at each site
BBS_tot_effect <- BBS_final_effect %>% group_by(GRIDREF,ENGLISH_NAME) %>% summarise(total = sum(TOT))
## calc mean and SD at each site
BBS_proxy_effect <- BBS_tot_effect %>% group_by(GRIDREF) %>% summarise(mean = mean(total), sd=sd(total))
BBS_proxy_effect$CV <- BBS_proxy_effect$sd/BBS_proxy_effect$mean ## calculate CV (=SD/mean)
BBS_proxy_effect$stability <- 1/BBS_proxy_effect$CV ## calculate stability (1/CV)
## remove NA values (where there is only one species at a site - can't calc SD or stability)
BBS_proxy_effect <- na.omit(BBS_proxy_effect) ## 199 sites

## response species (=38)
BBS_tot_response <- BBS_final_response %>% group_by(GRIDREF,ENGLISH_NAME) %>% summarise(total = sum(TOT))
## calc mean and SD at each site
BBS_proxy_response <- BBS_tot_response %>% group_by(GRIDREF) %>% summarise(mean = mean(total), sd=sd(total))
BBS_proxy_response$CV <- BBS_proxy_response$sd/BBS_proxy_response$mean ## calculate CV (=SD/mean)
BBS_proxy_response$stability <- 1/BBS_proxy_response$CV ## calculate stability (1/CV)
## remove NA values (where there is only one species at a site - can't calc SD or stability)
BBS_proxy_response <- na.omit(BBS_proxy_response) ## 199 sites

## all trait species (=28)
BBS_tot_all <- BBS_final_all %>% group_by(GRIDREF, ENGLISH_NAME) %>% summarise(total = sum(TOT))
## calc mean and SD at each site
BBS_proxy_all <- BBS_tot_all %>% group_by(GRIDREF) %>% summarise(mean = mean(total), sd=sd(total))
BBS_proxy_all$CV <- BBS_proxy_all$sd/BBS_proxy_all$mean ## calculate CV (=SD/mean)
BBS_proxy_all$stability <- 1/BBS_proxy_all$CV ## calculate stability (1/CV)
## remove NA values (where there is only one species at a site - can't calc SD or stability)
BBS_proxy_all <- na.omit(BBS_proxy_all) ## 199 sites

## remove SD and total columns from each
BBS_proxy_effect <- BBS_proxy_effect[,-c(3:4)]
BBS_proxy_response <- BBS_proxy_response[,-c(3:4)]
BBS_proxy_all <- BBS_proxy_all[,-c(3:4)]

## save files
write.csv(BBS_proxy_effect, file="../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_effect.csv", row.names=FALSE) ## 199 sites
write.csv(BBS_proxy_response, file="../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_response.csv", row.names=FALSE) ## 200 sites
write.csv(BBS_proxy_all, file="../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_all.csv", row.names=FALSE) ## 199 sites


