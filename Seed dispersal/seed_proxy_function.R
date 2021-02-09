###########################################################
## Title: Calculate mean and stability of seed dispersal function
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: June 2019
##########################################################

rm(list=ls()) # clear R
library(tidyverse)

## add data
#BBS_final <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_seed_det_comp.csv", header=TRUE)
BBS_final <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_seed_det_interpol.csv", header=TRUE)

length(unique(BBS_final$ENGLISH_NAME)) ## 45 species // 48 species
length(unique(BBS_final$GRIDREF)) ## 130 sites // 200 sites
length(unique(BBS_final$YEAR)) ## 15 years

## create species richness data frame for each site
## number of rows (i.e. number of species) for each site and year
spp_rich_seed <- BBS_final %>% group_by(GRIDREF,YEAR) %>% summarise(n=n())
## take average spp richness over iterations
spp_rich_seed <- spp_rich_seed %>% group_by(GRIDREF,YEAR) %>% summarise(average=mean(n))
## take average spp richness over time
spp_rich_seed <- spp_rich_seed %>% group_by(GRIDREF) %>% summarise(average_spp_richness=mean(average))
## save file
write.csv(spp_rich_seed, file="../Data/Analysis_data/Seed dispersal/mean_spp_richness_det_interpol.csv", row.names=FALSE) ## 200 sites

############################################################################################################
############################################################################################################
#### STABILITY AND MEAN OF FUNCTION
## calculate total abundance over time at each site
## the calculate mean and SD for each site, use this to calculate stability

## calculate this for all 45 // 48 species which are seed-eaters 
## tot_abund is the total community abundance at each site for each year
BBS_tot <- BBS_final %>% group_by(GRIDREF,YEAR) %>% summarise(tot_abund = sum(TOT_det))
## calculate the mean abund and SD over time for each site
BBS_tot2 <- BBS_tot %>% group_by(GRIDREF) %>% summarise(mean_abund=mean(tot_abund), sd=sd(tot_abund))
## calculate CV 
BBS_tot2$CV <- BBS_tot2$sd/BBS_tot2$mean_abund
## caculate stability
BBS_tot2$stability <- 1/BBS_tot2$CV
## keep mean stability and mean mean_abund for each site
BBS_proxy <- BBS_tot2[,c(1,2,5)] 

## save file
write.csv(BBS_proxy, file="../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_det_interpol.csv", row.names=FALSE) ## 200 sites
## 130 sites, 45 species for complete case detectability
## 200 sites, 48 species for interpolated detectability
