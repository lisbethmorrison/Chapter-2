###########################################################
## Title: Calculate mean and stability of pest control function
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: August 2019
##########################################################

rm(list=ls()) # clear R
library(tidyverse)

## add data
BBS_final <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_invert_interpol_repeat_det2.csv", header=TRUE)

length(unique(BBS_final$ENGLISH_NAME)) ## 59 species (55 for det1)
length(unique(BBS_final$GRIDREF)) ## 200 sites (180 for det1)
length(unique(BBS_final$YEAR)) ## 15 years

## create species richness data frame for each site
## number of rows (i.e. number of species) for each site, year and iteration
spp_rich_invert <- BBS_final %>% group_by(GRIDREF,YEAR,i) %>% summarise(n=n())
## take average spp richness over iterations
spp_rich_invert <- spp_rich_invert %>% group_by(GRIDREF,YEAR) %>% summarise(average=mean(n))
## take average spp richness over time
spp_rich_invert <- spp_rich_invert %>% group_by(GRIDREF) %>% summarise(average_spp_richness=mean(average))
## save file
write.csv(spp_rich_invert, file="../Data/Analysis_data/Pest Control/mean_spp_richness.csv", row.names=FALSE) ## 200 sites


# ############################################################################################################
# ############################################################################################################
# #### STABILITY AND MEAN OF FUNCTION
# ## calculate mean and SD for each site, and use this to calculate stability

## calculate this for all 40 species which are seed-eaters for each reptition of abundance interpolation (100 repeats)
BBS_tot <- BBS_final %>% group_by(GRIDREF,YEAR,i) %>% summarise(tot_abund = sum(TOT_det))
## calculate the mean and SD over time for each site
BBS_tot2 <- BBS_tot %>% group_by(GRIDREF, i) %>% summarise(mean_abund=mean(tot_abund), sd=sd(tot_abund))
## calculate CV 
BBS_tot2$CV <- BBS_tot2$sd/BBS_tot2$mean_abund
## caculate stability
BBS_tot2$stability <- 1/BBS_tot2$CV
## take mean stability and mean mean_abund for each site across repetitions
BBS_proxy <- BBS_tot2 %>% group_by(GRIDREF) %>% summarise(stability=mean(stability), mean_abund=mean(mean_abund))

## save file
write.csv(BBS_proxy, file="../Data/Analysis_data/Pest control/BBS_proxy_invert_det2.csv", row.names=FALSE) ## 200 sites
## 55 species, 180 sites for det1
