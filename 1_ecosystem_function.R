###########################################################
## Title: Calculate mean and stability of function
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: June 2019
##########################################################

rm(list=ls()) # clear R

## add data
BBS_final <- read.csv("../Data/BBS_2004_2018.csv", header=TRUE)
BBS_traits <- read.csv("../Data/BBS_species_trait_data.csv", header=TRUE)

## remove unecessary columns from trait data (only interested in yes/no for seed dispersal)
## only need english name and yes/no
BBS_traits <- BBS_traits[,c(1,13)]

## merge this in with BBS_final abundance data
BBS_final <- merge(BBS_traits, BBS_final, by=c("ENGLISH_NAME"), all=FALSE)
length(unique(BBS_final$ENGLISH_NAME)) ## 112 species (12 species don't have trait data)

## remove species which do not disperse seeds
BBS_final_seed <- BBS_final[BBS_final$yes.no.seed.dispersal=="yes",]
length(unique(BBS_final_seed$ENGLISH_NAME)) ## 73 species eat seeds

#### MEAN FUNCTION
BBS_mean <- BBS_final_seed %>% group_by(GRIDREF) %>% summarise(Mean_abund = mean(TOT))
## mean abundance of seed-eating birds at each site (200 sites)

#### STABILITY OF FUNCTION
## calculate mean and SD for each species and site
BBS_summ <- BBS_final_seed %>% group_by(GRIDREF,ENGLISH_NAME) %>% summarise(mean = mean(TOT), sd = sd(TOT))
BBS_summ$CV <- BBS_summ$sd/BBS_summ$mean ## calculate SD (=SD/mean)
BBS_summ$stability <- 1/BBS_summ$CV ## calculate stability (1/CV)
## remove inf values here? 10 species 
## these are species which have no variation in abundance over time (so CV=0, and stability cannot be calculated)
BBS_summ[sapply(BBS_summ, is.infinite)] <- NA
BBS_summ <- na.omit(BBS_summ)

## calculate mean stability across species for each site
BBS_stability <- BBS_summ %>% group_by(GRIDREF) %>% summarise(stability = mean(stability))

## cbind BBS_mean and BBS_stability
BBS_proxy <- cbind(BBS_mean, BBS_stability)
## remove duplicated gridref column
BBS_proxy[3] <- NULL  

## save file
write.csv(BBS_proxy, file="../Data/Seed dispersal/BBS_proxy_mean_stability.csv", row.names=FALSE)


