###########################################################
## Title: Prep BBS trait data for analysis 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2019
##########################################################

rm(list=ls()) # clear R

## load data
spp_list <- read.csv("../Data/BBS_abund_data/BBS_seed_species_list.csv", header=TRUE, fileEncoding="UTF-8-BOM") ## species with BBS data (124 species)
trait_data <- read.csv("../Data/Trait data/UK_bird_traits.csv", header=TRUE) ## species with trait data (600 species, lots of vagrants)
## additional trait data from external sources
bird_STI <- read.csv("../Data/Trait data/Bird_STI_144_spp.csv", header=TRUE)
# clutch_size <- read.csv("../Data/Trait data/Clutchsizdata_Jetz2008.csv", header=TRUE)
life_therm_hsi <- read.csv("../Data/Trait data/lifespan_thermalmax_hsi_traits.csv", header=TRUE)
brood_ssi <- read.csv("../Data/Trait data/Johnstone_etal_2014_traits.csv", header=TRUE)
clutch_brood_life <- read.csv("../Data/Trait data/Amniote_clutch_brood_longevity.csv", header=TRUE)

###################################################################################
###################################################################################

## match up trait data and BBS abundance data so we have species with both (114 species)

## change common name column name (reads in weirdly)
names(trait_data)[1]<-paste("Common_name")
## change column types to characters
spp_list$ENGLISH_NAME <- as.character(spp_list$ENGLISH_NAME)
trait_data$Common_name  <- as.character(trait_data$Common_name )

## remove apostrophes from trait_data file (can't remove them using gsub)
trait_data$Common_name[ trait_data$Common_name == "Cettiâ€™s Warbler" ] <- "Cetti's Warbler"
trait_data$Common_name[ trait_data$Common_name == "Montaguâ€™s Harrier" ] <- "Montagu's Harrier"
trait_data$Common_name[ trait_data$Common_name == "Saviâ€™s Warbler" ] <- "Savi's Warbler"

## merge trait data with BBS species list
bbs_trait_data <- merge(spp_list, trait_data, by.x="ENGLISH_NAME", by.y="Common_name", all=FALSE) ## 43 species (all seed-eating birds with BBS and trait data)

## complete the BBS_trait_data by adding more traits from external sources
## remove apostrophes from johnstone file (can't remove them using gsub)
names(brood_ssi)[1]<-paste("Species")
brood_ssi$Species  <- as.character(brood_ssi$Species )
brood_ssi$Species[ brood_ssi$Species == "Cettiâ€™s Warbler" ] <- "Cetti's Warbler"

## remove uneccessary columns
# clutch_size <- clutch_size[,c(3,4)] ## scientific name and clutch size
brood_ssi <- brood_ssi[,c(1,6)] ## common name and species specialisation index

## merge trait data with johnstone data (brood number and ssi)
bbs_trait_data <- merge(bbs_trait_data, brood_ssi, by.x="ENGLISH_NAME", by.y="Species", all.x=TRUE)
## merge trait data with clutch_brood_life
bbs_trait_data <- merge(bbs_trait_data, clutch_brood_life, by.x=c("ENGLISH_NAME", "Scientific_name"),
                        by.y=c("ENGLISH_NAME", "Scientific_name"), all.x=TRUE)
## merge trait data with STI 
bbs_trait_data <- merge(bbs_trait_data, bird_STI, by.x="ENGLISH_NAME", by.y="sp", all.x=TRUE)
#### bbs_trait_data has 43 species with 81 columns

## not all these columns are needed
bbs_trait_data <- bbs_trait_data[,-c(4:15,21,23:27,34,39:50,61:76)]
## 43 species with 34 columns
summary(bbs_trait_data)

## save file
write.csv(bbs_trait_data, file="../Data/Analysis_data/BBS_seed_trait_data.csv", row.names=FALSE)
## trait data for 114 species
## seeed-eating species which have trait AND abundance data (43 species)
bbs_trait_data <- read.csv("../Data/Analysis_data/BBS_seed_trait_data.csv", header=TRUE)

#########
## Effect traits
  # Bill length/width/depth
  # Gape width

effect_traits <- bbs_trait_data[,c(1,3,9,10,12:13,19)]
effect_traits <- na.omit(effect_traits)
## 32 species with effect traits (plus body mass)

#########
## Response traits
  # Species specialisation index
  # Species temperature index
  # Mean latitude  
  # Lifespan
  # Clutch size 
  # Number of broods

response_traits <- bbs_trait_data[,c(1,3,6,30:34)]
response_traits <- na.omit(response_traits)
## 38 species with complete response trait data

######
## Both traits
  # Wing length
  # Tail length
  # Kipp's distance
  # Tarsus length

effect_both_traits <- bbs_trait_data[,c(1,3,9,10,12:13,19,14:16,18)]
effect_both_traits <- na.omit(effect_both_traits)
## 32 species with complete effect and both trait data

response_both_traits <- bbs_trait_data[,c(1,3,6,30:34,14:16,18)]
response_both_traits <- na.omit(response_both_traits)
## 38 species with complete response and both trait data

##### merge effect and response to find subset of species which have ALL analysed traits (subset to analyse)
all_traits <- merge(effect_both_traits, response_both_traits, by=c("ENGLISH_NAME", "CBC_CODE"), all=FALSE)
## 28 seed-eating species which have ALL traits (and abundance data)

## save these files
write.csv(effect_traits, file="../Data/Trait data/Seed dispersal/seed_effect_traits.csv", row.names=FALSE)
write.csv(response_traits, file="../Data/Trait data/Seed dispersal/seed_response_traits.csv", row.names=FALSE)
write.csv(effect_both_traits, file="../Data/Trait data/Seed dispersal/seed_effect_both_traits.csv", row.names=FALSE)
write.csv(response_both_traits, file="../Data/Trait data/Seed dispersal/seed_response_both_traits.csv", row.names=FALSE)
write.csv(all_traits, file="../Data/Trait data/Seed dispersal/seed_all_traits.csv", row.names=FALSE)

##################################################################################################################
##################################################################################################################

# ###### merge with BBS abundance data to match number of species (was 124, now needs to be 114)
# spp_list <- bbs_trait_data[,c(1)]
# spp_list <- as.data.frame(spp_list)
# spp_list <- unique(spp_list) ## 114
# ## read in BBS 2004-2018 data
# BBS_final <- read.csv("../Data/BBS_abund_data/BBS_2004_2018.csv", header=TRUE) ## 124 species (needs to be 114)
# BBS_final$ENGLISH_NAME <- as.character(BBS_final$ENGLISH_NAME)
# BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Pied/White Wagtail" ] <- "Pied Wagtail"
# BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Common Crossbill" ] <- "Crossbill"
# 
# ## merge together
# BBS_final <- merge(BBS_final, spp_list, by.x="ENGLISH_NAME", by.y="spp_list", all=FALSE)
# length(unique(BBS_final$ENGLISH_NAME)) ## 114 species
# length(unique(BBS_final$GRIDREF)) ## 200 sites
# write.csv(BBS_final, file="../Data/Analysis_data/BBS_2004_2018_final.csv", row.names=FALSE) ## 200 sites and 114 species from 2004-2018
# ## species which have both abundance data and trait data (used to calc mean and stability of function)
