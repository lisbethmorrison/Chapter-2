###########################################################
## Title: Prep BBS trait data for analysis 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2019
##########################################################

rm(list=ls()) # clear R

## load data
seed_spp_list <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_seed_species_list.csv", header=TRUE, fileEncoding="UTF-8-BOM") ## species with BBS data (43 species)
invert_spp_list <- read.csv("../Data/Analysis_data/Pest control/BBS_invert_species_list.csv", header=TRUE, fileEncoding="UTF-8-BOM") ## species with BBS data (74 species)
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
seed_spp_list$ENGLISH_NAME <- as.character(seed_spp_list$ENGLISH_NAME)
invert_spp_list$ENGLISH_NAME <- as.character(invert_spp_list$ENGLISH_NAME)
trait_data$Common_name  <- as.character(trait_data$Common_name )

## remove apostrophes from trait_data file (can't remove them using gsub)
trait_data$Common_name[ trait_data$Common_name == "Cettiâ€™s Warbler" ] <- "Cetti Warbler"
trait_data$Common_name[ trait_data$Common_name == "Montaguâ€™s Harrier" ] <- "Montagu's Harrier"
trait_data$Common_name[ trait_data$Common_name == "Saviâ€™s Warbler" ] <- "Savi's Warbler"

#######################################################
################### SEED DISPERSAL ####################
#######################################################

## merge trait data with BBS species list
seed_bbs_trait_data <- merge(seed_spp_list, trait_data, by.x="ENGLISH_NAME", by.y="Common_name", all=FALSE) ## 40 species (all seed-eating birds with BBS and trait data)

## complete the BBS_trait_data by adding more traits from external sources
## remove apostrophes from johnstone file (can't remove them using gsub)
names(brood_ssi)[1]<-paste("Species")
brood_ssi$Species  <- as.character(brood_ssi$Species )
brood_ssi$Species[ brood_ssi$Species == "Cettiâ€™s Warbler" ] <- "Cetti Warbler"

## remove uneccessary columns
# clutch_size <- clutch_size[,c(3,4)] ## scientific name and clutch size
brood_ssi <- brood_ssi[,c(1,6)] ## common name and species specialisation index

## merge trait data with johnstone data (brood number and ssi)
seed_bbs_trait_data <- merge(seed_bbs_trait_data, brood_ssi, by.x="ENGLISH_NAME", by.y="Species", all.x=TRUE)
## merge trait data with clutch_brood_life
seed_bbs_trait_data <- merge(seed_bbs_trait_data, clutch_brood_life, by.x=c("ENGLISH_NAME", "Scientific_name"),
                        by.y=c("ENGLISH_NAME", "Scientific_name"), all.x=TRUE)
## merge trait data with STI 
seed_bbs_trait_data <- merge(seed_bbs_trait_data, bird_STI, by.x="ENGLISH_NAME", by.y="sp", all.x=TRUE)
## merge trait data with life_therm_hsi
seed_bbs_trait_data <- merge(seed_bbs_trait_data, life_therm_hsi, by.x="ENGLISH_NAME", by.y="common_name", all.x=TRUE)
#### seed_bbs_trait_data has 40 species with 84 columns

## not all these columns are needed
seed_bbs_trait_data <- seed_bbs_trait_data[,-c(4:17,19,21,23:27,34,39:76)]
## 40 species with 24 columns
summary(seed_bbs_trait_data)

## save file
write.csv(seed_bbs_trait_data, file="../Data/Analysis_data/Seed dispersal/BBS_seed_trait_data.csv", row.names=FALSE)
## seeed-eating species which have trait AND abundance data (40 species)
seed_bbs_trait_data <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_seed_trait_data.csv", header=TRUE)

#########
## Effect traits
  # Bill length/width/depth
  # Gape width

effect_traits <- seed_bbs_trait_data[,c(1,3,9,10,12:13,19)]
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

response_traits <- seed_bbs_trait_data[,c(1,3,6,30:34)]
response_traits <- na.omit(response_traits)
## 38 species with complete response trait data

######
## Both traits
  # Wing length
  # Tail length
  # Kipp's distance
  # Tarsus length

effect_both_traits <- seed_bbs_trait_data[,c(1,3,9,10,12:13,19,14:16,18)]
effect_both_traits <- na.omit(effect_both_traits)
## 32 species with complete effect and both trait data

response_both_traits <- seed_bbs_trait_data[,c(1,3,6,30:34,14:16,18)]
response_both_traits <- na.omit(response_both_traits)
## 38 species with complete response and both trait data

##### merge effect and response to find subset of species which have ALL analysed traits (subset to analyse)
all_traits <- merge(effect_both_traits, response_both_traits, by=c("ENGLISH_NAME", "CBC_CODE", "Tarsus_Length", 
                  "Kipp.s_Distance", "Wing_Chord", "Tail_Length"), all=FALSE)
## 28 seed-eating species which have ALL traits (and abundance data)

## save these files
write.csv(effect_traits, file="../Data/Trait data/Seed dispersal/seed_effect_traits.csv", row.names=FALSE) # 32
write.csv(response_traits, file="../Data/Trait data/Seed dispersal/seed_response_traits.csv", row.names=FALSE) # 38
write.csv(effect_both_traits, file="../Data/Trait data/Seed dispersal/seed_effect_both_traits.csv", row.names=FALSE) # 32
write.csv(response_both_traits, file="../Data/Trait data/Seed dispersal/seed_response_both_traits.csv", row.names=FALSE) # 38
write.csv(all_traits, file="../Data/Trait data/Seed dispersal/seed_all_traits.csv", row.names=FALSE) # 28

#####################################################
################### PEST CONTROL ####################
#####################################################

## merge trait data with BBS species list
invert_bbs_trait_data <- merge(invert_spp_list, trait_data, by.x="ENGLISH_NAME", by.y="Common_name", all=FALSE) ## 74 species (all invertivores with BBS and trait data)

## complete the BBS_trait_data by adding more traits from external sources
## merge trait data with johnstone data (brood number and ssi)
invert_bbs_trait_data <- merge(invert_bbs_trait_data, brood_ssi, by.x="ENGLISH_NAME", by.y="Species", all.x=TRUE)
## merge trait data with clutch_brood_life
invert_bbs_trait_data <- merge(invert_bbs_trait_data, clutch_brood_life, by.x=c("ENGLISH_NAME", "Scientific_name"),
                             by.y=c("ENGLISH_NAME", "Scientific_name"), all.x=TRUE)
## merge trait data with STI 
invert_bbs_trait_data <- merge(invert_bbs_trait_data, bird_STI, by.x="ENGLISH_NAME", by.y="sp", all.x=TRUE)

## merge trait data with life_therm_hsi
invert_bbs_trait_data <- merge(invert_bbs_trait_data, life_therm_hsi, by.x="ENGLISH_NAME", by.y="common_name", all.x=TRUE)

#### invert_bbs_trait_data has 66 species with 84 columns

## not all these columns are needed
invert_bbs_trait_data <- invert_bbs_trait_data[,-c(4:17,19,21,23:27,34,39:76)]
## 66 species with 24 columns
summary(invert_bbs_trait_data)

## save file
write.csv(invert_bbs_trait_data, file="../Data/Analysis_data/Pest control/BBS_invert_trait_data.csv", row.names=FALSE)
## invert-eating species which have trait AND abundance data (66 species)
## also merge both seed and invert trait databases (complete data for ALL species analysed)
# seed_bbs_trait_data$eco_function <- "seed"
# invert_bbs_trait_data$eco_function <- "invert"
trait_data <- rbind(seed_bbs_trait_data, invert_bbs_trait_data)
trait_data <- unique(trait_data) ## 81 unique species (25 species which eat seeds and inverts)
trait_data_info  <- colSums(!is.na(trait_data))
trait_data_info <- as.data.frame(trait_data_info)
trait_data_info$trait <- row.names(trait_data_info)
rownames(trait_data_info) <- 1:nrow(trait_data_info)
## 9 out of 24 traits are missing data (two traits are duplicates - lifespan and clutch size)
## save this file
write.csv(trait_data, file="../Data/Trait data/ALL_trait_data.csv", row.names=FALSE)

invert_bbs_trait_data <- read.csv("../Data/Analysis_data/Pest control/BBS_invert_trait_data.csv", header=TRUE)

#########
## Effect traits
# Bill length/width/depth
# Gape width

effect_traits <- invert_bbs_trait_data[,c(1,3,9,10,12:13,19)]
effect_traits <- na.omit(effect_traits)
## 48 species with effect traits (plus body mass)

#########
## Response traits
# Species specialisation index
# Species temperature index
# Mean latitude  
# Lifespan
# Clutch size 
# Number of broods

response_traits <- invert_bbs_trait_data[,c(1,3,6,30:34)]
response_traits <- na.omit(response_traits)
## 63 species with complete response trait data

######
## Both traits
# Wing length
# Tail length
# Kipp's distance
# Tarsus length

effect_both_traits <- invert_bbs_trait_data[,c(1,3,9,10,12:13,19,14:16,18)]
effect_both_traits <- na.omit(effect_both_traits)
## 48 species with complete effect and both trait data

response_both_traits <- invert_bbs_trait_data[,c(1,3,6,30:34,14:16,18)]
response_both_traits <- na.omit(response_both_traits)
## 63 species with complete response and both trait data

##### merge effect and response to find subset of species which have ALL analysed traits (subset to analyse)
all_traits <- merge(effect_both_traits, response_both_traits, by=c("ENGLISH_NAME", "CBC_CODE", "Tarsus_Length", 
                  "Kipp.s_Distance", "Wing_Chord", "Tail_Length"), all=FALSE)
## 42 seed-eating species which have ALL traits (and abundance data)

## save these files
write.csv(effect_traits, file="../Data/Trait data/Pest control/invert_effect_traits.csv", row.names=FALSE) # 48
write.csv(response_traits, file="../Data/Trait data/Pest control/invert_response_traits.csv", row.names=FALSE) # 63
write.csv(effect_both_traits, file="../Data/Trait data/Pest control/invert_effect_both_traits.csv", row.names=FALSE) # 48
write.csv(response_both_traits, file="../Data/Trait data/Pest control/invert_response_both_traits.csv", row.names=FALSE) # 63
write.csv(all_traits, file="../Data/Trait data/Pest control/invert_all_traits.csv", row.names=FALSE) # 42

