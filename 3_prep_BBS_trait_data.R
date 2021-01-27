###########################################################
## Title: Prep BBS trait data for analysis 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2019
##########################################################

rm(list=ls()) # clear R

library(dplyr)

## load data
seed_spp_list <- read.csv("../Data/Trait data/Seed dispersal/BBS_seed_species_list.csv", header=TRUE, fileEncoding="UTF-8-BOM") ## species with BBS data (40 species)
invert_spp_list <- read.csv("../Data/Trait data/Pest control/BBS_invert_species_list.csv", header=TRUE, fileEncoding="UTF-8-BOM") ## species with BBS data (67 species)
trait_data <- read.csv("../Data/Trait data/UK_bird_traits.csv", header=TRUE) ## species with trait data (600 species, lots of vagrants)
## additional trait data from external sources
bird_STI <- read.csv("../Data/Trait data/Devictor2012_STI.csv", header=TRUE)
clutch_size <- read.csv("../Data/Trait data/Jetz2008_clutchsize.csv", header=TRUE)
therm_max <- read.csv("../Data/Trait data/Jiguet_thermalmax.csv", header=TRUE)
brood_ssi <- read.csv("../Data/Trait data/Johnstone2014_brood_ssi.csv", header=TRUE)
clutch_brood_long <- read.csv("../Data/Trait data/Myhrvold2015_clutch_brood_maxlongevity.csv", header=TRUE)
gape_width <- read.csv("../Data/Trait data/Extra_gape_width_trait.csv", header=TRUE)
surrogate_species <- read.csv("../Data/Trait data/Surrogate_trait_info.csv", header=TRUE)

###################################################################################
###################################################################################

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

## merge trait data with extra gape width data
trait_data <- merge(trait_data, gape_width, by.x="Common_name", by.y="ENGLISH_NAME", all=TRUE)
## mege the two gape width columns together
trait_data$Gape_width <- coalesce(trait_data$Gape_width, trait_data$gape_width)
## remove extra columns
trait_data <- trait_data[,-c(76:77)]

## complete the BBS_trait_data by adding more traits from external sources
## edit species names to match (cetti's warbler and pied wagtail)
names(brood_ssi)[1]<-paste("Species")
brood_ssi$Species  <- as.character(brood_ssi$Species )
brood_ssi$Species[ brood_ssi$Species == "Cettiâ€™s Warbler" ] <- "Cetti Warbler"
clutch_brood_long$common_name <- gsub("White Wagtail", "Pied Wagtail", clutch_brood_long$common_name)
## merge genus and species name to create scientific name to merge with
clutch_brood_long$scientific_name <- paste(clutch_brood_long$ï..genus, clutch_brood_long$species, sep=" ")
bird_STI$sp <- gsub("Cetti's Warbler", "Cetti Warbler", bird_STI$sp )
bird_STI$sp <- gsub("Carrion/Hooded Crow", "Carrion Crow", bird_STI$sp )
bird_STI$sp <- gsub("Rock Dove/Feral Pigeon", "Rock Dove", bird_STI$sp )
## remove uneccessary columns
brood_ssi <- brood_ssi[,c(1,2,5)] ## common name and species specialisation index
names(brood_ssi)[2] <- paste("Johnstone_broodsize")
names(clutch_size)[4] <- paste("Jetz_clutchsize")

## merge trait data with johnstone data (brood number and ssi)
trait_data <- merge(trait_data, brood_ssi, by.x="Common_name", by.y="Species", all.x=TRUE)
## merge trait data with clutch_brood_long
trait_data <- merge(trait_data, clutch_brood_long, by.x="Unique_Scientific_Name",
                        by.y="scientific_name", all.x=TRUE)
trait_data <- merge(trait_data, clutch_brood_long, by.x="Scientific_name",
                             by.y="scientific_name", all.x=TRUE)
trait_data$litter_or_clutch_size_n <- coalesce(trait_data$litter_or_clutch_size_n.x, trait_data$litter_or_clutch_size_n.y)
trait_data$litters_or_clutches_per_y <- coalesce(trait_data$litters_or_clutches_per_y.x, trait_data$litters_or_clutches_per_y.y)
trait_data$maximum_longevity_y <- coalesce(trait_data$maximum_longevity_y.x, trait_data$maximum_longevity_y.y)

## merge trait data with STI 
trait_data <- merge(trait_data, bird_STI, by.x="Common_name", by.y="sp", all.x=TRUE)
## merge trait data with therm_max
trait_data <- merge(trait_data, therm_max, by.x=c("Common_name", "Scientific_name"), by.y=c("ENGLISH_NAME", "Scientific_name"), all.x=TRUE)
## merge trait data with clutch_size
trait_data <- merge(trait_data, clutch_size, by.x="Scientific_name", by.y="Species", all.x=TRUE)

## not all these columns are needed
trait_data <- trait_data[,c(1:3,17,19,21,27,29:32,34:37,76:77,90:93,95,98)]
## 40 species with 24 columns
summary(trait_data)

## merge the two brood_size datasets to fill NAs
trait_data$brood_size <- coalesce(trait_data$Johnstone_broodsize, trait_data$litters_or_clutches_per_y)
trait_data <- trait_data[,-c(16,19)] ## remove extra brood size columns
trait_data <- trait_data[,-c(5,21)] ## remove extra clutch size columns 

###################################################################################################################################################
###################################################################################################################################################

## merge trait data with BBS SEED DISPERSAL species list
seed_bbs_trait_data <- merge(seed_spp_list, trait_data, by.x="ENGLISH_NAME", by.y="Common_name", all=FALSE) ## 40 species (all seed-eating birds with BBS and trait data)
## use Sardinian Warbler longevity for Dartford Warbler
longevity <- 8.333333333
seed_bbs_trait_data[,18][is.na(seed_bbs_trait_data[,18])] <- longevity ## fill in Dartford Warbler NA using Sardinian Warbler

## merge trait data with BBS PEST CONTROL species list
invert_bbs_trait_data <- merge(invert_spp_list, trait_data, by.x="ENGLISH_NAME", by.y="Common_name", all=FALSE) ## 67 species (all invert-eating birds with BBS and trait data)
invert_bbs_trait_data[,18][is.na(invert_bbs_trait_data[,18])] <- longevity ## fill in Dartford Warbler NA using Sardinian Warbler
## use gape width for Western reef-Egret
gape_width <- 14.4
invert_bbs_trait_data[,15][is.na(invert_bbs_trait_data[,15])] <- gape_width ## fill in Little Egret NA using Western reef-Egret

## save files
write.csv(seed_bbs_trait_data, file="../Data/Trait data/Seed dispersal/BBS_seed_trait_data.csv", row.names=FALSE)
## seed-eating species which have trait AND abundance data (40 species)
write.csv(invert_bbs_trait_data, file="../Data/Trait data/Pest control/BBS_invert_trait_data.csv", row.names=FALSE)
## invert-eating species which have trait AND abundance data (59 species)

## compile both together with traits for all 81 unique species
trait_data <- rbind(seed_bbs_trait_data, invert_bbs_trait_data)
trait_data <- unique(trait_data) ## 75 unique species (becomes 74 when crane is removed in seed_FD script)
trait_data_info  <- colSums(is.na(trait_data))
trait_data_info <- as.data.frame(trait_data_info)
trait_data_info$trait <- row.names(trait_data_info)
rownames(trait_data_info) <- 1:nrow(trait_data_info)
colnames(trait_data_info)[1] <- "number_missing"
trait_data_info$percent_missing <- (trait_data_info$number_missing/81)*100 ## calculate percentage of missing data
trait_data_info <- trait_data_info[c(2,1,3)] ## re-order columns
trait_data_info <- trait_data_info[-c(1:4), ] ## remove first 4 rows - english name, CBC code and scientific name 
## 17 traits
## 2 out of 17 traits are missing data
## merge in surrogate species info
trait_data <- merge(trait_data, surrogate_species, by.x=c("ENGLISH_NAME", "Scientific_name", by.y=c("ENGLISH_NAME", "Scientific_name")))
## save files
write.csv(trait_data, file="../Data/Trait data/ALL_trait_data.csv", row.names=FALSE)
write.csv(trait_data_info, file="../Data/Trait data/ALL_trait_data_info.csv", row.names=FALSE)

# #########
# ## Effect traits
#   # Bill length/width/depth
#   # Gape width
# 
# #########
# ## Response traits
#   # Species specialisation index
#   # Species temperature index
#   # Mean latitude  
#   # Lifespan
#   # Clutch size 
#   # Number of broods
#   # Thermal maximum
# 
# ######
# ## Both traits
#   # Wing length
#   # Tail length
#   # Kipp's distance
#   # Tarsus length
#   # Body size
#   # HWI


