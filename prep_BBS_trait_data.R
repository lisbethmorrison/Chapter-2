###########################################################
## Title: Prep BBS trait data for analysis 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2019
##########################################################

rm(list=ls()) # clear R

## load data
spp_list <- read.csv("../Data/BBS_species_list.csv", header=TRUE, fileEncoding="UTF-8-BOM") ## this is the list of species names in BBS data (198 species)
trait_data <- read.csv("../Data/Trait data/UK_bird_traits.csv", header=TRUE) 

names(trait_data)[1]<-paste("Common_name")
spp_list$ENGLISH_NAME <- as.character(spp_list$ENGLISH_NAME)
trait_data$Common_name  <- as.character(trait_data$Common_name )

## 2 species have different names which can be changed
spp_list$ENGLISH_NAME[ spp_list$ENGLISH_NAME == "Pied/White Wagtail" ] <- "Pied Wagtail"
spp_list$ENGLISH_NAME[ spp_list$ENGLISH_NAME == "Common Crossbill" ] <- "Crossbill"

## remove apostrophes from trait_data file (can't remove them using gsub)
trait_data$Common_name[ trait_data$Common_name == "Cettiâ€™s Warbler" ] <- "Cetti's Warbler"
trait_data$Common_name[ trait_data$Common_name == "Montaguâ€™s Harrier" ] <- "Montagu's Harrier"
trait_data$Common_name[ trait_data$Common_name == "Saviâ€™s Warbler" ] <- "Savi's Warbler"

## 11 species do not have trait data
## 4 species need double checking:
## domestic greylag, domestic mallard, feral pigeon and lesser redpoll
## so total should = 194

## now merge these datasets by species name
bbs_trait_data <- merge(spp_list, trait_data, by.x="ENGLISH_NAME", by.y="Common_name", all=FALSE) ## 194 species

write.csv(bbs_trait_data, file="../Data/BBS_species_trait_data.csv", row.names=FALSE)

###### merge with BBS site data to match number of species
spp_list <- bbs_trait_data[,c(1)]
spp_list <- as.data.frame(spp_list)
spp_list <- unique(spp_list) ## 194
## read in BBS 2004-2018 data
BBS_final <- read.csv("../Data/BBS_2004_2018.csv", header=TRUE)
BBS_final$ENGLISH_NAME <- as.character(BBS_final$ENGLISH_NAME)
BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Pied/White Wagtail" ] <- "Pied Wagtail"
BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Common Crossbill" ] <- "Crossbill"

## merge together
BBS_final <- merge(BBS_final, spp_list, by.x="ENGLISH_NAME", by.y="spp_list", all=FALSE)
length(unique(BBS_final$ENGLISH_NAME))
write.csv(BBS_final, file="../Data/BBS_2004_2018.csv", row.names=FALSE) ## 200 sites and 194 species from 2004-2018

############################################################################################################################
### read in trait data again
BBS_traits <- read.csv("../Data/BBS_species_trait_data.csv", header=TRUE)
### manually imported some trait data (BTO lifespan, habitat specialisation index and thermal maximum)

## import other trait data
clutch_size <- read.csv("../Data/Trait data/Clutchsizdata_Jetz2008.csv", header=TRUE)
johnstone_data <- read.csv("../Data/Trait data/Johnstone_etal_2014_traits.csv", header=TRUE)
STI_data <- read.csv("../Data/Trait data/Bird_STI_144_spp.csv", header=TRUE)

## remove apostrophes from johnstone file (can't remove them using gsub)
names(johnstone_data)[1]<-paste("Species")
johnstone_data$Species  <- as.character(johnstone_data$Species )
johnstone_data$Species[ johnstone_data$Species == "Cettiâ€™s Warbler" ] <- "Cetti's Warbler"

## remove uneccessary columns
clutch_size <- clutch_size[,c(5,8)] ## scientific name and clutch size
johnstone_data <- johnstone_data[,c(1,2,6)] ## common name, brood number and species specialisation index

## merge trait data with clutch size data
BBS_traits2 <- merge(BBS_traits, clutch_size, by.x="Scientific_name", by.y="Species", all.x=TRUE)
## merge trait data with johnstone data
BBS_traits3 <- merge(BBS_traits2, johnstone_data, by.x="ENGLISH_NAME", by.y="Species", all.x = TRUE)

sum(is.na(BBS_traits3$Clutch.size))
sum(is.na(BBS_traits3$Clutch)) ## Jetz dataset has less NAs than JT
colSums(is.na(BBS_traits3)) ## Thermal maximum and Habitat specialisation index has most NAs (134)

## merge STI data by common name
BBS_traits4 <- merge(BBS_traits3, STI_data, by.x="ENGLISH_NAME", by.y="sp", all.x=TRUE)

## save file again
write.csv(BBS_traits4, file="../Data/BBS_species_trait_data.csv", row.names=FALSE)
