###########################################################
## Title: Calculate functional diversity of traits
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: July 2019
##########################################################

rm(list=ls()) # clear R
options(scipen=999)

library(factoextra)
library(vegan)
library(reshape2)
library(tidyverse)

## read in trait data 
BBS_traits <- read.csv("../Data/BBS_species_trait_data.csv", header=TRUE)

## Effect traits:
  # Bill length/width/depth
  # Tarsus length
  # Gape width

## Response traits
  # Habitat specialisation index XXXXX
  # Species specialisation index
  # Thermal maximum XXXXX
  # Species temperature index
  # Mean latitude  
  # Lifespan
  # Clutch size XXXX
  # Number of broods

## Both traits:
  # Kipp's distance
  # Wing/tail length

####################################################################
####################### EFFECT TRAITS ##############################
####################################################################

## select effect only traits 
effect_traits <- BBS_traits[, c(1,16:20)]
effect_traits[c(2:6)] <- lapply(effect_traits[c(2:6)], function(effect_traits) c(scale(effect_traits, center = TRUE, scale = TRUE))) ## centre and scale trait data
head(effect_traits)
str(effect_traits)
effect_traits <- na.omit(effect_traits) ## 159 species
## correlation matrix 
effect_corr <- cor(effect_traits[c(2:6)])
effect_corr

## PCA of all effect traits
effect <- prcomp(effect_traits[ ,2:6]) ## bill length, width and depth, gape width and tarsus length
summary(effect) # PC1 accounts for 73% of the total variation and PC2 14% of the total variation
## use these two axes to calculatse FD of effect traits

## check eigenvalues to determine the number of principal components to be considered
## An eigenvalue > 1 indicates that PCs account for more variance than accounted by one of the original variables in standardized data
## This is commonly used as a cutoff point for which PCs are retained
## You can also limit the number of component to that number that accounts for a certain fraction of the total variance
## For example, if you are satisfied with 70% of the total variance explained then use the number of components to achieve that

effect_eigen <- get_eigenvalue(effect)
effect_eigen ## PCA1 eigenvalue = 3.8, PCA2 eigenvalue = 0.74
## Technically only need PCA1 to explain maximum amount of variation

x <- predict(effect) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2)]) ## not correlated (r = -0.0000000000000002)
## use these as two independent axes of effect traits

effect_traits <- cbind(effect_traits, x) ## use PCA1 and PCA2 as two traits to calculate effect trait diversity
effect_traits <- effect_traits[-c(9:11)]

## need to merge in site data, so each site has a separate community of species 
BBS_data <- read.csv("../Data/BBS_2004_2018.csv", header=TRUE)

## merge effect traits with BBS species data to match number of species
spp_list2 <- data.frame(ENGLISH_NAME = unique(effect_traits$ENGLISH_NAME)) ## 124 species
BBS_data <- merge(spp_list2, BBS_data, by="ENGLISH_NAME")
BBS_data <- droplevels(BBS_data)
length(unique(BBS_data$ENGLISH_NAME)) ## 83 species

## make this into a dataframe with species as columns and sites as rows filled with abundance at each site (including zeros if species isn't present)
BBS_data2 <- BBS_data %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(n=n()) %>% 
  spread(ENGLISH_NAME, n, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_data2[,(1)]
## add in numbers 1:200
site_match$site <- 1:200
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data2 <- BBS_data2[,-(1)]

## merge BBS_data with effect traits to make same number of species
spp_list <- data.frame(ENGLISH_NAME = unique(BBS_data$ENGLISH_NAME)) ## 83 species
effect_traits <- merge(spp_list, effect_traits, by="ENGLISH_NAME") ## 83 species
effect_traits <- effect_traits[c(1,7,8)] ## english name, PCA1 and PCA2

## code to confirm that there are 83 species with effect trait and abundance data  
# matchingNames <- spp_list$ENGLISH_NAME %in% spp_list2$ENGLISH_NAME
# spp_list$result <- matchingNames

## change row numbers to species names
rownames(effect_traits) <- effect_traits[, 1]
effect_traits <- effect_traits[, -1]

# find distance matrix 
d <- dist(as.matrix(effect_traits))
# apply hirarchical clustering
cluster <- hclust(d, method="average")
# plot the dendrogram
plot(cluster)
# calculate total branch length for each site
FD <- treedive(BBS_data2, cluster, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 83 species
functional_div <- data.frame(FD_effect=FD, gridref=site_match$GRIDREF)


######################################################################
####################### RESPONSE TRAITS ##############################
######################################################################

## select response only traits 
response_traits <- BBS_traits[, c(1,22,26,31:33)]
response_traits[c(2:6)] <- lapply(response_traits[c(2:6)], function(response_traits) c(scale(response_traits, center = TRUE, scale = TRUE))) ## centre and scale trait data
head(response_traits)
str(response_traits)
response_traits <- na.omit(response_traits) ## 122 species without thermal maximum, clutch size and HSI
## correlation matrix 
response_corr <- cor(response_traits[c(2:6)])
response_corr

## PCA of all response traits
response <- prcomp(response_traits[ ,2:6]) ## all response traits
summary(response) ## PC1 accounts for 65% of the total variation and PC2 33% of the total variation

response_eigen <- get_eigenvalue(response)
response_eigen ## PCA1 eigenvalue = 1.39, PCA2 eigenvalue = 0.84

x <- predict(response) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2)]) ## not correlated (r = -0.0000000000000002)
## use these as two independent axes of response traits

response_traits <- cbind(response_traits, x) ## use PCA1 and PCA2 as two traits to calculate effect trait diversity
response_traits <- response_traits[-c(9:11)]

## need to merge in site data, so each site has a separate community of species 
BBS_data <- read.csv("../Data/BBS_2004_2018.csv", header=TRUE)

## merge effect traits with BBS species data to match number of species
spp_list2 <- data.frame(ENGLISH_NAME = unique(response_traits$ENGLISH_NAME)) ## 122 species
BBS_data <- merge(spp_list2, BBS_data, by="ENGLISH_NAME")
BBS_data <- droplevels(BBS_data)
length(unique(BBS_data$ENGLISH_NAME)) ## 93 species with trait and abundance data

## make this into a dataframe with species as columns and sites as rows filled with abundance at each site (including zeros if species isn't present)
BBS_data2 <- BBS_data %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(n=n()) %>% 
  spread(ENGLISH_NAME, n, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_data2[,(1)]
## add in numbers 1:200
site_match$site <- 1:200
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data2 <- BBS_data2[,-(1)]

## merge BBS_data with effect traits to make same number of species
spp_list <- data.frame(ENGLISH_NAME = unique(BBS_data$ENGLISH_NAME)) ## 93 species
response_traits <- merge(spp_list, response_traits, by="ENGLISH_NAME") ## 93 species
response_traits <- response_traits[c(1,7,8)] ## english name, PCA1 and PCA2

## code to confirm that there are 93 species with response trait and abundance data  
# matchingNames <- spp_list$ENGLISH_NAME %in% spp_list2$ENGLISH_NAME
# spp_list$result <- matchingNames

## change row numbers to species names
rownames(response_traits) <- response_traits[, 1]
response_traits <- response_traits[, -1]

# find distance matrix 
d <- dist(as.matrix(response_traits))
# apply hirarchical clustering
cluster <- hclust(d, method="average")
# plot the dendrogram
plot(cluster)
# calculate total branch length for each site
FD <- treedive(BBS_data2, cluster, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 93 species
functional_div2 <- data.frame(FD_response=FD, gridref=site_match$GRIDREF)

##################################################################
####################### BOTH TRAITS ##############################
##################################################################

## select both traits only
both_traits <- BBS_traits[, c(1,27:29)]
both_traits[c(2:4)] <- lapply(both_traits[c(2:4)], function(both_traits) c(scale(both_traits, center = TRUE, scale = TRUE))) ## centre and scale trait data
head(both_traits)
str(both_traits)
both_traits <- na.omit(both_traits) ## 194 species
## correlation matrix 
both_corr <- cor(both_traits[c(2:4)])
both_corr

## PCA of all both traits
both <- prcomp(both_traits[ ,2:4]) ## all response traits
summary(both) ## PC1 accounts for 94% of the total variation 
## just use this one axis to represent both traits

both_eigen <- get_eigenvalue(both)
both_eigen ## PCA1 eigenvalue = 2.6, PCA2 eigenvalue = 0.364





############## find species which have all effect and response traits plus abundance data
spp_effect <- data.frame(ENGLISH_NAME = rownames(effect_traits)) ## 83
spp_response <- data.frame(ENGLISH_NAME=rownames(response_traits)) ## 93
BBS_data <- read.csv("../Data/BBS_2004_2018.csv", header=TRUE)
spp_abund <- data.frame(ENGLISH_NAME=unique(BBS_data$ENGLISH_NAME)) ## 124

total_spp <- merge(spp_effect, spp_response, by="ENGLISH_NAME")
total_spp <- merge(total_spp, spp_abund, by="ENGLISH_NAME") ## 67 SPECIES 
