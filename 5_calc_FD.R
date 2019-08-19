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
library(FD)

## read in trait data 
effect_traits <- read.csv("../Data/Analysis_data/Seed dispersal/seed_effect_traits.csv", header=TRUE)
response_traits <- read.csv("../Data/Analysis_data/Seed dispersal/seed_response_traits.csv", header=TRUE)
effect_both_traits <- read.csv("../Data/Analysis_data/Seed dispersal/seed_effect_both_traits.csv", header=TRUE)
response_both_traits <- read.csv("../Data/Analysis_data/Seed dispersal/seed_response_both_traits.csv", header=TRUE)
all_traits <- read.csv("../Data/Analysis_data/Seed dispersal/seed_all_traits.csv", header=TRUE)
BBS_data <- read.csv("../Data/BBS_abund_data/BBS_2004_2018.csv", header=TRUE)

## Effect traits:
  # Bill length/width/depth
  # Gape width

## Response traits
  # Species specialisation index
  # Species temperature index
  # Mean latitude  
  # Lifespan
  # Clutch size
  # Number of broods

## Both traits:
  # Kipp's distance
  # Wing/tail length
  # Tarsus length

####################################################################
####################### EFFECT TRAITS ##############################
####################################################################

## select effect only traits from effect and response file
all_effect_traits <- all_traits[,c(1,4:7)] ## keep name and effect traits
effect_traits <- effect_traits[,-c(2:3)] ## remove CBC code and body size

## centre and scale trait data
effect_traits[c(2:5)] <- lapply(effect_traits[c(2:5)], function(effect_traits) c(scale(effect_traits, center = TRUE, scale = TRUE))) 
all_effect_traits[c(2:5)] <- lapply(all_effect_traits[c(2:5)], function(all_effect_traits) c(scale(all_effect_traits, center = TRUE, scale = TRUE))) 

## correlation matrix 
effect_corr1 <- cor(effect_traits[c(2:5)])
effect_corr1
effect_corr2 <- cor(all_effect_traits[c(2:5)])
effect_corr2

#############################################################################################
## First run this on species with effect traits (n=32)
## PCA of all effect traits
effect <- prcomp(effect_traits[ ,2:5]) ## bill length, width and depth, gape width and tarsus length
summary(effect) # PC1 accounts for 76.7% of the total variation and PC2 15.5% of the total variation
## use these two axes to calculatse FD of effect traits

## check eigenvalues to determine the number of principal components to be considered
## An eigenvalue > 1 indicates that PCs account for more variance than accounted by one of the original variables in standardized data
## This is commonly used as a cutoff point for which PCs are retained
## You can also limit the number of component to that number that accounts for a certain fraction of the total variance
## For example, if you are satisfied with 70% of the total variance explained then use the number of components to achieve that

effect_eigen <- get_eigenvalue(effect)
effect_eigen ## PCA1 eigenvalue = 3.06, PCA2 eigenvalue = 0.62
## Technically only need PCA1 to explain maximum amount of variation

x <- predict(effect) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2)]) ## not correlated (r = 0.0000000000000004832362)
## use these as two independent axes of effect traits

effect_traits <- cbind(effect_traits, x) ## use PCA1 and PCA2 as two traits to calculate effect trait diversity
## correlation between PCA1 and effect traits
x_corr2 <- cor(effect_traits[c(2:6)]) 
x_corr2 ## strong positive correlation with all effect traits
## correlation between PCA2 and effect traits
x_corr3 <- cor(effect_traits[c(2:5,7)])  ## low correlations
x_corr3 ## strong negative correlation with all effect traits
effect_traits <- effect_traits[,c(1,6:7)] ## keep only name, PCA1 and PCA2 columns

## need to merge in site data and seed dispersing species, so each site has a separate community of species 

## merge effect traits with BBS species data to match number of species
spp_list <- data.frame(ENGLISH_NAME = unique(effect_traits$ENGLISH_NAME)) ## 32 species
BBS_data_effect <- merge(spp_list, BBS_data, by="ENGLISH_NAME")
BBS_data_effect <- droplevels(BBS_data_effect)
length(unique(BBS_data_effect$ENGLISH_NAME)) ## 32 species

## make this into a dataframe with species as columns and sites as rows filled no years of data (i.e. presence/absence at each site)
BBS_data2 <- BBS_data_effect %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(n=n()) %>% 
  spread(ENGLISH_NAME, n, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_data2[,(1)]
## add in numbers 1:199
site_match$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data2 <- BBS_data2[,-(1)]

## change row numbers to species names
rownames(effect_traits) <- effect_traits[, 1]
effect_traits <- effect_traits[, -1]

#############################################################################################
## Now run this on species with all traits (n=28)
## PCA of all effect traits
effect2 <- prcomp(all_effect_traits[ ,2:5]) ## bill length, width and depth, gape width and tarsus length
summary(effect2) # PC1 accounts for 72% of the total variation and PC2 18% of the total variation
## use these two axes to calculatse FD of effect traits

## check eigenvalues to determine the number of principal components to be considered
## An eigenvalue > 1 indicates that PCs account for more variance than accounted by one of the original variables in standardized data
## This is commonly used as a cutoff point for which PCs are retained
## You can also limit the number of component to that number that accounts for a certain fraction of the total variance
## For example, if you are satisfied with 70% of the total variance explained then use the number of components to achieve that

effect_eigen2 <- get_eigenvalue(effect2)
effect_eigen2 ## PCA1 eigenvalue = 2.89, PCA2 eigenvalue = 0.72
## Technically only need PCA1 to explain maximum amount of variation

x2 <- predict(effect2) ## save scores from PCA
x2 <- as.data.frame(x2, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr2 <- cor(x2[c(1,2)]) ## not correlated (r = -0.0000000000000006828997)
## use these as two independent axes of effect traits

all_effect_traits1 <- cbind(all_effect_traits, x2) ## use PCA1 and PCA2 as two traits to calculate effect trait diversity
all_effect_traits1 <- all_effect_traits1[,c(1,6:7)] ## keep only name, PCA1 and PCA2 columns

## need to merge in site data and seed dispersing species, so each site has a separate community of species 

## merge effect traits with BBS species data to match number of species
spp_list2 <- data.frame(ENGLISH_NAME = unique(all_effect_traits1$ENGLISH_NAME)) ## 28 species
BBS_data_effect_response <- merge(spp_list2, BBS_data, by="ENGLISH_NAME")
BBS_data_effect_response <- droplevels(BBS_data_effect_response)
length(unique(BBS_data_effect_response$ENGLISH_NAME)) ## 28 species
length(unique(BBS_data_effect_response$GRIDREF)) ## this removes one site from the analysis (199 sites)

## make this into a dataframe with species as columns and sites as rows filled no years of data (i.e. presence/absence at each site)
BBS_data3 <- BBS_data_effect_response %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(n=n()) %>% 
  spread(ENGLISH_NAME, n, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match2 <- BBS_data3[,(1)]
## add in numbers 1:200
site_match2$site <- 1:199 ## less species == one site dropped
## remove first column (which is actual gridrefs) and leave sites as numbers 1-199
BBS_data3 <- BBS_data3[,-(1)]

## change row numbers to species names
rownames(all_effect_traits1) <- all_effect_traits1[, 1]
all_effect_traits1 <- all_effect_traits1[, -1]

################################################
########## Total Branch Length FD ##############
################################################

######### 32 species which ONLY have effect traits

# find distance matrix 
d <- dist(as.matrix(effect_traits))
# apply hirarchical clustering
cluster <- hclust(d, method="average")
# plot the dendrogram
plot(cluster)
# calculate total branch length for each site
FD <- treedive(BBS_data2, cluster, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 53 species
functional_div_effect1 <- data.frame(gridref=site_match$GRIDREF, FD_effect=FD)
## functional diversity total branch length 
## 53 species at 200 sites
## these are species which have ONLY effect traits 

## save file
write.csv(functional_div_effect1, file="../Data/Analysis_data/Seed dispersal/FD_effect_32spp.csv", row.names=FALSE)

######### 28 species which have BOTH effect and response traits
# find distance matrix 
d2 <- dist(as.matrix(all_effect_traits1))
# apply hirarchical clustering
cluster2 <- hclust(d2, method="average")
# plot the dendrogram
plot(cluster2)
# calculate total branch length for each site
FD2 <- treedive(BBS_data3, cluster2, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (199 in total) for 47 species
functional_div_effect2 <- data.frame(gridref=site_match2$GRIDREF, FD_effect2=FD2)
## functional diversity total branch length 
## 47 species at 199 sites
## these are species which have BOTH effect and response traits 

## save file
write.csv(functional_div_effect2, file="../Data/Analysis_data/Seed dispersal/FD_effect_28spp.csv", row.names=FALSE)

################################################
########## Multivariate FD metrics #############
################################################

## First run this on species with effect traits (n=32)

## create abundance matrix
## make this into a dataframe with species as columns and sites as rows filled with average abundance over time at each site
BBS_data4 <- BBS_data_effect %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(mean_abund=mean(TOT)) %>% 
  spread(ENGLISH_NAME, mean_abund, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_data4[,(1)]
## add in numbers 1:200
site_match$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data4 <- BBS_data4[,-(1)]

FD_multi_effect2 <- dbFD(effect_traits, BBS_data4)
## one site (85 = SU0451) doesn't work as there is only one species in the site (carrion crow)
## other species present at the site had missing gape width data so removed
FD_multi_effect <- data.frame(gridref=site_match$GRIDREF, FDVe=FD_multi_effect2$FEve, FDis=FD_multi_effect2$FDis, 
                              FRic=FD_multi_effect2$FRic, FDiv=FD_multi_effect2$FDiv)
## this is calculated on the two PCA axes derived from 5 effect traits
write.csv(FD_multi_effect, file="../Data/Analysis_data/Seed dispersal/FD_multi_effect_32spp.csv", row.names=FALSE)

#################################################################################################################
## Now run this on species with effect and response traits (n=28)

## create abundance matrix
# make this into a dataframe with species as columns and sites as rows filled with average abundance over time at each site
BBS_data5 <- BBS_data_effect_response %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(mean_abund=mean(TOT)) %>% 
  spread(ENGLISH_NAME, mean_abund, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match2 <- BBS_data5[,(1)]
## add in numbers 1:200
site_match2$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data5 <- BBS_data5[,-(1)]

FD_multi_effect_response2 <- dbFD(all_effect_traits1, BBS_data5)
## one site (85 = SU0451) doesn't work as there is only one species in the site (carrion crow)
## other species present at the site had missing gape width data so removed
FD_multi_effect_response <- data.frame(gridref=site_match2$GRIDREF, FDVe=FD_multi_effect_response2$FEve, FDis=FD_multi_effect_response2$FDis, 
                              FRic=FD_multi_effect_response2$FRic, FDiv=FD_multi_effect_response2$FDiv)
## this is calculated on the two PCA axes derived from 5 effect traits
write.csv(FD_multi_effect_response, file="../Data/Analysis_data/Seed dispersal/FD_multi_effect_28spp.csv", row.names=FALSE)

######################################################################
####################### RESPONSE TRAITS ##############################
######################################################################

## select response only traits from effect and response file
all_response_traits <- all_traits[,c(1,8:13)] ## keep name and response traits
response_traits <- response_traits[,-c(2)] ## remove CBC code 

## centre and scale trait data
response_traits[c(2:7)] <- lapply(response_traits[c(2:7)], function(response_traits) c(scale(response_traits, center = TRUE, scale = TRUE))) 
all_response_traits[c(2:7)] <- lapply(all_response_traits[c(2:7)], function(all_response_traits) c(scale(all_response_traits, center = TRUE, scale = TRUE))) 

## correlation matrix 
response_corr1 <- cor(response_traits[c(2:7)])
response_corr1
response_corr2 <- cor(all_response_traits[c(2:7)])
response_corr2
## not very correlated, highest r = -0.47

################################################
## Sort out data on species with only response traits (n=38)
## PCA of all response traits
response <- prcomp(response_traits[ ,2:7]) ## all response traits
summary(response) ## PC1 accounts for 29% of the total variation and PC2 24% of the total variation, PCA3 15%

response_eigen <- get_eigenvalue(response)
response_eigen ## PCA1 eigenvalue = 1.76, PCA2 eigenvalue = 1.45, PCA3 eigenvalue = 1.12 (keep all 3 as they are above 1)

x <- predict(response) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe

## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2,3)]) ## not correlated
## use these as two independent axes of response traits

response_traits <- cbind(response_traits, x) ## use PCA1, PCA2 and PCA3 as two 3 to calculate response trait diversity
## correlation between PCA1, PCA2 and PCA3 with response traits
x_corr2 <- cor(response_traits[c(2:10)]) 
x_corr2 
## PCA1 correlates with clutch size, brood size and STI
## PCA2 correlates with mean latitude and SSI
## PCA3 correlates with maximum longevity
response_traits <- response_traits[,c(1,8:10)] ## keep only name, PCA1, PCA2 amd PCA3 columns

## merge response traits with BBS species data to match number of species
spp_list <- data.frame(ENGLISH_NAME = unique(response_traits$ENGLISH_NAME)) ## 38 species
BBS_data_response <- merge(spp_list, BBS_data, by="ENGLISH_NAME")
BBS_data_response <- droplevels(BBS_data_response)
length(unique(BBS_data_response$ENGLISH_NAME)) ## 38 species with trait and abundance data

## make this into a dataframe with species as columns and sites as rows filled no years of data (i.e. presence/absence at each site)
BBS_data6<- BBS_data_response %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(n=n()) %>% 
  spread(ENGLISH_NAME, n, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_data6[,(1)]
## add in numbers 1:200
site_match$site <- 1:200
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data6 <- BBS_data6[,-(1)]

## change row numbers to species names
rownames(response_traits) <- response_traits[, 1]
response_traits <- response_traits[, -1]

################################################
## Sort out data for species with effect and response traits (n=28)

## PCA of all response traits
response2 <- prcomp(all_response_traits[ ,2:7]) ## all response traits
summary(response2) ## PC1 accounts for 31% of the total variation and PC2 25% of the total variation, PCA3 18%

response_eigen2 <- get_eigenvalue(response2)
response_eigen2 ## PCA1 eigenvalue = 1.88, PCA2 eigenvalue = 1.51, PCA3 eigenvalue = 1.12

x2 <- predict(response2) ## save scores from PCA
x2 <- as.data.frame(x2, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr2 <- cor(x2[c(1,2,3)]) ## not correlated
## use these as two independent axes of response traits

all_response_traits2 <- cbind(all_response_traits, x2) ## use PCA1, PCA2 and PCA3 as 3 traits to calculate response trait diversity
all_response_traits2 <- all_response_traits2[,c(1,8:10)]

## merge effect traits with BBS species data to match number of species
spp_list2 <- data.frame(ENGLISH_NAME = unique(all_response_traits2$ENGLISH_NAME)) ## 28 species
BBS_data_response_effect <- merge(spp_list2, BBS_data, by="ENGLISH_NAME")
BBS_data_response_effect <- droplevels(BBS_data_response_effect)
length(unique(BBS_data_response_effect$ENGLISH_NAME)) ## 28 species with trait and abundance data

## make this into a dataframe with species as columns and sites as rows filled no years of data (i.e. presence/absence at each site)
BBS_data7 <- BBS_data_response_effect %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(n=n()) %>% 
  spread(ENGLISH_NAME, n, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match2 <- BBS_data7[,(1)]
## add in numbers 1:199 (47 species reduces number of sites to 199)
site_match2$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data7 <- BBS_data7[,-(1)]

## change row numbers to species names
rownames(all_response_traits2) <- all_response_traits2[, 1]
all_response_traits2 <- all_response_traits2[, -1]

################################################
########## Total Branch Length FD ##############
################################################

################################################
## First run this on species with only response traits (n=66)

# find distance matrix 
d <- dist(as.matrix(response_traits))
# apply hirarchical clustering
cluster <- hclust(d, method="average")
# plot the dendrogram
plot(cluster)
# calculate total branch length for each site
FD <- treedive(BBS_data6, cluster, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 93 species
functional_div_response <- data.frame(gridref=site_match$GRIDREF, FD_response=FD)
## functional diversity total branch length 
## 66 species at 200 sites
## these are species which have ONLY response traits 

## save file
write.csv(functional_div_response, file="../Data/Analysis_data/Seed dispersal/FD_response_38spp.csv", row.names=FALSE)

################################################
## First run this on species with both effect and response traits (n=28)

# find distance matrix 
d2 <- dist(as.matrix(all_response_traits2))
# apply hirarchical clustering
cluster2 <- hclust(d2, method="average")
# plot the dendrogram
plot(cluster2)
# calculate total branch length for each site
FD2 <- treedive(BBS_data7, cluster2, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 93 species
functional_div_response2 <- data.frame(gridref=site_match2$GRIDREF, FD_response=FD2)
## functional diversity total branch length 
## 28 species at 199 sites
## these are species which have BOTH response and effect traits 

## save file
write.csv(functional_div_response2, file="../Data/Analysis_data/Seed dispersal/FD_response_28spp.csv", row.names=FALSE)


################################################
########## Multivariate FD metrics #############
################################################

## First run this on species with response traits (n=38)

## create abundance matrix
## make this into a dataframe with species as columns and sites as rows filled with average abundance over time at each site
BBS_data8 <- BBS_data_response %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(mean_abund=mean(TOT)) %>% 
  spread(ENGLISH_NAME, mean_abund, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_data8[,(1)]
## add in numbers 1:200
site_match$site <- 1:200
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data8 <- BBS_data8[,-(1)]

FD_multi_response2 <- dbFD(response_traits, BBS_data8)
## one site (85 = SU0451) doesn't work as there is only one species in the site (carrion crow)
## other species present at the site had missing gape width data so removed
FD_multi_response <- data.frame(gridref=site_match$GRIDREF, FDVe=FD_multi_response2$FEve, FDis=FD_multi_response2$FDis, 
                              FRic=FD_multi_response2$FRic, FDiv=FD_multi_response2$FDiv)
## this is calculated on the two PCA axes derived from 5 effect traits
write.csv(FD_multi_response, file="../Data/Analysis_data/Seed dispersal/FD_multi_response_38spp.csv", row.names=FALSE)
## all 200 sites and 38 species

#################################################################################################################
## Now run this on species with effect and response traits (n=28)

## create abundance matrix
# make this into a dataframe with species as columns and sites as rows filled with average abundance over time at each site
BBS_data9 <- BBS_data_response_effect %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(mean_abund=mean(TOT)) %>% 
  spread(ENGLISH_NAME, mean_abund, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match2 <- BBS_data9[,(1)]
## add in numbers 1:200
site_match2$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data9 <- BBS_data9[,-(1)]

FD_multi_response_effect2 <- dbFD(all_response_traits2, BBS_data9)
## one site (85 = SU0451) doesn't work as there is only one species in the site (carrion crow)
## other species present at the site had missing gape width data so removed
FD_multi_response_effect <- data.frame(gridref=site_match2$GRIDREF, FDVe=FD_multi_response_effect2$FEve, FDis=FD_multi_response_effect2$FDis, 
                                       FRic=FD_multi_response_effect2$FRic, FDiv=FD_multi_response_effect2$FDiv)
## this is calculated on the two PCA axes derived from 5 effect traits
write.csv(FD_multi_response_effect, file="../Data/Analysis_data/Seed dispersal/FD_multi_response_28spp.csv", row.names=FALSE)
## 199 sites and 28 species


#############################################################################
####################### BOTH AND EFFECT TRAITS ##############################
#############################################################################

## We want to combine both traits and effect traits => calculate FD metrics using these
## Then combine both and response traits => calculate FD metrics using these
## Repeat this with subset of species (n=28)

## select both only traits from all traits file
all_effect_both_traits <- all_traits[,c(1,4:7,14:17)] ## keep name and both traits (28spp)
effect_both_traits <- effect_both_traits[,-c(2)] ## remove CBC code (32spp)

## centre and scale trait data
all_effect_both_traits[c(2:9)] <- lapply(all_effect_both_traits[c(2:9)], function(all_effect_both_traits) c(scale(all_effect_both_traits, center = TRUE, scale = TRUE))) 
effect_both_traits[c(3:10)] <- lapply(effect_both_traits[c(3:10)], function(effect_both_traits) c(scale(effect_both_traits, center = TRUE, scale = TRUE))) 

## correlation matrix 
both_corr1 <- cor(all_effect_both_traits[c(2:9)])
both_corr1
both_corr3 <- cor(effect_both_traits[c(3:10)])
both_corr3
## high correlations (range from 0.4 to 0.95)

################################################
## Sort out data on species with only effect AND both traits (n=32)
## PCA of all both traits
effect_both <- prcomp(effect_both_traits[ ,3:10]) ## all response traits
summary(effect_both) ## PC1 accounts for 73% of the total variation and PC2 11% of the total variation

effect_both_eigen <- get_eigenvalue(effect_both)
effect_both_eigen ## PCA1 eigenvalue = 5.87, PCA2 eigenvalue = 0.92 (take the first two axes)

x <- predict(effect_both) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe

## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2)]) ## not correlated
## use these as two independent axes of response traits

effect_both_traits <- cbind(effect_both_traits, x) ## use PCA1, PCA2 and PCA3 as two 3 to calculate response trait diversity
## correlation between PCA1 and PCA2 with effect traits
x_corr2 <- cor(effect_both_traits[c(2:12)]) 
x_corr2 
## PCA1 high correlations with all traits
## PCA2 low correlations with all traits
effect_both_traits <- effect_both_traits[,c(1,11:12)] ## keep only name, PCA1 and PCA2 columns

## merge response traits with BBS species data to match number of species
spp_list <- data.frame(ENGLISH_NAME = unique(effect_both_traits$ENGLISH_NAME)) ## 32 species
BBS_data_effect_both <- merge(spp_list, BBS_data, by="ENGLISH_NAME")
BBS_data_effect_both <- droplevels(BBS_data_effect_both)
length(unique(BBS_data_effect_both$ENGLISH_NAME)) ## 32 species with trait and abundance data

## make this into a dataframe with species as columns and sites as rows filled no years of data (i.e. presence/absence at each site)
BBS_data10<- BBS_data_effect_both %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(n=n()) %>% 
  spread(ENGLISH_NAME, n, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_data10[,(1)]
## add in numbers 1:200
site_match$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data10 <- BBS_data10[,-(1)]

## change row numbers to species names
rownames(effect_both_traits) <- effect_both_traits[, 1]
effect_both_traits <- effect_both_traits[, -1]

################################################
## Sort out data for species with all traits (n=28)

## PCA of all response traits
effect_both2 <- prcomp(all_effect_both_traits[ ,2:9]) ## all response traits
summary(effect_both2) ## PC1 accounts for 67% of the total variation and PC2 13% of the total variation

effect_both_eigen2 <- get_eigenvalue(effect_both2)
effect_both_eigen2 ## PCA1 eigenvalue = 5.4, PCA2 eigenvalue = 1.07, PCA3 eigenvalue = 0.81 (take first two)

x2 <- predict(effect_both2) ## save scores from PCA
x2 <- as.data.frame(x2, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr2 <- cor(x2[c(1,2)]) ## not correlated
## use these as two independent axes of response traits

all_effect_both_traits <- cbind(all_effect_both_traits, x2) ## use PCA1, PCA2 and PCA3 as 3 traits to calculate response trait diversity
all_effect_both_traits <- all_effect_both_traits[,c(1,10:11)]

## merge effect traits with BBS species data to match number of species
spp_list2 <- data.frame(ENGLISH_NAME = unique(all_effect_both_traits$ENGLISH_NAME)) ## 28 species
BBS_data_all_effect_both <- merge(spp_list2, BBS_data, by="ENGLISH_NAME")
BBS_data_all_effect_both <- droplevels(BBS_data_all_effect_both)
length(unique(BBS_data_all_effect_both$ENGLISH_NAME)) ## 28 species with trait and abundance data

## make this into a dataframe with species as columns and sites as rows filled no years of data (i.e. presence/absence at each site)
BBS_data11 <- BBS_data_all_effect_both %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(n=n()) %>% 
  spread(ENGLISH_NAME, n, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match2 <- BBS_data11[,(1)]
## add in numbers 1:199 
site_match2$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data11 <- BBS_data11[,-(1)]

## change row numbers to species names
rownames(all_effect_both_traits) <- all_effect_both_traits[, 1]
all_effect_both_traits <- all_effect_both_traits[, -1]

################################################
########## Total Branch Length FD ##############
################################################

################################################
## First run this on species with only response traits (n=66)

# find distance matrix 
d <- dist(as.matrix(effect_both_traits))
# apply hirarchical clustering
cluster <- hclust(d, method="average")
# plot the dendrogram
plot(cluster)
# calculate total branch length for each site
FD <- treedive(BBS_data10, cluster, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 93 species
functional_div_effect_both <- data.frame(gridref=site_match$GRIDREF, FD_effect_both=FD)
## functional diversity total branch length 
## 32 species at 199 sites
## these are species which have effect and both traits

## save file
write.csv(functional_div_effect_both, file="../Data/Analysis_data/Seed dispersal/FD_effect_both_32spp.csv", row.names=FALSE)

################################################
## Now run this on species with all traits (n=28)

# find distance matrix 
d2 <- dist(as.matrix(all_effect_both_traits))
# apply hirarchical clustering
cluster2 <- hclust(d2, method="average")
# plot the dendrogram
plot(cluster2)
# calculate total branch length for each site
FD2 <- treedive(BBS_data11, cluster2, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 93 species
functional_div_effect_both2 <- data.frame(gridref=site_match2$GRIDREF, FD_effect_both=FD2)
## functional diversity total branch length 
## 28 species at 199 sites
## these are species which have BOTH response and effect traits 

## save file
write.csv(functional_div_effect_both2, file="../Data/Analysis_data/Seed dispersal/FD_effect_both_28spp.csv", row.names=FALSE)


################################################
########## Multivariate FD metrics #############
################################################

## First run this on species with effect and both traits (n=32)

## create abundance matrix
## make this into a dataframe with species as columns and sites as rows filled with average abundance over time at each site
BBS_data12 <- BBS_data_effect_both %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(mean_abund=mean(TOT)) %>% 
  spread(ENGLISH_NAME, mean_abund, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_data12[,(1)]
## add in numbers 1:200
site_match$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data12 <- BBS_data12[,-(1)]

FD_multi_effect_both2 <- dbFD(effect_both_traits, BBS_data12)
## one site (85 = SU0451) doesn't work as there is only one species in the site (carrion crow)
## other species present at the site had missing gape width data so removed
FD_multi_effect_both <- data.frame(gridref=site_match$GRIDREF, FDVe=FD_multi_effect_both2$FEve, FDis=FD_multi_effect_both2$FDis, 
                                FRic=FD_multi_effect_both2$FRic, FDiv=FD_multi_effect_both2$FDiv)
## this is calculated on the two PCA axes derived from 5 effect traits
write.csv(FD_multi_effect_both, file="../Data/Analysis_data/Seed dispersal/FD_multi_effect_both_32spp.csv", row.names=FALSE)
## all 200 sites and 38 species

#################################################################################################################
## Now run this on species with all traits (n=28)

## create abundance matrix
# make this into a dataframe with species as columns and sites as rows filled with average abundance over time at each site
BBS_data13 <- BBS_data_all_effect_both %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(mean_abund=mean(TOT)) %>% 
  spread(ENGLISH_NAME, mean_abund, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match2 <- BBS_data13[,(1)]
## add in numbers 1:200
site_match2$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data13 <- BBS_data13[,-(1)]

FD_multi_all_effect_both2 <- dbFD(all_effect_both_traits, BBS_data13)
## one site (85 = SU0451) doesn't work as there is only one species in the site (carrion crow)
## other species present at the site had missing gape width data so removed
FD_multi_all_effect_both <- data.frame(gridref=site_match2$GRIDREF, FDVe=FD_multi_all_effect_both2$FEve, FDis=FD_multi_all_effect_both2$FDis, 
                                       FRic=FD_multi_all_effect_both2$FRic, FDiv=FD_multi_all_effect_both2$FDiv)
## this is calculated on the two PCA axes derived from 5 effect traits
write.csv(FD_multi_all_effect_both, file="../Data/Analysis_data/Seed dispersal/FD_multi_effect_both_28spp.csv", row.names=FALSE)
## 199 sites and 28 species

###############################################################################
####################### BOTH AND RESPONSE TRAITS ##############################
###############################################################################

## We want to combine both traits and effect traits => calculate FD metrics using these
## Then combine both and response traits => calculate FD metrics using these
## Repeat this with subset of species (n=28)

## select both only traits from all traits file
all_response_both_traits <- all_traits[,c(1,8:17)] ## keep name and both traits (28spp)
response_both_traits <- response_both_traits[,-c(2)] ## remove CBC code (38spp)

## centre and scale trait data
all_response_both_traits[c(2:11)] <- lapply(all_response_both_traits[c(2:11)], function(all_response_both_traits) c(scale(all_response_both_traits, center = TRUE, scale = TRUE))) 
response_both_traits[c(2:11)] <- lapply(response_both_traits[c(2:11)], function(response_both_traits) c(scale(response_both_traits, center = TRUE, scale = TRUE))) 

## correlation matrix 
both_corr <- cor(all_response_both_traits[c(2:11)])
both_corr
both_corr2 <- cor(response_both_traits[c(2:11)])
both_corr2
## high correlations (range from 0.4 to 0.95)

################################################
## Sort out data on species with only response AND both traits (n=38)
## PCA of all both traits
response_both <- prcomp(response_both_traits[ ,2:11]) ## all response traits
summary(response_both) ## PC1 accounts for 37% of the total variation and PC2 18.5% of the total variation, PCA3 14.6%

response_both_eigen <- get_eigenvalue(response_both)
response_both_eigen ## PCA1 eigenvalue = 3.69, PCA2 eigenvalue = 1.85, PCA3 = 1.46, PCA4 = 1.05 (take the first 4 axes)

x <- predict(response_both) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe

## correlation between PCA1, PCA2, PCA3 and PCA4
x_corr <- cor(x[c(1,2,3,4)]) ## not correlated
## use these as 4 independent axes of response traits

response_both_traits <- cbind(response_both_traits, x) ## use PCA1, PCA2 and PCA3 as two 3 to calculate response trait diversity
## correlation between PCA1 and PCA2 with effect traits
x_corr2 <- cor(response_both_traits[c(2:15)]) 
x_corr2 
response_both_traits <- response_both_traits[,c(1,12:15)] ## keep only name, PCA1 and PCA2 columns

## merge response traits with BBS species data to match number of species
spp_list <- data.frame(ENGLISH_NAME = unique(response_both_traits$ENGLISH_NAME)) ## 38 species
BBS_data_response_both <- merge(spp_list, BBS_data, by="ENGLISH_NAME")
BBS_data_response_both <- droplevels(BBS_data_response_both)
length(unique(BBS_data_response_both$ENGLISH_NAME)) ## 38 species with trait and abundance data

## make this into a dataframe with species as columns and sites as rows filled no years of data (i.e. presence/absence at each site)
BBS_data14<- BBS_data_response_both %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(n=n()) %>% 
  spread(ENGLISH_NAME, n, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_data14[,(1)]
## add in numbers 1:200
site_match$site <- 1:200
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data14 <- BBS_data14[,-(1)]

## change row numbers to species names
rownames(response_both_traits) <- response_both_traits[, 1]
response_both_traits <- response_both_traits[, -1]

################################################
## Sort out data for species with all traits (n=28)

## PCA of all response traits
response_both2 <- prcomp(all_response_both_traits[ ,2:11]) ## all response traits
summary(response_both2) ## PC1 accounts for 35% of the total variation and PC2 20% of the total variation, PCA3 15%

response_both_eigen2 <- get_eigenvalue(response_both2)
response_both_eigen2 ## PCA1 eigenvalue = 3.5, PCA2 eigenvalue = 2.06, PCA3 eigenvalue = 1.52, PCA4 = 1.01 (take first 4)

x2 <- predict(response_both2) ## save scores from PCA
x2 <- as.data.frame(x2, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr2 <- cor(x2[c(1,2,3,4)]) ## not correlated
## use these as two independent axes of response traits

all_response_both_traits <- cbind(all_response_both_traits, x2) ## use PCA1, PCA2 and PCA3 as 3 traits to calculate response trait diversity
all_response_both_traits <- all_response_both_traits[,c(1,12:15)]

## merge effect traits with BBS species data to match number of species
spp_list2 <- data.frame(ENGLISH_NAME = unique(all_response_both_traits$ENGLISH_NAME)) ## 28 species
BBS_data_all_response_both <- merge(spp_list2, BBS_data, by="ENGLISH_NAME")
BBS_data_all_response_both <- droplevels(BBS_data_all_response_both)
length(unique(BBS_data_all_response_both$ENGLISH_NAME)) ## 28 species with trait and abundance data

## make this into a dataframe with species as columns and sites as rows filled no years of data (i.e. presence/absence at each site)
BBS_data15 <- BBS_data_all_response_both %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(n=n()) %>% 
  spread(ENGLISH_NAME, n, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match2 <- BBS_data15[,(1)]
## add in numbers 1:199 
site_match2$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data15 <- BBS_data15[,-(1)]

## change row numbers to species names
rownames(all_response_both_traits) <- all_response_both_traits[, 1]
all_response_both_traits <- all_response_both_traits[, -1]

################################################
########## Total Branch Length FD ##############
################################################

################################################
## First run this on species with both and response traits (n=38)

# find distance matrix 
d <- dist(as.matrix(response_both_traits))
# apply hirarchical clustering
cluster <- hclust(d, method="average")
# plot the dendrogram
plot(cluster)
# calculate total branch length for each site
FD <- treedive(BBS_data14, cluster, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 93 species
functional_div_response_both <- data.frame(gridref=site_match$GRIDREF, FD_response_both=FD)
## functional diversity total branch length 
## 38 species at 200 sites
## these are species which have effect and both traits

## save file
write.csv(functional_div_response_both, file="../Data/Analysis_data/Seed dispersal/FD_response_both_38spp.csv", row.names=FALSE)

################################################
## Now run this on species with all traits (n=28)

# find distance matrix 
d2 <- dist(as.matrix(all_response_both_traits))
# apply hirarchical clustering
cluster2 <- hclust(d2, method="average")
# plot the dendrogram
plot(cluster2)
# calculate total branch length for each site
FD2 <- treedive(BBS_data15, cluster2, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 93 species
functional_div_response_both2 <- data.frame(gridref=site_match2$GRIDREF, FD_response_both=FD2)
## functional diversity total branch length 
## 28 species at 199 sites
## these are species which have BOTH response and effect traits 

## save file
write.csv(functional_div_response_both2, file="../Data/Analysis_data/Seed dispersal/FD_response_both_28spp.csv", row.names=FALSE)


################################################
########## Multivariate FD metrics #############
################################################

## First run this on species with response and both traits (n=38)

## create abundance matrix
## make this into a dataframe with species as columns and sites as rows filled with average abundance over time at each site
BBS_data16 <- BBS_data_response_both %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(mean_abund=mean(TOT)) %>% 
  spread(ENGLISH_NAME, mean_abund, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_data16[,(1)]
## add in numbers 1:200
site_match$site <- 1:200
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data16 <- BBS_data16[,-(1)]

FD_multi_response_both2 <- dbFD(response_both_traits, BBS_data16)
## one site (85 = SU0451) doesn't work as there is only one species in the site (carrion crow)
## other species present at the site had missing gape width data so removed
FD_multi_response_both <- data.frame(gridref=site_match$GRIDREF, FDVe=FD_multi_response_both2$FEve, FDis=FD_multi_response_both2$FDis, 
                                   FRic=FD_multi_response_both2$FRic, FDiv=FD_multi_response_both2$FDiv)
## this is calculated on the two PCA axes derived from 5 effect traits
write.csv(FD_multi_response_both, file="../Data/Analysis_data/Seed dispersal/FD_multi_response_both_38spp.csv", row.names=FALSE)
## all 200 sites and 38 species

#################################################################################################################
## Now run this on species with all traits (n=28)

## create abundance matrix
# make this into a dataframe with species as columns and sites as rows filled with average abundance over time at each site
BBS_data17 <- BBS_data_all_response_both %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(mean_abund=mean(TOT)) %>% 
  spread(ENGLISH_NAME, mean_abund, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match2 <- BBS_data17[,(1)]
## add in numbers 1:200
site_match2$site <- 1:199
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_data17 <- BBS_data17[,-(1)]

FD_multi_all_response_both2 <- dbFD(all_response_both_traits, BBS_data17)
## one site (85 = SU0451) doesn't work as there is only one species in the site (carrion crow)
## other species present at the site had missing gape width data so removed
FD_multi_all_response_both <- data.frame(gridref=site_match2$GRIDREF, FDVe=FD_multi_all_response_both2$FEve, FDis=FD_multi_all_response_both2$FDis, 
                                       FRic=FD_multi_all_response_both2$FRic, FDiv=FD_multi_all_response_both2$FDiv)
## this is calculated on the two PCA axes derived from 5 effect traits
write.csv(FD_multi_all_response_both, file="../Data/Analysis_data/Seed dispersal/FD_multi_response_both_28spp.csv", row.names=FALSE)
## 199 sites and 28 species
