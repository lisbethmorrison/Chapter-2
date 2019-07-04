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
## remove year and no_years
## take the species richness across all years for each site
## also remove TOT (don't need abundance to calculate FD dendrogram)
BBS_data <- BBS_data[-c(4:6)]
BBS_data <- unique(BBS_data)

BBS_final <- merge(BBS_data, effect_traits, by="ENGLISH_NAME", all=FALSE)
length(unique(BBS_final$GRIDREF)) # 200 sites
length(unique(BBS_final$ENGLISH_NAME)) # 83 species which have all trait and abundance data

## now calculate FD dendrogram for each site using PCA1 and PCA2
traits <- BBS_final[BBS_final$GRIDREF=="TQ7877",]
traits <- traits[c(1,9,10)] ## english name, PCA1 and PCA2

my.dist = as.matrix(dist(traits, method = "euclidean"))  ## distance matrix
rownames(my.dist) <- traits$ENGLISH_NAME
colnames(my.dist) <- traits$ENGLISH_NAME
my.dendro = hclust(as.dist(my.dist), method = "average") ## heirarchial clustering to divide clusters among species
plot(my.dendro)
FD <- treedive(traits, my.dendro, match.force=TRUE)



######################################################################
####################### RESPONSE TRAITS ##############################
######################################################################

## select response only traits 
response_traits <- BBS_traits[, c(1,22,26,31:33)]
response_traits[c(2:6)] <- lapply(response_traits[c(2:6)], function(response_traits) c(scale(response_traits, center = TRUE, scale = TRUE))) ## centre and scale trait data
head(response_traits)
str(response_traits)
response_traits <- na.omit(response_traits) ## 122 species without thermal maximum, clutch size and HHI
## correlation matrix 
response_corr <- cor(response_traits[c(2:6)])
response_corr

## PCA of all response traits
response <- prcomp(response_traits[ ,2:6]) ## all response traits
summary(response) ## PC1 accounts for 65% of the total variation and PC2 33% of the total variation

response_eigen <- get_eigenvalue(response)
response_eigen ## PCA1 eigenvalue = 1.39, PCA2 eigenvalue = 0.84

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
