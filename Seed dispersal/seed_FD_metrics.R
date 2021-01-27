###########################################################
## Title: Calculate functional diversity of traits for seed dispersers
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
library(zoo)

## read in trait data 
seed_traits <- read.csv("../Data/Trait data/Seed dispersal/BBS_seed_trait_data.csv", header=TRUE) ## 40 species
BBS_data <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_seed_interpol_repeat_det.csv", header=TRUE)

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
# Thermal maximum

## Both traits:
# Kipp's distance
# Wing/tail length
# Tarsus length
# Body size
# HWI

## remove Crane from analysis (no SSI data)
seed_traits <- seed_traits[!seed_traits$ENGLISH_NAME=="Crane",]
BBS_data <- BBS_data[!BBS_data$ENGLISH_NAME=="Crane",]
BBS_data <- droplevels(BBS_data)

####################################################################
####################### EFFECT TRAITS ##############################
####################################################################

## select effect only traits from effect and response file
seed_effect <- seed_traits[,c(1,7:9,15)] ## keep name and effect traits

## centre and scale trait data
seed_effect[c(2:5)] <- lapply(seed_effect[c(2:5)], function(seed_effect) c(scale(seed_effect, center = TRUE, scale = TRUE))) 

## correlation matrix 
effect_corr <- cor(seed_effect[c(2:5)])
effect_corr ## all r>0.47

#############################################################################################
## High correlation between traits so run PCA on all effect traits
pca_effect <- prcomp(seed_effect[ ,2:5]) ## bill length, width and depth and gape width
summary(pca_effect) # PC1 accounts for 77% of the total variation and PC2 15.9% of the total variation
## use these two axes to calculate FD of effect traits

## check eigenvalues to determine the number of principal components to be considered
## An eigenvalue > 1 indicates that PCs account for more variance than accounted by one of the original variables in standardized data
## This is commonly used as a cutoff point for which PCs are retained
## You can also limit the number of component to that number that accounts for a certain fraction of the total variance
## For example, if you are satisfied with 70% of the total variance explained then use the number of components to achieve that

effect_eigen <- get_eigenvalue(pca_effect)
effect_eigen ## PCA1 eigenvalue = 3.08, PCA2 eigenvalue = 0.62, PCA3 eigenvalue = 0.19
## Technically only need PCA1 to explain maximum amount of variation

pca_effect$rotation

x <- predict(pca_effect) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2,3)]) ## not correlated r = ~0
## use these as two independent axes of effect traits

seed_effect <- cbind(seed_effect, x) ## use PCA1 and PCA2 as two traits to calculate effect trait diversity
## correlation between PCA1 and effect traits
x_corr2 <- cor(seed_effect[c(2:6)]) 
x_corr2 ## strong positive correlation with all effect traits
## correlation between PCA2 and effect traits
x_corr3 <- cor(seed_effect[c(2:5,7)])  #
x_corr3 # r=0.6 with bill length, weak correlations with others
seed_effect <- seed_effect[,c(1,6:8)] ## keep only name and all 4 PCA axes
## change row numbers to species names
rownames(seed_effect) <- seed_effect[, 1]
seed_effect <- seed_effect[, -1]

## plot 3D plot 
library(plot3D)
x <- seed_effect$PC1
y <- seed_effect$PC2
z <- seed_effect$PC3
species <- row.names(seed_effect)
## save
png("../Graphs/seed_effect_traits.png", width=60, height=60, units="mm", res=600)
scatter3D(x,y,z,box=TRUE,type="h",pch=16,cex=0.5,d=2,col="black",bty="b2",axes=TRUE,label=species, nticks=10,theta=40, 
          phi=30, xlab="", ylab="", 
          zlab="",par(mar = c(0,0,0,0),
          par(oma = c(0,0,0,0))))
dev.off()
text3D(x,y,z, cex=0.6, add=TRUE)
## x = beak size 77%
## y = Bill length vs width vs depth 16%
## z = Bill width vs gape width 5%

## put BBS data into a dataframe with species as columns and sites as rows filled with average abundance across iterations and year
## this is used as presence/absence for FD
## and used to weight FDis 
BBS_av_abund <- BBS_data %>%
  group_by(GRIDREF, ENGLISH_NAME) %>% 
  summarise(mean_abund=mean(TOT_det)) %>% 
  spread(ENGLISH_NAME, mean_abund, drop = FALSE, fill = 0) 
## save first column as a site_match dataframe
site_match <- BBS_av_abund[,(1)]
## add in numbers 1:199
site_match$site <- 1:200
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_av_abund <- BBS_av_abund[,-(1)]
# 
# ################################################
# ########## Total Branch Length FD ##############
# ################################################
# 
# # find distance matrix 
# d <- dist(as.matrix(seed_effect))
# # apply hirarchical clustering
# cluster <- hclust(d, method="average")
# # plot the dendrogram
# plot(cluster)
# # calculate total branch length for each site
# FD <- treedive(BBS_av_abund, cluster, match.force=TRUE) ## works!!!
# # dataframe with FD values and gridrefs (200 in total) for 53 species
# functional_div_effect <- data.frame(GRIDREF=site_match$GRIDREF, FD_effect=FD)
# ## functional diversity total branch length 
# ## 40 species at 200 sites using effect traits
# 
# ## save file
# write.csv(functional_div_effect, file="../Data/Analysis_data/Seed dispersal/FD_seed_effect.csv", row.names=FALSE)

################################################
########## Multivariate FD metrics #############
################################################

## use function dbFD to calculate multivariate FD metrics
FD_multi_effect <- dbFD(seed_effect, BBS_av_abund)
FD_multi_effect <- data.frame(GRIDREF=site_match$GRIDREF, FDVe=FD_multi_effect$FEve, FDis=FD_multi_effect$FDis, 
                              FRic=FD_multi_effect$FRic, FDiv=FD_multi_effect$FDiv)
FD_multi_effect <- FD_multi_effect[!(FD_multi_effect$GRIDREF=="SU0451"),] ## 199 sites now
## multivariate FD metrics for 39 species at 199 sites using effect traits using 3 PC axes
## site SD0451 only has one species (skylark) 
write.csv(FD_multi_effect, file="../Data/Analysis_data/Seed dispersal/FD_multi_seed_effect_det2.csv", row.names=FALSE)
## 180 sites for det1

## check relationship between FDis effect PCA vs FDis effect raw (i.e. standardised trait data)
FDis_effect_PCA <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_effect_4axes.csv", header=TRUE)
FDis_effect_raw <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_effect_raw.csv", header=TRUE)

## remove other columns
FDis_effect_PCA <- FDis_effect_PCA[,-c(2,4:5)]
FDis_effect_raw <- FDis_effect_raw[,-c(2,4:5)]
colnames(FDis_effect_PCA)[2] <- "FDis_PCA" 
colnames(FDis_effect_raw)[2] <- "FDis_raw" 
FDis_effect <- merge(FDis_effect_PCA, FDis_effect_raw, by="GRIDREF")

seed_effect <- ggplot(data = FDis_effect, aes(x = FDis_PCA, y = FDis_raw)) + 
  geom_point(color='black') +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  labs(x = "FDis PCA", y = "FDis standardised traits") +
  ggtitle("Seed dispersal") +
  theme_classic() +
  theme(text = element_text(size = 16))
seed_effect
ggsave(filename="../Graphs/seed_effect_FDis.png", plot=seed_effect, width=10, height=6)

## relationship is linear as predicted

#
#
#
#
#
#
#

######################################################################
####################### RESPONSE TRAITS ##############################
######################################################################

## select response only traits from effect and response file
seed_response <- seed_traits[,c(1,5,16:21)] ## keep name and response traits

## centre and scale trait data
seed_response[c(2:8)] <- lapply(seed_response[c(2:8)], function(seed_response) c(scale(seed_response, center = TRUE, scale = TRUE))) 

## correlation matrix 
response_corr <- cor(seed_response[c(2:8)])
response_corr
## not very correlated, highest r = -0.66 (thermal max and mean latitude), next highest 0.61 (thermal max and brood size)

################################################
## PCA of all response traits
response_pca <- prcomp(seed_response[ ,2:8])
summary(response_pca) ## PC1 accounts for 37.7% of the total variation and PC2 20.5% of the total variation, PCA3 17.3%

response_eigen <- get_eigenvalue(response_pca)
response_eigen ## PCA1 eigenvalue = 2.63, PCA2 eigenvalue = 1.44, PCA3 eigenvalue = 1.21 (keep all 3 as they are above 1)

response_pca$rotation

x <- predict(response_pca) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe

## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2,3)]) ## not correlated
## use these as two independent axes of response traits

seed_response <- cbind(seed_response, x) ## use PCA1, PCA2 and PCA3 to calculate response trait diversity
## correlation between PCA1, PCA2 and PCA3 with response traits
x_corr2 <- cor(seed_response[c(2:11)]) 
x_corr2 

seed_response <- seed_response[,c(1,9:12)] ## keep only name, PCA1, PCA2, PCA3 and PCA4 columns
## change row numbers to species names
rownames(seed_response) <- seed_response[, 1]
seed_response <- seed_response[, -1]

## plot 3D plot 
library(plot3D)
x <- seed_response$PC1
y <- seed_response$PC2
z <- seed_response$PC3
species <- row.names(seed_response)
## save
png("../Graphs/seed_response_traits.png", width=60, height=60, units="mm", res=600)
scatter3D(x,y,z,box=TRUE,pch=16,cex=0.5,d=2,col="black",bty="b2",axes=TRUE,label=species, nticks=10,theta=40, 
          phi=30, xlab="", ylab="", zlab="",par(mar = c(0,0,0,0),
                                                par(oma = c(0,0,0,0))))
dev.off()
text3D(x,y,z, labels = species, add=TRUE)
## x = Mean latitude vs brood size and thermal maximum (37%)
## y = STI vs SSI and clutch size  (21%)
## z = Maximum longevity vs STI (17%)

### use BBS_av_abund created above in effect traits to weight FD metrics

# ################################################
# ########## Total Branch Length FD ##############
# ################################################
# 
# ################################################
# ## First run this on species with only response traits (n=63)
# 
# # find distance matrix 
# d <- dist(as.matrix(seed_response))
# # apply hirarchical clustering
# cluster <- hclust(d, method="average")
# # plot the dendrogram
# plot(cluster)
# # calculate total branch length for each site
# FD <- treedive(BBS_av_abund, cluster, match.force=TRUE) ## works!!!
# # dataframe with FD values and gridrefs (200 in total) for 67 species
# functional_div_response <- data.frame(GRIDREF=site_match$GRIDREF, FD_response=FD)
# ## functional diversity total branch length 
# ## 40 species at 200 sites
# 
# ## save file
# write.csv(functional_div_response, file="../Data/Analysis_data/Seed dispersal/FD_seed_response.csv", row.names=FALSE)

################################################
########## Multivariate FD metrics #############
################################################

## use function dbFD to calculate multivariate FD metrics
FD_multi_response <- dbFD(seed_response, BBS_av_abund)
FD_multi_response <- data.frame(GRIDREF=site_match$GRIDREF, FDVe=FD_multi_response$FEve, FDis=FD_multi_response$FDis, 
                                FRic=FD_multi_response$FRic, FDiv=FD_multi_response$FDiv)
FD_multi_response <- FD_multi_response[!(FD_multi_response$GRIDREF=="SU0451"),] ## 199 sites now
## this is calculated on the 4 PCA axes derived from 7 response traits
write.csv(FD_multi_response, file="../Data/Analysis_data/Seed dispersal/FD_multi_seed_response_det2.csv", row.names=FALSE)
## 199 sites and 39 species using 4 PC axes
## 180 sites and 39 species for det1

###########################
## correlation between FDis effect and FDis response

FDis_effect <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_effect_det2.csv", header=TRUE)
FDis_response <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_response_det2.csv", header=TRUE)

FDis_effect <- FDis_effect[,-c(2,4:5)]
FDis_response <- FDis_response[,-c(2,4:5)]
colnames(FDis_effect)[2] <- "FDis_effect" 
colnames(FDis_response)[2] <- "FDis_response" 
FDis_seed <- merge(FDis_effect, FDis_response, by="GRIDREF")

cor.test(FDis_seed$FDis_effect, FDis_seed$FDis_response) ## p<0.001, r = -0.33

## compare this to random sample of communities
## find species richness for each site and each year
BBS_info <- BBS_data[,c(1,2,4)]
BBS_info <- unique(BBS_info)
BBS_spp_richness2 <- BBS_info %>% group_by(GRIDREF, YEAR) %>% summarise(no_spp = n()) ## number of species at each site each year
BBS_spp_richness <- BBS_spp_richness2 %>% group_by(GRIDREF) %>% summarise(spp_rich = mean(no_spp)) ## average spp richness over time (same for each year)
## use this to randomly sample different communities of same species richness x 100 for now
## use average abundance across all 200 sites for 'fake' communities (from BBS_av_abund file)

## create lookup i.e average abundance over time and sites to use as abundance for 'fake' communities
BBS_abund <- BBS_data %>% group_by(ENGLISH_NAME) %>% summarise(mean_abund=mean(TOT_det))

## first loop: loop through each row in BBS_abund to get the species richness value
## second loop will be i in 100 - i.e. 100 random simulations
sample_final <- NULL

for(j in 1:nrow(BBS_spp_richness)){
  print(j)
  for (i in 1:100){ ## repeat random sampling 10 times
    sp_rich <- BBS_spp_richness$spp_rich[j]
    gridref <- BBS_spp_richness$GRIDREF[j]
    
    sample <- BBS_abund %>% sample_n(sp_rich, replace=FALSE)
    sample_temp <- data.frame(sample$ENGLISH_NAME, sample$mean_abund, gridref, sp_rich,i)
    sample_final <- as.data.frame(rbind(sample_final, sample_temp))
    
    
  }
}

############################ RESPONSE FDis #################################
FDis_response_sample <- NULL
## loop through 1:100 (each simulation)
for (i in 1:100){
  print(i)
  BBS_subset <- sample_final[sample_final$i==i,]
  BBS_subset <- BBS_subset[,c(1:3)]
  BBS_sample <- BBS_subset %>%
    group_by(gridref, sample.ENGLISH_NAME) %>% 
    spread(sample.ENGLISH_NAME, sample.mean_abund, drop = FALSE, fill = 0)
  site_match <- BBS_sample[,(1)]
  ## add in numbers 1:199
  site_match$site <- 1:180
  ## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
  BBS_sample <- BBS_sample[,-(1)]
  
  FD_multi_response <- dbFD(seed_response, BBS_sample, calc.FRic = FALSE, calc.FDiv = FALSE)
  FD_multi_response <- data.frame(GRIDREF=site_match$gridref, FDis=FD_multi_response$FDis, i)
  FDis_response_sample <- rbind(FDis_response_sample, FD_multi_response)
  
}

colnames(FDis_response_sample)[2] <- "FDis_response"
## remove site SU4051 where species richness == 1 so FD cannot be calculated 
FDis_response_sample <- FDis_response_sample[!(FDis_response_sample$GRIDREF=="SU0451"),] ## 199 sites now
## save file 
write.csv(FDis_response_sample, file="../Data/Analysis_data/Seed dispersal/FDis_response_sample_det1.csv", row.names=FALSE)

############################ EFFECT FDis #################################
FDis_effect_sample <- NULL
## loop through 1:100 (each repitition)
for (i in 1:100){
  print(i)
  BBS_subset <- sample_final[sample_final$i==i,]
  BBS_subset <- BBS_subset[,c(1:3)]
  BBS_sample <- BBS_subset %>%
    group_by(gridref, sample.ENGLISH_NAME) %>% 
    spread(sample.ENGLISH_NAME, sample.mean_abund, drop = FALSE, fill = 0)
  site_match <- BBS_sample[,(1)]
  ## add in numbers 1:199
  site_match$site <- 1:180
  ## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
  BBS_sample <- BBS_sample[,-(1)]
  
  FD_multi_effect <- dbFD(seed_effect, BBS_sample, calc.FRic = FALSE, calc.FDiv = FALSE)
  FD_multi_effect <- data.frame(GRIDREF=site_match$gridref, FDis=FD_multi_effect$FDis, i)
  FDis_effect_sample <- rbind(FDis_effect_sample, FD_multi_effect)
  
}
colnames(FDis_effect_sample)[2] <- "FDis_effect"
## remove site SU4051 where species richness == 1 so FD cannot be calculated 
FDis_effect_sample <- FDis_effect_sample[!(FDis_effect_sample$GRIDREF=="SU0451"),] ## 199 sites now
write.csv(FDis_effect_sample, file="../Data/Analysis_data/Seed dispersal/FDis_effect_sample_det1.csv", row.names=FALSE)

FDis_effect_sample <- read.csv("../Data/Analysis_data/Seed dispersal/FDis_effect_sample_det2.csv", header=TRUE)
FDis_response_sample <- read.csv("../Data/Analysis_data/Seed dispersal/FDis_response_sample_det2.csv", header=TRUE)

FDis_effect$i <- 0
FDis_response$i <- 0

FDis_effect_all <- rbind(FDis_effect, FDis_effect_sample)
FDis_response_all <- rbind(FDis_response, FDis_response_sample)

FDis_seed_sample <- merge(FDis_effect_all, FDis_response_all, by=c("GRIDREF", "i"))
## correlation for each repetition (0 = true data)
FDis_seed_cor <- FDis_seed_sample %>%
  group_by(i) %>%
  summarize(COR=cor(FDis_effect, FDis_response))
## take out true correlation(r=-0.34)
FDis_seed_cor2 <- FDis_seed_cor[!(FDis_seed_cor$i==0),]
mean_fake_cor <- mean(FDis_seed_cor2$COR)
mean_fake_cor ## -0.06

## plot true vs simulated correlations
true <- filter(FDis_seed_sample, i==0)
sim <- FDis_seed_sample[!(FDis_seed_sample$i==0),]

## plot
seed_cor <- ggplot(sim, aes(x=FDis_effect, y=FDis_response, group=i)) +
  geom_ribbon(stat='smooth', method = "lm", se=TRUE, alpha=0.1, 
              fill = "gray60") +
  geom_line(stat='smooth', method = "lm", alpha=0.8, lwd=0.3, se=FALSE, colour="black") +
  geom_ribbon(data=true, stat='smooth', method = "lm",
              se=TRUE, alpha=0.3, fill = "gray60") +
  geom_line(data=true, stat='smooth', method = "lm",
            alpha=1.5, lwd=0.3,se=FALSE, colour="red") +
  labs(x = expression("F"[DIS]*" effect traits"), y = expression("F"[DIS]*" response traits")) +
  scale_x_continuous(breaks=seq(0,3,0.5)) +
  theme_classic() +
  theme(text = element_text(size = 6), legend.position="none")
seed_cor
ggsave(filename="../Graphs/seed_FDis_correlation_det1.png", plot=seed_cor, width = 150, 
       height = 100, dpi = 600, units = "mm", device='png')


# #############################################################################
# ####################### BOTH AND EFFECT TRAITS ##############################
# #############################################################################
# 
# ## now run the analysis on both and effect traits (10 traits in total)
# ## run a PCA on both traits, then combine the both PCA and effect PCAs (3 PCAs in total)
# ## and use these 4 PCA axes to calculate FD metrics
# 
# ## select both only traits from all traits file
# seed_both <- seed_traits[,c(1,6,10:14)] ## keep name and both traits (42spp)
# 
# ## first check each trait for normal distribution
# hist(seed_both$Body_Mass) ## right skew
# seed_both$Body_Mass <- log10(seed_both$Body_Mass) ## bit better, still right skew
# hist(seed_both$Tarsus_Length) ## right skew
# seed_both$Tarsus_Length <- log10(seed_both$Tarsus_Length)
# hist(seed_both$Kipp.s_Distance)
# seed_both$Kipp.s_Distance <- log(seed_both$Kipp.s_Distance)
# hist(seed_both$Wing_Chord)
# seed_both$Wing_Chord <- log10(seed_both$Wing_Chord)
# hist(seed_both$Hand.Wing.Index..Claramunt.2011.)
# seed_both$Hand.Wing.Index..Claramunt.2011. <- log(seed_both$Hand.Wing.Index..Claramunt.2011.)
# hist(seed_both$Tail_Length)
# seed_both$Tail_Length <- log(seed_both$Tail_Length)
# 
# ## centre and scale trait data
# seed_both[c(2:7)] <- lapply(seed_both[c(2:7)], function(seed_both) c(scale(seed_both, center = TRUE, scale = TRUE))) 
# 
# ## correlation matrix 
# both_corr <- cor(seed_both[c(2:7)])
# both_corr
# ## high correlations (range from 0.38 to 0.97)
# 
# ################################################
# ## PCA of all both traits
# both_pca <- prcomp(seed_both[ ,2:7]) ## all response traits
# summary(both_pca) ## PC1 accounts for 77% of the total variation and PC2 14% of the total variation
# 
# both_eigen <- get_eigenvalue(both_pca)
# both_eigen ## PCA1 eigenvalue = 4.6, PCA2 eigenvalue = 0.83 (just take first axis - never calc FD for both only)
# 
# x <- predict(both_pca) ## save scores from PCA
# x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
# 
# seed_both <- cbind(seed_both, x) 
# ## correlation between PCA1 and PCA2 with effect traits
# x_corr2 <- cor(seed_both[c(2:13)]) 
# x_corr2 
# ## PCA1 high correlations with all traits
# seed_both <- seed_both[,c(1,8)] ## keep only name and PCA1
# 
# ################################################
# ########## Total Branch Length FD ##############
# ################################################
# 
# ################################################
# ## Run this on effect and both traits
# 
# ## first cbind effect PCA and both PCA axes together
# seed_effect_both <- cbind(seed_effect, seed_both)
# ## remove duplicated english name column
# seed_effect_both <- seed_effect_both[,-3] ## PC1.1 is both traits PCA axis
# 
# # find distance matrix 
# d <- dist(as.matrix(seed_effect_both))
# # apply hirarchical clustering
# cluster <- hclust(d, method="average")
# # plot the dendrogram
# plot(cluster)
# # calculate total branch length for each site
# FD <- treedive(BBS_av_abund, cluster, match.force=TRUE) ## works!!!
# functional_div_effect_both <- data.frame(GRIDREF=site_match$GRIDREF, FD_effect_both=FD)
# ## functional diversity total branch length 
# ## 40 species at 200 sites
# ## using effect and both trait PCAs
# 
# ## save file
# write.csv(functional_div_effect_both, file="../Data/Analysis_data/Seed dispersal/FD_seed_effect_both.csv", row.names=FALSE)
# 
# ################################################
# ########## Multivariate FD metrics #############
# ################################################
# 
# ## First run this on effect and both traits
# FD_multi_effect_both <- dbFD(seed_effect_both, BBS_av_abund)
# FD_multi_effect_both <- data.frame(GRIDREF=site_match$GRIDREF, FDVe=FD_multi_effect_both$FEve, FDis=FD_multi_effect_both$FDis, 
#                                    FRic=FD_multi_effect_both$FRic, FDiv=FD_multi_effect_both$FDiv)
# write.csv(FD_multi_effect_both, file="../Data/Analysis_data/Seed dispersal/FD_multi_seed_effect_both.csv", row.names=FALSE)
# ## 199 sites and 40 species using 3 PCA axes (two from effect, one from both traits)
# 
# ###############################################################################
# ####################### BOTH AND RESPONSE TRAITS ##############################
# ###############################################################################
# 
# ################################################
# ########## Total Branch Length FD ##############
# ################################################
# 
# ################################################
# ## Run this on response and both traits
# 
# ## first cbind response PCA and both PCA axes together
# seed_response_both <- cbind(seed_response, seed_both)
# ## remove duplicated english name column
# seed_response_both <- seed_response_both[,-4] ## PC1.1 is both traits PCA axis
# 
# # find distance matrix 
# d <- dist(as.matrix(seed_response_both))
# # apply hirarchical clustering
# cluster <- hclust(d, method="average")
# # plot the dendrogram
# plot(cluster)
# # calculate total branch length for each site
# FD <- treedive(BBS_av_abund, cluster, match.force=TRUE) ## works!!!
# # dataframe with FD values and gridrefs (200 in total) for 67 species
# functional_div_response_both <- data.frame(GRIDREF=site_match$GRIDREF, FD_response_both=FD)
# ## functional diversity total branch length 
# ## 40 species at 200 sites
# ## using response and both trait PCA axes
# ## save file
# write.csv(functional_div_response_both, file="../Data/Analysis_data/Seed dispersal/FD_seed_response_both.csv", row.names=FALSE)
# 
# ################################################
# ########## Multivariate FD metrics #############
# ################################################
# 
# FD_multi_response_both <- dbFD(seed_response_both, BBS_av_abund) ## warnings with FRic (but this isn't used, just FDis)
# FD_multi_response_both <- data.frame(GRIDREF=site_match$GRIDREF, FDVe=FD_multi_response_both$FEve, FDis=FD_multi_response_both$FDis, 
#                                      FRic=FD_multi_response_both$FRic, FDiv=FD_multi_response_both$FDiv)
# write.csv(FD_multi_response_both, file="../Data/Analysis_data/Seed dispersal/FD_multi_seed_response_both.csv", row.names=FALSE)
# ## 199 sites and 40 species


#################################################################
####################### ALL TRAITS ##############################
#################################################################

seed_all <- seed_traits[,-c(2:4,14)]
## centre and scale trait data
seed_all[c(2:17)] <- lapply(seed_all[c(2:17)], function(seed_all) c(scale(seed_all, center = TRUE, scale = TRUE))) 
#seed_all <- seed_all[c(1,4,5,6,12,2,13,14,15,16,17,8,9,11,7,3,10)]

## lattice correlation graphs
library(lattice)

pairs(seed_all[2:12], pch=19, lower.panel = NULL)

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

png("../Graphs/seed_trait_cor.png", width=800, height=500, units="mm", res=500)
pairs(seed_all[2:17],
      lower.panel=panel.smooth, upper.panel=panel.cor, 
      pch=20)
dev.off()

## effect and response only
library(GGally)
ggpairs(seed_all[2:12])

## correlation matrix 
all_corr <- cor(seed_all[c(2:18)])
all_corr ## range from 0.03-0.97

#############################################################################################
## High correlation between traits so run PCA on all effect traits
pca_all <- prcomp(seed_all[ ,2:17]) ## all 17 traits
summary(pca_all) ## First 5 axes explain 86% of total variation

all_eigen <- get_eigenvalue(pca_all)
all_eigen ## use first 5 axes - they are all above 1

x <- pca_all$rotation

x <- predict(pca_all) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1:5)]) ## not correlated r = ~0
## use these as two independent axes of all traits

seed_all <- cbind(seed_all, x) 
seed_all <- seed_all[,c(1,18:22)] ## keep only name, PCA1 to PCA5 columns
## change row numbers to species names
rownames(seed_all) <- seed_all[, 1]
seed_all <- seed_all[, -1]

################################################
########## Total Branch Length FD ##############
################################################

# find distance matrix 
d <- dist(as.matrix(seed_all))
# apply hirarchical clustering
cluster <- hclust(d, method="average")
# plot the dendrogram
plot(cluster)
# calculate total branch length for each site
FD <- treedive(BBS_av_abund, cluster, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 67 species
functional_div_all <- data.frame(GRIDREF=site_match$GRIDREF, FD_all=FD)
## functional diversity total branch length 
## 40 species at 200 sites
## using effect, response and both trait PCA axes
## save file
write.csv(functional_div_all, file="../Data/Analysis_data/Seed dispersal/FD_seed_all.csv", row.names=FALSE)

################################################
########## Multivariate FD metrics #############
################################################

FD_multi_all <- dbFD(seed_all, BBS_av_abund) ## warnings with FRic (but this isn't used, just FDis)
FD_multi_all <- data.frame(GRIDREF=site_match$GRIDREF, FDVe=FD_multi_all$FEve, FDis=FD_multi_all$FDis, 
                                FRic=FD_multi_all$FRic, FDiv=FD_multi_all$FDiv)
FD_multi_all <- FD_multi_all[!(FD_multi_all$GRIDREF=="SU0451"),] ## 199 sites now
write.csv(FD_multi_all, file="../Data/Analysis_data/Seed dispersal/FD_multi_seed_all_det2.csv", row.names=FALSE)
## 199 sites and 40 species
## 180 sites for det1

