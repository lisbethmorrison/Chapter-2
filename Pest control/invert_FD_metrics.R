###########################################################
## Title: Calculate functional diversity of traits for invertivores
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: August 2019
##########################################################

rm(list=ls()) # clear R
options(scipen=999)

library(factoextra)
library(vegan)
library(reshape2)
library(tidyverse)
library(FD)

## read in trait data 
invert_traits <- read.csv("../Data/Trait data/Pest control/BBS_invert_trait_data.csv", header=TRUE)
BBS_data <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_invert_interpol_repeat_det2.csv", header=TRUE)

#### for detectability 1 there are only 55 species with data - need to remove these from invert_traits
species <- unique(BBS_data$ENGLISH_NAME)
invert_traits <- invert_traits[invert_traits$ENGLISH_NAME %in% species, ]

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

####################################################################
####################### EFFECT TRAITS ##############################
####################################################################

## select effect only traits from effect and response file
invert_effect <- invert_traits[,c(1,7:9,15)] ## keep name and effect traits

## centre and scale trait data
invert_effect[c(2:5)] <- lapply(invert_effect[c(2:5)], function(invert_effect) c(scale(invert_effect, center = TRUE, scale = TRUE))) 

## correlation matrix 
effect_corr <- cor(invert_effect[c(2:5)])
effect_corr ## all r>0.55

#############################################################################################
## High correlation between traits so run PCA on all effect traits
pca_effect <- prcomp(invert_effect[ ,2:5]) ## bill length, width and depth and gape width
summary(pca_effect) # PC1 accounts for 79.9% of the total variation and PC2 13.5% of the total variation
## use these two axes to calculate FD of effect traits

## check eigenvalues to determine the number of principal components to be considered
## An eigenvalue > 1 indicates that PCs account for more variance than accounted by one of the original variables in standardized data
## This is commonly used as a cutoff point for which PCs are retained
## You can also limit the number of component to that number that accounts for a certain fraction of the total variance
## For example, if you are satisfied with 70% of the total variance explained then use the number of components to achieve that

effect_eigen <- get_eigenvalue(pca_effect)
effect_eigen ## PCA1 eigenvalue = 3.19, PCA2 eigenvalue = 0.54
## Technically only need PCA1 to explain maximum amount of variation

x <- predict(pca_effect) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2)]) ## not correlated r = ~0
## use these as two independent axes of effect traits

invert_effect <- cbind(invert_effect, x) ## use PCA1 and PCA2 as two traits to calculate effect trait diversity
## correlation between PCA1 and effect traits
x_corr2 <- cor(invert_effect[c(2:6)]) 
x_corr2 ## strong negative correlation with all effect traits
## correlation between PCA2 and effect traits
x_corr3 <- cor(invert_effect[c(2:5,7)])  #
x_corr3 # r=0.6 with bill length, weak correlations with others
invert_effect <- invert_effect[,c(1,6:8)] ## keep only name and all 4 PCA axes
## change row numbers to species names
rownames(invert_effect) <- invert_effect[, 1]
invert_effect <- invert_effect[, -1]

## plot 3D plot 
library(plot3D)
x <- invert_effect$PC1
y <- invert_effect$PC2
z <- invert_effect$PC3
species <- row.names(invert_effect)
## save
png("../Graphs/invert_effect_traits.png", width=60, height=60, units="mm", res=600)
scatter3D(x,y,z,box=TRUE,pch=16,cex=0.5,d=2,col="black",bty="b2",axes=TRUE,label=species, nticks=10,theta=40, 
          phi=30, xlab="", ylab="", zlab="",par(mar = c(0,0,0,0),
                                                 par(oma = c(0,0,0,0))))
dev.off()
text3D(x,y,z, labels = species, add=TRUE)
## x = Beak size (92%)
## y = Bill length vs gape width (4%)
## z = Bill width and depth vs gape width (2%)

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

# ################################################
# ########## Total Branch Length FD ##############
# ################################################
# 
# # find distance matrix 
# d <- dist(as.matrix(invert_effect))
# # apply hirarchical clustering
# cluster <- hclust(d, method="average")
# # plot the dendrogram
# plot(cluster)
# # calculate total branch length for each site
# FD <- treedive(BBS_av_abund, cluster, match.force=TRUE) ## works!!!
# # dataframe with FD values and gridrefs (200 in total) for 53 species
# functional_div_effect <- data.frame(GRIDREF=site_match$GRIDREF, FD_effect=FD)
# ## functional diversity total branch length 
# ## 67 species at 200 sites using effect traits
# 
# ## save file
# write.csv(functional_div_effect, file="../Data/Analysis_data/Pest control/FD_invert_effect.csv", row.names=FALSE)

################################################
########## Multivariate FD metrics #############
################################################

## use function dbFD to calculate multivariate FD metrics
FD_multi_effect <- dbFD(invert_effect, BBS_av_abund)
FD_multi_effect <- data.frame(GRIDREF=site_match$GRIDREF, FDVe=FD_multi_effect$FEve, FDis=FD_multi_effect$FDis, 
                              FRic=FD_multi_effect$FRic, FDiv=FD_multi_effect$FDiv)
## multivariate FD metrics for 67 species at 200 sites using effect traits using 3 PC axes
write.csv(FD_multi_effect, file="../Data/Analysis_data/Pest control/FD_multi_invert_effect_det2.csv", row.names=FALSE)
## 180 sites, 55 species for det1
## 200 sites, 59 species for det2

## check relationship between FDis effect PCA vs FDis effect raw (i.e. standardised trait data)
FDis_effect_PCA <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_effect_4axes.csv", header=TRUE)
FDis_effect_raw <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_effect_raw.csv", header=TRUE)

## remove other columns
FDis_effect_PCA <- FDis_effect_PCA[,-c(2,4:5)]
FDis_effect_raw <- FDis_effect_raw[,-c(2,4:5)]
colnames(FDis_effect_PCA)[2] <- "FDis_PCA" 
colnames(FDis_effect_raw)[2] <- "FDis_raw" 
FDis_effect <- merge(FDis_effect_PCA, FDis_effect_raw, by="GRIDREF")

invert_effect <- ggplot(data = FDis_effect, aes(x = FDis_PCA, y = FDis_raw)) + 
  geom_point(color='black') +
  geom_smooth(method = "lm", se = FALSE, color="red") +
  ggtitle("Pest control") +
  labs(x = "FDis PCA", y = "FDis standardised traits") +
  theme_classic() +
  theme(text = element_text(size = 16))
invert_effect
ggsave(filename="../Graphs/invert_effect_FDis.png", plot=invert_effect, width=10, height=6)

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
invert_response <- invert_traits[,c(1,5,16:21)] ## keep name and response traits

## centre and scale trait data
invert_response[c(2:8)] <- lapply(invert_response[c(2:8)], function(invert_response) c(scale(invert_response, center = TRUE, scale = TRUE))) 

## correlation matrix 
response_corr <- cor(invert_response[c(2:8)])
response_corr
## not very correlated, highest r = -0.62 (STI/thermal max and mean latitude), next highest 0.37

################################################
## PCA of all response traits
response_pca <- prcomp(invert_response[ ,2:8]) ## bill length, width and depth and gape width
summary(response_pca) ## PC1 accounts for 35% of the total variation and PC2 21% of the total variation, PCA3 17%

response_eigen <- get_eigenvalue(response_pca)
response_eigen ## PCA1 eigenvalue = 2.46, PCA2 eigenvalue = 1.45, PCA3 eigenvalue = 1.19 (keep all 3 as they are above 1)

x <- predict(response_pca) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe

## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2,3)]) ## not correlated
## use these as two independent axes of response traits

invert_response <- cbind(invert_response, x) ## use PCA1, PCA2 and PCA3 to calculate response trait diversity
## correlation between PCA1, PCA2 and PCA3 with response traits
x_corr2 <- cor(invert_response[c(2:11)]) 
x_corr2 
## PCA1 correlates with mean latitude, STI and thermal max
## PCA2 correlates with SSI and max longevity
## PCA3 correlates with clutch and brood size
invert_response <- invert_response[,c(1,9:12)] ## keep only name, PCA1, PCA2, PCA3 and PCA4 columns
## change row numbers to species names
rownames(invert_response) <- invert_response[, 1]
invert_response <- invert_response[, -1]

## plot 3D plot 
library(plot3D)
x <- invert_response$PC1
y <- invert_response$PC2
z <- invert_response$PC3
species <- row.names(invert_response)
## save
png("../Graphs/invert_response_traits.png", width=60, height=60, units="mm", res=600)
scatter3D(x,y,z,box=TRUE,pch=16,cex=0.5,d=2,col="black",bty="b2",axes=TRUE,label=species, nticks=10,theta=40, 
          phi=30, xlab="", ylab="", zlab="",par(mar = c(0,0,0,0),
                                                par(oma = c(0,0,0,0))))
dev.off()
text3D(x,y,z, labels = species, add=TRUE)
## x = Mean latitude vs STI and thermal maximum (35%)
## y = Clutch size vs maximum longevity (20%)
## z = Clutch size vs brood size (18%)

### use BBS_av_abund created above in effect traits to weight FD metrics

# ################################################
# ########## Total Branch Length FD ##############
# ################################################
# 
# ################################################
# ## First run this on species with only response traits (n=63)
# 
# # find distance matrix 
# d <- dist(as.matrix(invert_response))
# # apply hirarchical clustering
# cluster <- hclust(d, method="average")
# # plot the dendrogram
# plot(cluster)
# # calculate total branch length for each site
# FD <- treedive(BBS_av_abund, cluster, match.force=TRUE) ## works!!!
# # dataframe with FD values and gridrefs (200 in total) for 67 species
# functional_div_response <- data.frame(GRIDREF=site_match$GRIDREF, FD_response=FD)
# ## functional diversity total branch length 
# ## 67 species at 200 sites
# 
# ## save file
# write.csv(functional_div_response, file="../Data/Analysis_data/Pest control/FD_invert_response.csv", row.names=FALSE)

################################################
########## Multivariate FD metrics #############
################################################

## use function dbFD to calculate multivariate FD metrics
FD_multi_response <- dbFD(invert_response, BBS_av_abund)
FD_multi_response <- data.frame(GRIDREF=site_match$GRIDREF, FDVe=FD_multi_response$FEve, FDis=FD_multi_response$FDis, 
                                FRic=FD_multi_response$FRic, FDiv=FD_multi_response$FDiv)
## this is calculated on the 4 PCA axes derived from 7 response traits 
write.csv(FD_multi_response, file="../Data/Analysis_data/Pest control/FD_multi_invert_response_det2.csv", row.names=FALSE)
## all 200 sites and 59 species
## 180 sites and 55 species for det1

###########################
## correlation between FDis effect and FDis response

FDis_effect <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_effect_det2.csv", header=TRUE)
FDis_response <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_response_det2.csv", header=TRUE)

FDis_effect <- FDis_effect[,-c(2,4:5)]
FDis_response <- FDis_response[,-c(2,4:5)]
colnames(FDis_effect)[2] <- "FDis_effect" 
colnames(FDis_response)[2] <- "FDis_response" 
FDis_invert <- merge(FDis_effect, FDis_response, by="GRIDREF")

cor.test(FDis_invert$FDis_effect, FDis_invert$FDis_response) ## p=0.07, r = -0.13

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
## second loop: i in 100
sample_final <- NULL

for(j in 1:nrow(BBS_spp_richness)){
  print(j)
  for (i in 1:100){ ## repeat random sampling 100 times
    sp_rich <- BBS_spp_richness$spp_rich[j]
    gridref <- BBS_spp_richness$GRIDREF[j]

    sample <- BBS_abund %>% sample_n(sp_rich, replace=FALSE)
    sample_temp <- data.frame(sample$ENGLISH_NAME, sample$mean_abund, gridref, sp_rich,i)
    sample_final <- as.data.frame(rbind(sample_final, sample_temp))
    
    
  }
}

############################ RESPONSE FDis #################################
FDis_response_sample <- NULL
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
  site_match$site <- 1:200
  ## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
  BBS_sample <- BBS_sample[,-(1)]
  
  FD_multi_response <- dbFD(invert_response, BBS_sample)
  FD_multi_response <- data.frame(GRIDREF=site_match$gridref, FDis=FD_multi_response$FDis, i)
  FDis_response_sample <- rbind(FDis_response_sample, FD_multi_response)
  
}

colnames(FDis_response_sample)[2] <- "FDis_response"
## save file 
write.csv(FDis_response_sample, file="../Data/Analysis_data/Pest control/FDis_response_sample_det2.csv", row.names=FALSE)

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
  site_match$site <- 1:200
  ## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
  BBS_sample <- BBS_sample[,-(1)]
  
  FD_multi_effect <- dbFD(invert_effect, BBS_sample)
  FD_multi_effect <- data.frame(GRIDREF=site_match$gridref, FDis=FD_multi_effect$FDis, i)
  FDis_effect_sample <- rbind(FDis_effect_sample, FD_multi_effect)
  
}
colnames(FDis_effect_sample)[2] <- "FDis_effect"
write.csv(FDis_effect_sample, file="../Data/Analysis_data/Pest control/FDis_effect_sample_det2.csv", row.names=FALSE)

FDis_effect_sample <- read.csv("../Data/Analysis_data/Pest control/FDis_effect_sample_det2.csv", header=TRUE)
FDis_response_sample <- read.csv("../Data/Analysis_data/Pest control/FDis_response_sample_det2.csv", header=TRUE)

FDis_effect$i <- 0
FDis_response$i <- 0

FDis_effect_all <- rbind(FDis_effect, FDis_effect_sample)
FDis_response_all <- rbind(FDis_response, FDis_response_sample)

FDis_invert_sample <- merge(FDis_effect_all, FDis_response_all, by=c("GRIDREF", "i"))
## correlation for each repetition (0 = true data)
FDis_invert_cor <- FDis_invert_sample %>%
  group_by(i) %>%
  summarize(COR=cor(FDis_effect, FDis_response))
## true community has negative correlation
## all 'fake' communities have positive correlations 
## weak negative correlation between FDis effect and FDis response must be due to particular communities we have observed
## take out true correlation(r=-0.12)
FDis_invert_cor2 <- FDis_invert_cor[!(FDis_invert_cor$i==0),]
mean_fake_cor <- mean(FDis_invert_cor2$COR)
mean_fake_cor ## 0.34

true <- filter(FDis_invert_sample, i==0)
sim <- FDis_invert_sample[!(FDis_invert_sample$i==0),]

## plot
invert_cor <- ggplot(sim, aes(x=FDis_effect, y=FDis_response, group=i)) +
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
invert_cor
ggsave(filename="../Graphs/invert_FDis_correlation_det2.png", plot=invert_cor, width=10, height=6)

## plot both seed and invert together
library(cowplot)
fig3 <- plot_grid(seed_cor, invert_cor, labels=c("(a)", "(b)"), label_size=6, ncol = 2, nrow = 1, hjust=0)
fig3
ggsave(filename="../Graphs/fig3_det2.png", plot=fig3, width = 80, 
       height = 40, dpi = 600, units = "mm", device='png')

# #############################################################################
# ####################### BOTH AND EFFECT TRAITS ##############################
# #############################################################################
# 
# ## now run the analysis on both and effect traits (10 traits in total)
# ## run a PCA on both traits, then combine the both PCA and effect PCAs (3 PCAs in total)
# ## and use these 4 PCA axes to calculate FD metrics
# 
# ## select both only traits from all traits file
# invert_both <- invert_traits[,c(1,6,10:14)] ## keep name and both traits (42spp)
# 
# ## first check each trait for normal distribution
# hist(invert_both$Body_Mass) ## right skew
# invert_both$Body_Mass <- log10(invert_both$Body_Mass) ## better
# hist(invert_both$Tarsus_Length) ## right skew
# invert_both$Tarsus_Length <- log(invert_both$Tarsus_Length)
# hist(invert_both$Kipp.s_Distance)
# invert_both$Kipp.s_Distance <- log(invert_both$Kipp.s_Distance)
# hist(invert_both$Wing_Chord)
# invert_both$Wing_Chord <- log10(invert_both$Wing_Chord)
# hist(invert_both$Hand.Wing.Index..Claramunt.2011.)
# invert_both$Hand.Wing.Index..Claramunt.2011. <- log(invert_both$Hand.Wing.Index..Claramunt.2011.)
# hist(invert_both$Tail_Length)
# invert_both$Tail_Length <- log(invert_both$Tail_Length)
# 
# ## centre and scale trait data
# invert_both[c(2:7)] <- lapply(invert_both[c(2:7)], function(invert_both) c(scale(invert_both, center = TRUE, scale = TRUE))) 
# 
# ## correlation matrix 
# both_corr <- cor(invert_both[c(2:7)])
# both_corr
# ## high correlations (range from 0.38 to 0.94)
# 
# ################################################
# ## PCA of all both traits
# both_pca <- prcomp(invert_both[ ,2:7]) ## all response traits
# summary(both_pca) ## PC1 accounts for 72% of the total variation and PC2 11% of the total variation
# 
# both_eigen <- get_eigenvalue(both_pca)
# both_eigen ## PCA1 eigenvalue = 4.36, PCA2 eigenvalue = 0.69 (just take first axis - never calc FD for both only)
# 
# x <- predict(both_pca) ## save scores from PCA
# x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
# 
# invert_both <- cbind(invert_both, x) 
# ## correlation between PCA1 and PCA2 with effect traits
# x_corr2 <- cor(invert_both[c(2:13)]) 
# x_corr2 
# ## PCA1 high correlations with all traits
# invert_both <- invert_both[,c(1,8)] ## keep only name and PCA1
# 
# ################################################
# ########## Total Branch Length FD ##############
# ################################################
# 
# ################################################
# ## Run this on effect and both traits
# 
# ## first cbind effect PCA and both PCA axes together
# invert_effect_both <- cbind(invert_effect, invert_both)
# ## remove duplicated english name column
# invert_effect_both <- invert_effect_both[,-3] ## PC1.1 is both traits PCA axis
# 
# # find distance matrix 
# d <- dist(as.matrix(invert_effect_both))
# # apply hirarchical clustering
# cluster <- hclust(d, method="average")
# # plot the dendrogram
# plot(cluster)
# # calculate total branch length for each site
# FD <- treedive(BBS_av_abund, cluster, match.force=TRUE) ## works!!!
# functional_div_effect_both <- data.frame(GRIDREF=site_match$GRIDREF, FD_effect_both=FD)
# ## functional diversity total branch length 
# ## 67 species at 200 sites
# ## using effect and both trait PCAs
# 
# ## save file
# write.csv(functional_div_effect_both, file="../Data/Analysis_data/Pest control/FD_invert_effect_both.csv", row.names=FALSE)
# 
# ################################################
# ########## Multivariate FD metrics #############
# ################################################
# 
# ## First run this on effect and both traits
# FD_multi_effect_both <- dbFD(invert_effect_both, BBS_av_abund)
# FD_multi_effect_both <- data.frame(GRIDREF=site_match$GRIDREF, FDVe=FD_multi_effect_both$FEve, FDis=FD_multi_effect_both$FDis, 
#                                    FRic=FD_multi_effect_both$FRic, FDiv=FD_multi_effect_both$FDiv)
# write.csv(FD_multi_effect_both, file="../Data/Analysis_data/Pest control/FD_multi_invert_effect_both.csv", row.names=FALSE)
# ## all 200 sites and 67 species using 3 PCA axes (two from effect, one from both traits)
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
# invert_response_both <- cbind(invert_response, invert_both)
# ## remove duplicated english name column
# invert_response_both <- invert_response_both[,-4] ## PC1.1 is both traits PCA axis
# 
# # find distance matrix 
# d <- dist(as.matrix(invert_response_both))
# # apply hirarchical clustering
# cluster <- hclust(d, method="average")
# # plot the dendrogram
# plot(cluster)
# # calculate total branch length for each site
# FD <- treedive(BBS_av_abund, cluster, match.force=TRUE) ## works!!!
# # dataframe with FD values and gridrefs (200 in total) for 67 species
# functional_div_response_both <- data.frame(GRIDREF=site_match$GRIDREF, FD_response_both=FD)
# ## functional diversity total branch length 
# ## 67 species at 200 sites
# ## using response and both trait PCA axes
# ## save file
# write.csv(functional_div_response_both, file="../Data/Analysis_data/Pest control/FD_invert_response_both.csv", row.names=FALSE)
# 
# ################################################
# ########## Multivariate FD metrics #############
# ################################################
# 
# FD_multi_response_both <- dbFD(invert_response_both, BBS_av_abund) ## warnings with FRic (but this isn't used, just FDis)
# FD_multi_response_both <- data.frame(GRIDREF=site_match$GRIDREF, FDVe=FD_multi_response_both$FEve, FDis=FD_multi_response_both$FDis, 
#                                      FRic=FD_multi_response_both$FRic, FDiv=FD_multi_response_both$FDiv)
# write.csv(FD_multi_response_both, file="../Data/Analysis_data/Pest control/FD_multi_invert_response_both.csv", row.names=FALSE)
# ## all 200 sites and 67 species


#################################################################
####################### ALL TRAITS ##############################
#################################################################

invert_all <- invert_traits[,-c(2:4,14)]
## centre and scale trait data
invert_all[c(2:17)] <- lapply(invert_all[c(2:17)], function(invert_all) c(scale(invert_all, center = TRUE, scale = TRUE))) 
#invert_all <- invert_all[c(1,4,5,6,12,2,13,14,15,16,17,18,8,9,11,7,3,10)]

library(lattice)

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


png("../Graphs/invert_trait_cor.png", width=800, height=500, units="mm", res=500)
pairs(invert_all[2:18],
      lower.panel=panel.smooth, upper.panel=panel.cor, 
      pch=20)
dev.off()

## correlation matrix 
all_corr <- cor(invert_all[c(2:18)])
all_corr ## range from 0.009-0.9

#############################################################################################
## High correlation between traits so run PCA on all effect traits
pca_all <- prcomp(invert_all[ ,2:17]) ## all 17 traits
summary(pca_all) ## First 5 axes explain 87% of total variation

all_eigen <- get_eigenvalue(pca_all)
all_eigen ## use first 5 axes - they are all above 1

x <- pca_all$rotation
x <- predict(pca_all) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1:5)]) ## not correlated r = ~0
## use these as two independent axes of all traits

invert_all <- cbind(invert_all, x) ## use PCA1 and PCA2 as two traits to calculate effect trait diversity
invert_all <- invert_all[,c(1,18:22)] ## keep only name, PCA1 to PCA5 columns
## change row numbers to species names
rownames(invert_all) <- invert_all[, 1]
invert_all <- invert_all[, -1]

################################################
########## Total Branch Length FD ##############
################################################

################################################
## Run this on effect, response and both traits

# find distance matrix 
d <- dist(as.matrix(invert_all))
# apply hirarchical clustering
cluster <- hclust(d, method="average")
# plot the dendrogram
plot(cluster)
# calculate total branch length for each site
FD <- treedive(BBS_av_abund, cluster, match.force=TRUE) ## works!!!
# dataframe with FD values and gridrefs (200 in total) for 67 species
functional_div_all <- data.frame(GRIDREF=site_match$GRIDREF, FD_all=FD)
## functional diversity total branch length 
## 67 species at 200 sites
## using effect, response and both trait PCA axes
## save file
write.csv(functional_div_all, file="../Data/Analysis_data/Pest control/FD_invert_all.csv", row.names=FALSE)

################################################
########## Multivariate FD metrics #############
################################################

FD_multi_all <- dbFD(invert_all, BBS_av_abund) ## warnings with FRic (but this isn't used, just FDis)
FD_multi_all <- data.frame(GRIDREF=site_match$GRIDREF, FDVe=FD_multi_all$FEve, FDis=FD_multi_all$FDis, 
                                     FRic=FD_multi_all$FRic, FDiv=FD_multi_all$FDiv)
write.csv(FD_multi_all, file="../Data/Analysis_data/Pest control/FD_multi_invert_all_det2.csv", row.names=FALSE)
## all 200 sites and 59 species


