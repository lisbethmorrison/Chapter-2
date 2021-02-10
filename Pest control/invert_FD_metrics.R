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
#BBS_data <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_invert_det_comp.csv", header=TRUE) ## 80 species + 130 sites
BBS_data <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_invert_det_interpol.csv", header=TRUE) ## 87 species + 200 sites

## remove species from seed_traits which are not in BBS complete case data (7 species)
species <- unique(BBS_data$ENGLISH_NAME)
invert_traits <- invert_traits[invert_traits$ENGLISH_NAME %in% species, ]
## only do this for complete case data!

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
summary(pca_effect) # PC1 accounts for 77% of the total variation and PC2 16% of the total variation
## use these two axes to calculate FD of effect traits

## check eigenvalues to determine the number of principal components to be considered
## An eigenvalue > 1 indicates that PCs account for more variance than accounted by one of the original variables in standardized data
## This is commonly used as a cutoff point for which PCs are retained
## You can also limit the number of component to that number that accounts for a certain fraction of the total variance
## For example, if you are satisfied with 70% of the total variance explained then use the number of components to achieve that

effect_eigen <- get_eigenvalue(pca_effect)
effect_eigen ## PCA1 eigenvalue = 3.1, PCA2 eigenvalue = 0.63
## Technically only need PCA1 to explain maximum amount of variation

x <- predict(pca_effect) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe
## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2)]) ## not correlated r = ~0
## use these as two independent axes of effect traits

invert_effect <- cbind(invert_effect, x) ## use PCA1, PCA2 and PCA3 as two traits to calculate effect trait diversity
## correlation between PCA1 and effect traits
x_corr2 <- cor(invert_effect[c(2:6)]) 
x_corr2 ## strong negative correlation with all effect traits
## correlation between PCA2 and effect traits
x_corr3 <- cor(invert_effect[c(2:5,7)])  #
x_corr3 # r=0.6 with bill length, weak correlations with others
invert_effect <- invert_effect[,c(1,6:8)] ## keep only name and first 3 PCA axes
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
## add in numbers 1:130 // 1:200
site_match$site <- 1:200
## remove first column (which is actual gridrefs) and leave sites as numbers 1-200
BBS_av_abund <- BBS_av_abund[,-(1)]

################################################
########## Multivariate FD metrics #############
################################################

## use function dbFD to calculate multivariate FD metrics
FD_multi_effect <- dbFD(invert_effect, BBS_av_abund)
FD_multi_effect <- data.frame(GRIDREF=site_match$GRIDREF, FDis=FD_multi_effect$FDis)
## multivariate FD metrics for 80 species at 130 sites using effect traits using 3 PC axes [complete case det]
## multivariate FD metrics for 87 species at 200 sites using effect traits using 3 PC axes [interpolate det]
write.csv(FD_multi_effect, file="../Data/Analysis_data/Pest control/FD_multi_invert_effect_det_interpol.csv", row.names=FALSE)


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
summary(response_pca) ## PC1 accounts for 40% of the total variation and PC2 18% of the total variation, PCA3 15%

response_eigen <- get_eigenvalue(response_pca)
response_eigen ## PCA1 eigenvalue = 2.8, PCA2 eigenvalue = 1.3, PCA3 eigenvalue = 1.1 (keep top 4 as they are above 1)

x <- predict(response_pca) ## save scores from PCA
x <- as.data.frame(x, stringsAsFactors = FALSE) ## convert to dataframe

## correlation between PCA1 and PCA2
x_corr <- cor(x[c(1,2,3)]) ## not correlated
## use these as two independent axes of response traits

invert_response <- cbind(invert_response, x) ## use PCA1, PCA2, PCA3 and PCA4 to calculate response trait diversity
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

################################################
########## Multivariate FD metrics #############
################################################

## use function dbFD to calculate multivariate FD metrics
FD_multi_response <- dbFD(invert_response, BBS_av_abund)
FD_multi_response <- data.frame(GRIDREF=site_match$GRIDREF,FDis=FD_multi_response$FDis)
## multivariate FD metrics for 80 species at 130 sites using response traits using 4 PC axes [complete case det]
## multivariate FD metrics for 87 species at 200 sites using response traits using 4 PC axes [interpolate det]
write.csv(FD_multi_response, file="../Data/Analysis_data/Pest control/FD_multi_invert_response_det_interpol.csv", row.names=FALSE)

#################################################################
####################### ALL TRAITS ##############################
#################################################################

invert_all <- invert_traits[,-c(2:4,14)]
## centre and scale trait data
invert_all[c(2:17)] <- lapply(invert_all[c(2:17)], function(invert_all) c(scale(invert_all, center = TRUE, scale = TRUE))) 

#############################################################################################
## High correlation between traits so run PCA on all effect traits
pca_all <- prcomp(invert_all[ ,2:17]) ## all 16 traits
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
########## Multivariate FD metrics #############
################################################

FD_multi_all <- dbFD(invert_all, BBS_av_abund) ## warnings with FRic (but this isn't used, just FDis)
FD_multi_all <- data.frame(GRIDREF=site_match$GRIDREF, FDis=FD_multi_all$FDis)
## multivariate FD metrics for 80 species at 130 sites using all traits using 5 PC axes [complete case det]
## multivariate FD metrics for 87 species at 200 sites using all traits using 5 PC axes [interpolate det]
write.csv(FD_multi_all, file="../Data/Analysis_data/Pest control/FD_multi_invert_all_det_interpol.csv", row.names=FALSE)





###########################
## correlation between FDis effect and FDis response

FDis_effect <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_effect_det_interpol.csv", header=TRUE)
FDis_response <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_response_det_interpol.csv", header=TRUE)

colnames(FDis_effect)[2] <- "FDis_effect" 
colnames(FDis_response)[2] <- "FDis_response" 
FDis_invert <- merge(FDis_effect, FDis_response, by="GRIDREF")

cor.test(FDis_invert$FDis_effect, FDis_invert$FDis_response) ## p=0.08, r = -0.12

## compare this to random sample of communities
## find species richness for each site and each year
BBS_info <- BBS_data[,c(1:3)]
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
  
  FD_multi_response <- dbFD(invert_response, BBS_sample, calc.FRic = FALSE, calc.FDiv = FALSE)
  FD_multi_response <- data.frame(GRIDREF=site_match$gridref, FDis=FD_multi_response$FDis, i)
  FDis_response_sample <- rbind(FDis_response_sample, FD_multi_response)
  
}

colnames(FDis_response_sample)[2] <- "FDis_response"
## save file 
write.csv(FDis_response_sample, file="../Data/Analysis_data/Pest control/FDis_response_sample_det_interpol.csv", row.names=FALSE)

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
  
  FD_multi_effect <- dbFD(invert_effect, BBS_sample, calc.FRic = FALSE, calc.FDiv = FALSE)
  FD_multi_effect <- data.frame(GRIDREF=site_match$gridref, FDis=FD_multi_effect$FDis, i)
  FDis_effect_sample <- rbind(FDis_effect_sample, FD_multi_effect)
  
}
colnames(FDis_effect_sample)[2] <- "FDis_effect"
write.csv(FDis_effect_sample, file="../Data/Analysis_data/Pest control/FDis_effect_sample_det_interpol.csv", row.names=FALSE)

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

