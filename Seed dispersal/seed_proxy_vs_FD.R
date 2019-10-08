###########################################################
## Title: Analyse relationship between EF proxy and FD metrics
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: July 2019
##########################################################

rm(list=ls()) # clear R
options(scipen=999)

library(rcompanion)
library(ggplot2)

## read in data
proxy_effect <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_effect.csv", header=TRUE) ## 32 spp 199 sites
proxy_response <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_response.csv", header=TRUE) ## 38 spp 200 sites
proxy_all <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_all.csv", header=TRUE) ## 28 spp 199 sites

## FD dendrogram results
FD_effect <- read.csv("../Data/Analysis_data/Seed dispersal/FD_effect_32spp.csv", header=TRUE) ## full effect traits
FD_effect2 <- read.csv("../Data/Analysis_data/Seed dispersal/FD_effect_28spp.csv", header=TRUE) ## subset effect traits
FD_response <- read.csv("../Data/Analysis_data/Seed dispersal/FD_response_38spp.csv", header=TRUE) ## full response traits
FD_response2 <- read.csv("../Data/Analysis_data/Seed dispersal/FD_response_28spp.csv", header=TRUE) ## subset response traits
FD_effect_both <- read.csv("../Data/Analysis_data/Seed dispersal/FD_effect_both_32spp.csv", header=TRUE) ## full effect and both traits
FD_effect_both2 <- read.csv("../Data/Analysis_data/Seed dispersal/FD_effect_both_28spp.csv", header=TRUE) ## subset effect and both traits
FD_response_both <- read.csv("../Data/Analysis_data/Seed dispersal/FD_response_both_38spp.csv", header=TRUE) ## full response and both traits
FD_response_both2 <- read.csv("../Data/Analysis_data/Seed dispersal/FD_response_both_28spp.csv", header=TRUE) ## subset response and both traits

## multivariate FD results
FD_multi_effect <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_effect_32spp.csv", header=TRUE) ## full effect traits
FD_multi_effect2 <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_effect_28spp.csv", header=TRUE) ## subset effect traits
FD_multi_response <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_response_38spp.csv", header=TRUE) ## full response traits
FD_multi_response2 <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_response_28spp.csv", header=TRUE) ## subset response traits
FD_multi_effect_both <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_effect_both_32spp.csv", header=TRUE) ## full effect and both traits
FD_multi_effect_both2 <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_effect_both_28spp.csv", header=TRUE) ## subset effect and both traits
FD_multi_response_both <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_response_both_38spp.csv", header=TRUE) ## full response and both traits
FD_multi_response_both2 <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_response_both_28spp.csv", header=TRUE) ## subset response and both traits

##################################################################################################################################################################

##################### DATA COLLATION FOR ANALYSIS ##################### 

## effect ONLY species and sites (n=32)
## merge FD dendro and multivariate together
effect_32 <- merge(FD_effect, FD_multi_effect, by="gridref")
## merge with proxy for 32 species
effect_32 <- merge(effect_32, proxy_effect, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(effect_32)[2] <- "FD"
## save file
write.csv(effect_32, file="../Data/Analysis_data/Seed dispersal/Effect_32spp_final.csv", row.names=FALSE)
## 32 species at 199 sites

## response ONLY species and sites (n=38)
## merge FD dendro and multivariate together
response_38 <- merge(FD_response, FD_multi_response, by="gridref")
## merge with proxy for 32 species
response_38 <- merge(response_38, proxy_response, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(response_38)[2] <- "FD"
## save file
write.csv(response_38, file="../Data/Analysis_data/Seed dispersal/Response_38spp_final.csv", row.names=FALSE)
## 38 species at 200 sites

## effect and both only species and sites (n=32)
## merge FD dendro and multivariate together
effect_both_32 <- merge(FD_effect_both, FD_multi_effect_both, by="gridref")
## merge with proxy for 32 species
effect_both_32 <- merge(effect_both_32, proxy_effect, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(effect_both_32)[2] <- "FD"
## save file
write.csv(effect_both_32, file="../Data/Analysis_data/Seed dispersal/Effect_both_32spp_final.csv", row.names=FALSE)
## 32 species at 199 sites

## response and both only species and sites (n=38)
## merge FD dendro and multivariate together
response_both_38 <- merge(FD_response_both, FD_multi_response, by="gridref")
## merge with proxy for 32 species
response_both_38 <- merge(response_both_38, proxy_response, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(response_both_38)[2] <- "FD"
## save file
write.csv(response_both_38, file="../Data/Analysis_data/Seed dispersal/Response_both_38spp_final.csv", row.names=FALSE)
## 38 species at 200 sites

############## same as above for subset of species

## effect ONLY species and sites (n=28)
## merge FD dendro and multivariate together
effect_32_sub <- merge(FD_effect2, FD_multi_effect2, by="gridref")
## merge with proxy for 28 species
effect_32_sub <- merge(effect_32_sub, proxy_effect, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(effect_32_sub)[2] <- "FD"
## save file
write.csv(effect_32_sub, file="../Data/Analysis_data/Seed dispersal/Effect_28spp_final.csv", row.names=FALSE)
## 28 species at 199 sites

## response ONLY species and sites (n=28)
## merge FD dendro and multivariate together
response_38_sub <- merge(FD_response2, FD_multi_response2, by="gridref")
## merge with proxy for 28 species
response_38_sub <- merge(response_38_sub, proxy_response, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(response_38_sub)[2] <- "FD"
## save file
write.csv(response_38_sub, file="../Data/Analysis_data/Seed dispersal/Response_28spp_final.csv", row.names=FALSE)
## 28 species at 200 sites

## effect and both only species and sites (n=28)
## merge FD dendro and multivariate together
effect_both_32_sub <- merge(FD_effect_both2, FD_multi_effect_both2, by="gridref")
## merge with proxy for 28 species
effect_both_32_sub <- merge(effect_both_32_sub, proxy_effect, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(effect_both_32_sub)[2] <- "FD"
## save file
write.csv(effect_both_32_sub, file="../Data/Analysis_data/Seed dispersal/Effect_both_28spp_final.csv", row.names=FALSE)
## 28 species at 199 sites

## response and both only species and sites (n=28)
## merge FD dendro and multivariate together
response_both_38_sub <- merge(FD_response_both2, FD_multi_response2, by="gridref")
## merge with proxy for 28 species
response_both_38_sub <- merge(response_both_38_sub, proxy_response, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(response_both_38_sub)[2] <- "FD"
## save file
write.csv(response_both_38_sub, file="../Data/Analysis_data/Seed dispersal/Response_both_28spp_final.csv", row.names=FALSE)
## 28 species at 200 sites

####################################################################################################################################

rm(list=ls()) # clear R

## analysis on full dataset (NOT subset of 28 species)
## read in data
effect_32 <- read.csv("../Data/Analysis_data/Seed dispersal/Effect_32spp_final.csv", header=TRUE)
response_38 <- read.csv("../Data/Analysis_data/Seed dispersal/Response_38spp_final.csv", header=TRUE)
effect_both_32 <- read.csv("../Data/Analysis_data/Seed dispersal/Effect_both_32spp_final.csv", header=TRUE)
response_both_38 <- read.csv("../Data/Analysis_data/Seed dispersal/Response_both_38spp_final.csv", header=TRUE)

## keep only gridref, FD, FDis, mean and stability from each dataset
effect_32 <- effect_32[,c(1,2,4,7:8)]
response_38 <- response_38[,c(1,2,4,7:8)]
effect_both_32 <- effect_both_32[,c(1,2,4,7:8)]
response_both_38 <- response_both_38[,c(1,2,4,7:8)]

#######################################################
#### Hypothesis 1: effect trait diversity correlates with mean function, not stability
str(effect_32)
hist(effect_32$mean) ## right skew
effect_32$mean_log <- log(effect_32$mean)
hist(effect_32$mean_log) ## much better
hist(effect_32$stability) ## slight right skew
effect_32$stability_log <- log(effect_32$stability)
hist(effect_32$stability_log) ## better

## run models
mean_FD <- lm(mean_log ~ FD + FDis, data=effect_32)
summary(mean_FD) ## mean has positive relationship with FD and negative relationship with FDis

par(mfrow=c(1,1))
hist(residuals(mean_FD))
par(mfrow=c(2,2))
plot(mean_FD)

## plot result mean ~ FD
ggplot(effect_32, aes(x = FD, y = mean_log)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black") +
  xlab("FD effect traits") +
  ylab("(log) Mean") +
  theme_classic()

stability_FD <- lm(stability_log ~ FD + FDis, data=effect_32)
summary(stability_FD) ## stability has negative relationship with FD and positive relationship with FDis

## plot result stability ~ FD
ggplot(effect_32, aes(x = FD, y = stability_log)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black") +
  xlab("FD effect traits") +
  ylab("(log) Stability") +
  theme_classic()

par(mfrow=c(1,1))
hist(residuals(stability_FD))
par(mfrow=c(2,2))
plot(stability_FD) 

#######################################################
#### Hypothesis 2: response trait diversity correlates with stability of function, not mean
str(response_38)
hist(response_38$mean) ## right skew
response_38$mean_log <- log(response_38$mean)
hist(response_38$mean_log) ## much better
hist(response_38$stability) ## slight right skew
response_38$stability_log <- log(response_38$stability)
hist(response_38$stability_log) ## better

## run models
mean_FD2 <- lm(mean_log ~ FD + FDis, data=response_38)
summary(mean_FD2) ## mean has positive relationship FDis

par(mfrow=c(1,1))
hist(residuals(mean_FD2))
par(mfrow=c(2,2))
plot(mean_FD2) 

stability_FD2 <- lm(stability_log ~ FD + FDis, data=response_38)
summary(stability_FD2) ## stability has negative relationship with FD and positive relationship with FDis

par(mfrow=c(1,1))
hist(residuals(stability_FD2))
par(mfrow=c(2,2))
plot(stability_FD2) 

#######################################################
#### Hypothesis 3: effect and both trait diversity correlates more strongly with mean compared to stability
str(effect_both_32)
hist(effect_both_32$mean) ## right skew
effect_both_32$mean_log <- log(effect_both_32$mean)
hist(effect_both_32$mean_log) ## much better
hist(effect_both_32$stability) ## slight right skew
effect_both_32$stability_log <- log(effect_both_32$stability)
hist(effect_both_32$stability_log) ## better

## run models
mean_FD3 <- lm(mean_log ~ FD + FDis, data=effect_both_32)
summary(mean_FD3) ## mean has positive relationship FD and negative relationship with FDis

par(mfrow=c(1,1))
hist(residuals(mean_FD3))
par(mfrow=c(2,2))
plot(mean_FD3) 

## plot result mean ~ FD
ggplot(effect_both_32, aes(x = FD, y = mean_log)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black")

stability_FD3 <- lm(stability_log ~ FD + FDis, data=effect_both_32)
summary(stability_FD3) ## stability has negative relationship with FD and positive relationship with FDis

## plot result stability ~ FD
ggplot(effect_both_32, aes(x = FD, y = stability_log)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black")

par(mfrow=c(1,1))
hist(residuals(stability_FD3))
par(mfrow=c(2,2))
plot(stability_FD3) 

#######################################################
#### Hypothesis 4: response and both trait diversity correlates more strongly with stability compared to mean
str(response_both_38)
hist(response_both_38$mean) ## right skew
response_both_38$mean_log <- log(response_both_38$mean)
hist(response_both_38$mean_log) ## much better
hist(response_both_38$stability) ## slight right skew
response_both_38$stability_log <- log(response_both_38$stability)
hist(response_both_38$stability_log) ## better

## run models
mean_FD4 <- lm(mean_log ~ FD + FDis, data=response_both_38)
summary(mean_FD4) ## mean has positive relationship FD and FDis

par(mfrow=c(1,1))
hist(residuals(mean_FD4))
par(mfrow=c(2,2))
plot(mean_FD4) ## row 136 is an outlier - high FDis compared to other sites

stability_FD4 <- lm(stability_log ~ FD + FDis, data=response_both_38)
summary(stability_FD4) ## stability has negative relationship with FD and positive relationship with FDis

par(mfrow=c(1,1))
hist(residuals(mean_FD4))
par(mfrow=c(2,2))
plot(mean_FD4) ## row 136 is an outlier - high FDis compared to other sites


############################################################
######## RUN SAME ANALYSIS ON SUBSET OF SPECIES ############
############################################################

rm(list=ls()) # clear R

## analysis on full dataset (NOT subset of 28 species)
## read in data
effect_28 <- read.csv("../Data/Analysis_data/Seed dispersal/Effect_28spp_final.csv", header=TRUE)
response_28 <- read.csv("../Data/Analysis_data/Seed dispersal/Response_28spp_final.csv", header=TRUE)
effect_both_28 <- read.csv("../Data/Analysis_data/Seed dispersal/Effect_both_28spp_final.csv", header=TRUE)
response_both_28 <- read.csv("../Data/Analysis_data/Seed dispersal/Response_both_28spp_final.csv", header=TRUE)

## keep only gridref, FD, FDis, mean and stability from each dataset
effect_28 <- effect_28[,c(1,2,4,7:8)]
response_28 <- response_28[,c(1,2,4,7:8)]
effect_both_28 <- effect_both_28[,c(1,2,4,7:8)]
response_both_28 <- response_both_28[,c(1,2,4,7:8)]

#######################################################
#### Hypothesis 1: effect trait diversity correlates with mean function, not stability
str(effect_28)
hist(effect_28$mean) ## right skew
effect_28$mean_log <- log(effect_28$mean)
hist(effect_28$mean_log) ## much better
hist(effect_28$stability) ## slight right skew
effect_28$stability_log <- log(effect_28$stability)
hist(effect_28$stability_log) ## better

## run models
mean_FD <- lm(mean_log ~ FD + FDis, data=effect_28)
summary(mean_FD) ## mean has positive relationship with FD and negative relationship with FDis

par(mfrow=c(1,1))
hist(residuals(mean_FD))
par(mfrow=c(2,2))
plot(mean_FD) 

stability_FD <- lm(stability_log ~ FD + FDis, data=effect_28)
summary(stability_FD) ## stability has negative relationship with FD and positive relationship with FDis

par(mfrow=c(1,1))
hist(residuals(stability_FD))
par(mfrow=c(2,2))
plot(stability_FD) 

#######################################################
#### Hypothesis 2: response trait diversity correlates with stability of function, not mean
str(response_28)
hist(response_28$mean) ## right skew
response_28$mean_log <- log(response_28$mean)
hist(response_28$mean_log) ## much better
hist(response_28$stability) ## slight right skew
response_28$stability_log <- log(response_28$stability)
hist(response_28$stability_log) ## better

## run models
mean_FD2 <- lm(mean_log ~ FD + FDis, data=response_28)
summary(mean_FD2) ## mean has positive relationship FD

par(mfrow=c(1,1))
hist(residuals(mean_FD2))
par(mfrow=c(2,2))
plot(mean_FD2) 

stability_FD2 <- lm(stability_log ~ FD + FDis, data=response_28)
summary(stability_FD2) ## stability has negative relationship with FD and positive relationship with FDis

par(mfrow=c(1,1))
hist(residuals(stability_FD2))
par(mfrow=c(2,2))
plot(stability_FD2) 

#######################################################
#### Hypothesis 3: effect and both trait diversity correlates more strongly with mean compared to stability
str(effect_both_28)
hist(effect_both_28$mean) ## right skew
effect_both_28$mean_log <- log(effect_both_28$mean)
hist(effect_both_28$mean_log) ## much better
hist(effect_both_28$stability) ## slight right skew
effect_both_28$stability_log <- log(effect_both_28$stability)
hist(effect_both_28$stability_log) ## better

## run models
mean_FD3 <- lm(mean_log ~ FD + FDis, data=effect_both_28)
summary(mean_FD3) ## mean has positive relationship FD and negative relationship with FDis

par(mfrow=c(1,1))
hist(residuals(mean_FD3))
par(mfrow=c(2,2))
plot(mean_FD3) 

stability_FD3 <- lm(stability_log ~ FD + FDis, data=effect_both_28)
summary(stability_FD3) ## stability has negative relationship with FD and positive relationship with FDis

par(mfrow=c(1,1))
hist(residuals(stability_FD3))
par(mfrow=c(2,2))
plot(stability_FD3) 

#######################################################
#### Hypothesis 4: response and both trait diversity correlates more strongly with stability compared to mean
str(response_both_28)
hist(response_both_28$mean) ## right skew
response_both_28$mean_log <- log(response_both_28$mean)
hist(response_both_28$mean_log) ## much better
hist(response_both_28$stability) ## slight right skew
response_both_28$stability_log <- log(response_both_28$stability)
hist(response_both_28$stability_log) ## better

## run models
mean_FD4 <- lm(mean_log ~ FD + FDis, data=response_both_28)
summary(mean_FD4) ## mean has positive relationship FD

par(mfrow=c(1,1))
hist(residuals(mean_FD4))
par(mfrow=c(2,2))
plot(mean_FD4) ## row 136 is an outlier - high FDis compared to other sites

stability_FD4 <- lm(stability_log ~ FD + FDis, data=response_both_28)
summary(stability_FD4) ## stability has negative relationship with FD and positive relationship with FDis

par(mfrow=c(1,1))
hist(residuals(mean_FD4))
par(mfrow=c(2,2))
plot(mean_FD4) ## row 136 is an outlier - high FDis compared to other sites
