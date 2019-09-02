###########################################################
## Title: Analyse relationship between EF proxy and FD metrics
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: August 2019
##########################################################

rm(list=ls()) # clear R
options(scipen=999)

library(rcompanion)

## read in data
proxy_effect <- read.csv("../Data/Analysis_data/Pest control/BBS_proxy_invert_effect.csv", header=TRUE) ## 48 spp 199 sites
proxy_response <- read.csv("../Data/Analysis_data/Pest control/BBS_proxy_invert_response.csv", header=TRUE) ## 63 spp 200 sites
proxy_all <- read.csv("../Data/Analysis_data/Pest control/BBS_proxy_invert_all.csv", header=TRUE) ## 42 spp 199 sites

## FD dendrogram results
FD_effect <- read.csv("../Data/Analysis_data/Pest control/FD_effect_48spp.csv", header=TRUE) ## full effect traits
FD_effect2 <- read.csv("../Data/Analysis_data/Pest control/FD_effect_42spp.csv", header=TRUE) ## subset effect traits
FD_response <- read.csv("../Data/Analysis_data/Pest control/FD_response_63spp.csv", header=TRUE) ## full response traits
FD_response2 <- read.csv("../Data/Analysis_data/Pest control/FD_response_42spp.csv", header=TRUE) ## subset response traits
FD_effect_both <- read.csv("../Data/Analysis_data/Pest control/FD_effect_both_48spp.csv", header=TRUE) ## full effect and both traits
FD_effect_both2 <- read.csv("../Data/Analysis_data/Pest control/FD_effect_both_42spp.csv", header=TRUE) ## subset effect and both traits
FD_response_both <- read.csv("../Data/Analysis_data/Pest control/FD_response_both_63spp.csv", header=TRUE) ## full response and both traits
FD_response_both2 <- read.csv("../Data/Analysis_data/Pest control/FD_response_both_42spp.csv", header=TRUE) ## subset response and both traits

## multivariate FD results
FD_multi_effect <- read.csv("../Data/Analysis_data/Pest control/FD_multi_effect_48spp.csv", header=TRUE) ## full effect traits
FD_multi_effect2 <- read.csv("../Data/Analysis_data/Pest control/FD_multi_effect_42spp.csv", header=TRUE) ## subset effect traits
FD_multi_response <- read.csv("../Data/Analysis_data/Pest control/FD_multi_response_63spp.csv", header=TRUE) ## full response traits
FD_multi_response2 <- read.csv("../Data/Analysis_data/Pest control/FD_multi_response_42spp.csv", header=TRUE) ## subset response traits
FD_multi_effect_both <- read.csv("../Data/Analysis_data/Pest control/FD_multi_effect_both_48spp.csv", header=TRUE) ## full effect and both traits
FD_multi_effect_both2 <- read.csv("../Data/Analysis_data/Pest control/FD_multi_effect_both_42spp.csv", header=TRUE) ## subset effect and both traits
FD_multi_response_both <- read.csv("../Data/Analysis_data/Pest control/FD_multi_response_both_63spp.csv", header=TRUE) ## full response and both traits
FD_multi_response_both2 <- read.csv("../Data/Analysis_data/Pest control/FD_multi_response_both_42spp.csv", header=TRUE) ## subset response and both traits

##################################################################################################################################################################

##################### DATA COLLATION FOR ANALYSIS ##################### 

## effect ONLY species and sites (n=48)
## merge FD dendro and multivariate together
effect_48 <- merge(FD_effect, FD_multi_effect, by="gridref")
## merge with proxy for 32 species
effect_48 <- merge(effect_48, proxy_effect, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(effect_48)[2] <- "FD"
## save file
write.csv(effect_48, file="../Data/Analysis_data/Pest control/Effect_48spp_final.csv", row.names=FALSE)
## 48 species at 199 sites

## response ONLY species and sites (n=63)
## merge FD dendro and multivariate together
response_63 <- merge(FD_response, FD_multi_response, by="gridref")
## merge with proxy for 32 species
response_63 <- merge(response_63, proxy_response, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(response_63)[2] <- "FD"
## save file
write.csv(response_63, file="../Data/Analysis_data/Pest control/Response_63spp_final.csv", row.names=FALSE)
## 63 species at 200 sites

## effect and both only species and sites (n=48)
## merge FD dendro and multivariate together
effect_both_48 <- merge(FD_effect_both, FD_multi_effect_both, by="gridref")
## merge with proxy for 48 species
effect_both_48 <- merge(effect_both_48, proxy_effect, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(effect_both_48)[2] <- "FD"
## save file
write.csv(effect_both_48, file="../Data/Analysis_data/Pest control/Effect_both_48spp_final.csv", row.names=FALSE)
## 48 species at 199 sites

## response and both only species and sites (n=63)
## merge FD dendro and multivariate together
response_both_63 <- merge(FD_response_both, FD_multi_response, by="gridref")
## merge with proxy for 63 species
response_both_63 <- merge(response_both_63, proxy_response, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(response_both_63)[2] <- "FD"
## save file
write.csv(response_both_63, file="../Data/Analysis_data/Pest control/Response_both_63spp_final.csv", row.names=FALSE)
## 63 species at 200 sites

############## same as above for subset of species

## effect ONLY species and sites (n=42)
## merge FD dendro and multivariate together
effect_42_sub <- merge(FD_effect2, FD_multi_effect2, by="gridref")
## merge with proxy for 42 species
effect_42_sub <- merge(effect_42_sub, proxy_effect, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(effect_42_sub)[2] <- "FD"
## save file
write.csv(effect_42_sub, file="../Data/Analysis_data/Pest control/Effect_42spp_final.csv", row.names=FALSE)
## 42 species at 199 sites

## response ONLY species and sites (n=42)
## merge FD dendro and multivariate together
response_42_sub <- merge(FD_response2, FD_multi_response2, by="gridref")
## merge with proxy for 42 species
response_42_sub <- merge(response_42_sub, proxy_response, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(response_42_sub)[2] <- "FD"
## save file
write.csv(response_42_sub, file="../Data/Analysis_data/Pest control/Response_42spp_final.csv", row.names=FALSE)
## 42 species at 200 sites

## effect and both only species and sites (n=42)
## merge FD dendro and multivariate together
effect_both_42_sub <- merge(FD_effect_both2, FD_multi_effect_both2, by="gridref")
## merge with proxy for 42 species
effect_both_42_sub <- merge(effect_both_42_sub, proxy_effect, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(effect_both_42_sub)[2] <- "FD"
## save file
write.csv(effect_both_42_sub, file="../Data/Analysis_data/Pest control/Effect_both_42spp_final.csv", row.names=FALSE)
## 42 species at 199 sites

## response and both only species and sites (n=42)
## merge FD dendro and multivariate together
response_both_42_sub <- merge(FD_response_both2, FD_multi_response2, by="gridref")
## merge with proxy for 42 species
response_both_42_sub <- merge(response_both_42_sub, proxy_response, by.x="gridref", by.y="GRIDREF")
## re-name FD_effect to just FD (dendrogram metric)
names(response_both_42_sub)[2] <- "FD"
## save file
write.csv(response_both_42_sub, file="../Data/Analysis_data/Pest control/Response_both_42spp_final.csv", row.names=FALSE)
## 42 species at 200 sites

####################################################################################################################################

rm(list=ls()) # clear R

## analysis on full dataset (NOT subset of 28 species)
## read in data
effect_48 <- read.csv("../Data/Analysis_data/Pest control/Effect_48spp_final.csv", header=TRUE)
response_63 <- read.csv("../Data/Analysis_data/Pest control/Response_63spp_final.csv", header=TRUE)
effect_both_48 <- read.csv("../Data/Analysis_data/Pest control/Effect_both_48spp_final.csv", header=TRUE)
response_both_63 <- read.csv("../Data/Analysis_data/Pest control/Response_both_63spp_final.csv", header=TRUE)

## keep only gridref, FD, FDis, mean and stability from each dataset
effect_48 <- effect_48[,c(1,2,4,7:8)]
response_63 <- response_63[,c(1,2,4,7:8)]
effect_both_48 <- effect_both_48[,c(1,2,4,7:8)]
response_both_63 <- response_both_63[,c(1,2,4,7:8)]

#######################################################
#### Hypothesis 1: effect trait diversity correlates with mean function, not stability
str(effect_48)
hist(effect_48$mean) ## right skew
effect_48$mean_log <- log(effect_48$mean)
hist(effect_48$mean_log) ## much better
hist(effect_48$stability) ## looks good

## run models
mean_FD <- lm(mean_log ~ FD + FDis, data=effect_48)
summary(mean_FD) ## non-significant

par(mfrow=c(1,1))
hist(residuals(mean_FD))
par(mfrow=c(2,2))
plot(mean_FD) 

stability_FD <- lm(stability ~ FD + FDis, data=effect_48)
summary(stability_FD) ## stability has negative relationship with FD 

par(mfrow=c(1,1))
hist(residuals(stability_FD))
par(mfrow=c(2,2))
plot(stability_FD) 
par(mfrow=c(1,1))

#######################################################
#### Hypothesis 2: response trait diversity correlates with stability of function, not mean
str(response_63)
hist(response_63$mean) ## right skew
response_63$mean_log <- log(response_63$mean)
hist(response_63$mean_log) ## much better
hist(response_63$stability) ## good

## run models
mean_FD2 <- lm(mean_log ~ FD + FDis, data=response_63)
summary(mean_FD2) ## non-significant

par(mfrow=c(1,1))
hist(residuals(mean_FD2))
par(mfrow=c(2,2))
plot(mean_FD2) 

stability_FD2 <- lm(stability ~ FD + FDis, data=response_63)
summary(stability_FD2) ## stability has negative relationship with FD

par(mfrow=c(1,1))
hist(residuals(stability_FD2))
par(mfrow=c(2,2))
plot(stability_FD2) 
par(mfrow=c(1,1))

#######################################################
#### Hypothesis 3: effect and both trait diversity correlates more strongly with mean compared to stability
str(effect_both_48)
hist(effect_both_48$mean) ## right skew
effect_both_48$mean_log <- log(effect_both_48$mean)
hist(effect_both_48$mean_log) ## much better
hist(effect_both_48$stability) ## slight right skew
effect_both_48$stability_log <- log(effect_both_48$stability)
hist(effect_both_48$stability_log) ## better

## run models
mean_FD3 <- lm(mean_log ~ FD + FDis, data=effect_both_48)
summary(mean_FD3) ## non-significant (FDis negative p=0.06)

par(mfrow=c(1,1))
hist(residuals(mean_FD3))
par(mfrow=c(2,2))
plot(mean_FD3) 

stability_FD3 <- lm(stability_log ~ FD + FDis, data=effect_both_48)
summary(stability_FD3) ## stability has negative relationship with FD 

par(mfrow=c(1,1))
hist(residuals(stability_FD3))
par(mfrow=c(2,2))
plot(stability_FD3) 
par(mfrow=c(1,1))

#######################################################
#### Hypothesis 4: response and both trait diversity correlates more strongly with stability compared to mean
str(response_both_63)
hist(response_both_63$mean) ## right skew
response_both_63$mean_log <- log(response_both_63$mean)
hist(response_both_63$mean_log) ## much better
hist(response_both_63$stability) ## good

## run models
mean_FD4 <- lm(mean_log ~ FD + FDis, data=response_both_63)
summary(mean_FD4) ## non-significant

par(mfrow=c(1,1))
hist(residuals(mean_FD4))
par(mfrow=c(2,2))
plot(mean_FD4) ## row 136 is an outlier - high FDis compared to other sites

stability_FD4 <- lm(stability ~ FD + FDis, data=response_both_63)
summary(stability_FD4) ## stability has negative relationship with FD

par(mfrow=c(1,1))
hist(residuals(mean_FD4))
par(mfrow=c(2,2))
plot(mean_FD4) ## row 136 is an outlier - high FDis compared to other sites


##########################################################
######## REPEAT ANALYSIS ON SUBSET OF SPECIES ############
##########################################################

rm(list=ls()) # clear R

## read in data
effect_42 <- read.csv("../Data/Analysis_data/Pest control/Effect_42spp_final.csv", header=TRUE)
response_42 <- read.csv("../Data/Analysis_data/Pest control/Response_42spp_final.csv", header=TRUE)
effect_both_42 <- read.csv("../Data/Analysis_data/Pest control/Effect_both_42spp_final.csv", header=TRUE)
response_both_42 <- read.csv("../Data/Analysis_data/Pest control/Response_both_42spp_final.csv", header=TRUE)

## keep only gridref, FD, FDis, mean and stability from each dataset
effect_42 <- effect_42[,c(1,2,4,7:8)]
response_42 <- response_42[,c(1,2,4,7:8)]
effect_both_42 <- effect_both_42[,c(1,2,4,7:8)]
response_both_42 <- response_both_42[,c(1,2,4,7:8)]

#######################################################
#### Hypothesis 1: effect trait diversity correlates with mean function, not stability
str(effect_42)
hist(effect_42$mean) ## right skew
effect_42$mean_log <- log(effect_42$mean)
hist(effect_42$mean_log) ## much better
hist(effect_42$stability) ## looks good

## run models
mean_FD <- lm(mean_log ~ FD + FDis, data=effect_42)
summary(mean_FD) ## non-significant

par(mfrow=c(1,1))
hist(residuals(mean_FD))
par(mfrow=c(2,2))
plot(mean_FD) 

stability_FD <- lm(stability ~ FD + FDis, data=effect_42)
summary(stability_FD) ## stability has negative relationship with FD 

par(mfrow=c(1,1))
hist(residuals(stability_FD))
par(mfrow=c(2,2))
plot(stability_FD) 
par(mfrow=c(1,1))

#######################################################
#### Hypothesis 2: response trait diversity correlates with stability of function, not mean
str(response_42)
hist(response_42$mean) ## right skew
response_42$mean_log <- log(response_42$mean)
hist(response_42$mean_log) ## much better
hist(response_42$stability) ## good

## run models
mean_FD2 <- lm(mean_log ~ FD + FDis, data=response_42)
summary(mean_FD2) ## non-significant

par(mfrow=c(1,1))
hist(residuals(mean_FD2))
par(mfrow=c(2,2))
plot(mean_FD2) 

stability_FD2 <- lm(stability ~ FD + FDis, data=response_42)
summary(stability_FD2) ## Non-significant

par(mfrow=c(1,1))
hist(residuals(stability_FD2))
par(mfrow=c(2,2))
plot(stability_FD2) 
par(mfrow=c(1,1))

#######################################################
#### Hypothesis 3: effect and both trait diversity correlates more strongly with mean compared to stability
str(effect_both_42)
hist(effect_both_42$mean) ## right skew
effect_both_42$mean_log <- log(effect_both_42$mean)
hist(effect_both_42$mean_log) ## much better
hist(effect_both_42$stability) ## slight right skew
effect_both_42$stability_log <- log(effect_both_42$stability)
hist(effect_both_42$stability_log) ## better

## run models
mean_FD3 <- lm(mean_log ~ FD + FDis, data=effect_both_42)
summary(mean_FD3) ## non-significant (FDis negative p=0.06)

par(mfrow=c(1,1))
hist(residuals(mean_FD3))
par(mfrow=c(2,2))
plot(mean_FD3) 

stability_FD3 <- lm(stability_log ~ FD + FDis, data=effect_both_42)
summary(stability_FD3) ## Non-significant 

par(mfrow=c(1,1))
hist(residuals(stability_FD3))
par(mfrow=c(2,2))
plot(stability_FD3) 
par(mfrow=c(1,1))

#######################################################
#### Hypothesis 4: response and both trait diversity correlates more strongly with stability compared to mean
str(response_both_42)
hist(response_both_42$mean) ## right skew
response_both_42$mean_log <- log(response_both_42$mean)
hist(response_both_42$mean_log) ## much better
hist(response_both_42$stability) ## good

## run models
mean_FD4 <- lm(mean_log ~ FD + FDis, data=response_both_42)
summary(mean_FD4) ## non-significant

par(mfrow=c(1,1))
hist(residuals(mean_FD4))
par(mfrow=c(2,2))
plot(mean_FD4) 

stability_FD4 <- lm(stability ~ FD + FDis, data=response_both_42)
summary(stability_FD4) ## non-significant

par(mfrow=c(1,1))
hist(residuals(mean_FD4))
par(mfrow=c(2,2))
plot(mean_FD4) 
