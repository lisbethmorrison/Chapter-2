###########################################################
## Title: Analyse relationship between EF proxy and FD metrics
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: July 2019
##########################################################

rm(list=ls()) # clear R
options(scipen=999)

library(rcompanion)

## read in data
proxy_effect <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_effect.csv", header=TRUE)
proxy_response <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_response.csv", header=TRUE)
proxy_effect_response <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_effect_response.csv", header=TRUE)

## FD dendrogram results
FD_effect <- read.csv("../Data/Analysis_data/Seed dispersal/FD_effect_53spp.csv", header=TRUE)
FD_effect_response <- read.csv("../Data/Analysis_data/Seed dispersal/FD_effect_47spp.csv", header=TRUE)
FD_response <- read.csv("../Data/Analysis_data/Seed dispersal/FD_response_66spp.csv", header=TRUE)
FD_response_effect <- read.csv("../Data/Analysis_data/Seed dispersal/FD_response_47spp.csv", header=TRUE)
## multivariate FD results
FD_multi_effect <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_effect_53spp.csv", header=TRUE)
FD_multi_effect_response <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_effect_47spp.csv", header=TRUE)
FD_multi_response <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_response_66spp.csv", header=TRUE)
FD_multi_response_effect <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_response_47spp.csv", header=TRUE)

#################################
## merge files together

## effect ONLY proxy and FD
effect_53 <- merge(proxy_effect, FD_effect, by.x="GRIDREF", by.y="gridref")
## 200 sites for 53 species (effect traits only)
effect_53_2 <- merge(proxy_effect, FD_multi_effect, by.x="GRIDREF", by.y="gridref")
effect_53_2 <- na.omit(effect_53_2) ## one site doesn't have FD values
## 199 sites for 53 species (effect traits only)

## response ONLY proxy and FD
response_66 <- merge(proxy_response, FD_response, by.x="GRIDREF", by.y="gridref")
## 200 sites for 66 species (response traits only)
response_66_2 <- merge(proxy_response, FD_multi_response, by.x="GRIDREF", by.y="gridref")
## 200 sites for 66 species (response traits only)

## both effect AND response proxy and FD
## first rename FD column
names(FD_multi_effect_response)[2:5] <- c("FDVe_effect", "FDis_effect", "FRic_effect", "FDiv_effect")
names(FD_multi_response_effect)[2:5] <- c("FDVe_response", "FDis_response", "FRic_response", "FDiv_response")

colnames(FD_multi_effect_response)[2] <- "FD_effect"
colnames(FD_multi_response_effect)[2] <- "FD_response"

effect_response_47 <- merge(proxy_effect_response, FD_effect_response, by.x="GRIDREF", by.y="gridref")
effect_response_47 <- merge(effect_response_47, FD_response_effect, by.x="GRIDREF", by.y="gridref")
## 199 sites for 47 species which have response AND effect traits
effect_response_47_2 <- merge(proxy_effect_response, FD_multi_effect_response, by.x="GRIDREF", by.y="gridref")
effect_response_47_2 <- merge(effect_response_47_2, FD_multi_response_effect, by.x="GRIDREF", by.y="gridref")
## 199 sites for 47 species which have response AND effect traits

#################################

## FD dendrogram
#################################
## test for relationship between proxy and FD effect on species with ONLY effect traits (53 spp)
str(effect_53)
par(mfrow=c(1,1))
hist(effect_53$Mean_abund) ## right skew
qqnorm(effect_53$Mean_abund)
qqline(effect_53$Mean_abund,col="red")
## try sqrt transformation
effect_53$Mean_abund_sqrt = sqrt(effect_53$Mean_abund)
plotNormalHistogram(effect_53$Mean_abund_sqrt) ## much better
qqnorm(effect_53$Mean_abund_sqrt)
qqline(effect_53$Mean_abund_sqrt, col="red")
## try log transformation
effect_53$Mean_abund_log = log(effect_53$Mean_abund)
hist(effect_53$Mean_abund_log) ## even better
qqnorm(effect_53$Mean_abund_log)
qqline(effect_53$Mean_abund_log, col="red")

plotNormalHistogram(effect_53$stability) ## looks good
plotNormalHistogram(effect_53$FD_effect) ## left skew
## try square transformation
effect_53$FD_effect_sqr = (effect_53$FD_effect)^2
hist(effect_53$FD_effect_sqr) ## much better
qqnorm(effect_53$FD_effect_sqr)
qqline(effect_53$FD_effect_sqr, col="red")

model1 <- lm(Mean_abund ~ FD_effect, data=effect_53)
summary(model1)
par(mfrow=c(2,2))
plot(model1)

x = residuals(model1)
plotNormalHistogram(x)

model2 <- lm(Mean_abund_sqrt ~ FD_effect_sqr, data=effect_53)
summary(model2) 
par(mfrow=c(2,2))
plot(model2)

par(mfrow=c(1,1))
x2 = residuals(model2)
plotNormalHistogram(x2)

model3 <- lm(Mean_abund_log ~ FD_effect_sqr, data=effect_53)
summary(model3) ## positive significant (p=0.0085) higher mean = higher effect diversity
par(mfrow=c(2,2))
plot(model3) ## looks good

par(mfrow=c(1,1))
x3 = residuals(model3)
plotNormalHistogram(x3) ## residuals are normally distributed

## plot result
ggplot(effect_53, aes(x = FD_effect_sqr, y = Mean_abund_log)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black") +
  labs(x = "Effect trait functional diversity", y = "Mean proxy of function") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=5) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


model4 <- lm(stability ~ FD_effect_sqr, data=effect_53)
summary(model4) ## non-significant (p=0.671)
par(mfrow=c(2,2))
plot(model4) ## good model fit

par(mfrow=c(1,1))
x4 = residuals(model4)
plotNormalHistogram(x4) ## residuals are normally distributed

## FD multivariate metrics
str(effect_53_2)
par(mfrow=c(1,1))
hist(effect_53_2$Mean_abund) ## right skew
## try log transformation
effect_53_2$Mean_abund_log = log(effect_53_2$Mean_abund)
hist(effect_53_2$Mean_abund_log) 
qqnorm(effect_53_2$Mean_abund_log)
qqline(effect_53_2$Mean_abund_log, 
       col="red")
hist(effect_53_2$stability) ## looks good
hist(effect_53_2$FDVe) ## slight left skew
## try square transformation
effect_53_2$FDVe_sqr = (effect_53_2$FDVe)^2
hist(effect_53_2$FDVe_sqr) ## better
hist(effect_53_2$FDis) ## slight right skew
effect_53_2$FDis_log = log(effect_53_2$FDis)
hist(effect_53_2$FDis_log) ## bit better
hist(effect_53_2$FRic) ## definite right skew
effect_53_2$FRic_log = log(effect_53_2$FRic)
hist(effect_53_2$FRic_log) ## bit better (FRic probably won't be used anyway)
hist(effect_53_2$FDiv) ## looks good

## effect trait models

multi_model1 <- lm(Mean_abund_log ~ FDVe + FDis_log + FRic_log + FDiv, data=effect_53_2)
summary(multi_model1) ## FRic significant (positive)
par(mfrow=c(2,2))
plot(multi_model1)

multi_model2 <- lm(stability ~ FDVe + FDis_log + FRic_log + FDiv, data=effect_53_2)
summary(multi_model2) ## FDVe and FDis significant (both negative)
par(mfrow=c(2,2))
plot(multi_model2)

#################################
## test for relationship between proxy and FD response on species with ONLY response traits (66 spp)
str(response_66)
par(mfrow=c(1,1))
plotNormalHistogram(response_66$Mean_abund) ## right skew
qqnorm(response_66$Mean_abund,
       ylab="Sample Quantiles for Turbidity")
qqline(response_66$Mean_abund, 
       col="red")
## try sqrt transformation
response_66$Mean_abund_sqrt = sqrt(response_66$Mean_abund)
plotNormalHistogram(response_66$Mean_abund_sqrt) ## bit better
qqnorm(response_66$Mean_abund_sqrt,
       ylab="Sample Quantiles for Turbidity")
qqline(response_66$Mean_abund_sqrt, 
       col="red")
## try log transformation
response_66$Mean_abund_log = log(response_66$Mean_abund)
plotNormalHistogram(response_66$Mean_abund_log) ## much better
qqnorm(response_66$Mean_abund_log,
       ylab="Sample Quantiles for Turbidity")
qqline(response_66$Mean_abund_log, 
       col="red")

plotNormalHistogram(response_66$stability) ## looks good
plotNormalHistogram(response_66$FD_response) ## looks good

model5 <- lm(Mean_abund_log ~ FD_response, data=response_66)
summary(model5) ## just significant (p=0.045)
par(mfrow=c(2,2))
plot(model5) ## looks good

par(mfrow=c(1,1))
x5 = residuals(model5)
plotNormalHistogram(x5) ## residuals are normally distributed

## plot result
ggplot(response_66, aes(x = FD_response, y = Mean_abund_log)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "black") +
  labs(x = "Response trait functional diversity", y = "Mean proxy of function") +
  theme_bw() +
  theme(text = element_text(size = 16)) +
  labs(size=5) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


model6 <- lm(stability ~ FD_response, data=response_66)
summary(model6) ## non-significant (p=0.46)
par(mfrow=c(2,2))
plot(model6) ## good model fit

par(mfrow=c(1,1))
x6 = residuals(model6)
plotNormalHistogram(x6) ## residuals are normally distributed


## FD multivariate metrics
str(response_66_2)
par(mfrow=c(1,1))
plotNormalHistogram(response_66_2$Mean_abund) ## right skew
## try log transformation
response_66_2$Mean_abund_log = log(response_66_2$Mean_abund)
plotNormalHistogram(response_66_2$Mean_abund_log) 
qqnorm(response_66_2$Mean_abund_log,
       ylab="Sample Quantiles for Turbidity")
qqline(response_66_2$Mean_abund_log, 
       col="red")
plotNormalHistogram(response_66_2$stability) ## looks good
plotNormalHistogram(response_66_2$FDVe) ## looks good
plotNormalHistogram(response_66_2$FDis) ## right skew
response_66_2$FDis_log = log(response_66_2$FDis)
plotNormalHistogram(response_66_2$FDis_log) ## bit better
plotNormalHistogram(response_66_2$FRic) ## definite right skew
response_66_2$FRic_log = log(response_66_2$FRic)
plotNormalHistogram(response_66_2$FRic_log) ## much better (FRic probably won't be used anyway)
plotNormalHistogram(response_66_2$FDiv) ## looks good

## response trait models

multi_model3 <- lm(Mean_abund_log ~ FDVe + FDis_log + FRic_log + FDiv, data=response_66_2)
summary(multi_model3) ## FDve, FDis and FRic significant (FRic and FDVe negative, FDis positive)
par(mfrow=c(2,2))
plot(multi_model3)

multi_model4 <- lm(stability ~ FDVe + FDis_log + FRic_log + FDiv, data=response_66_2)
summary(multi_model4) ## FDve significant (negative)
par(mfrow=c(2,2))
plot(multi_model4)



#################################
## test for relationship between proxy and FD effect and FD response on species with BOTH effect and response traits (47 spp)
str(effect_response_47)
par(mfrow=c(1,1))
plotNormalHistogram(effect_response_47$Mean_abund) ## right skew
qqnorm(effect_response_47$Mean_abund,
       ylab="Sample Quantiles for Turbidity")
qqline(effect_response_47$Mean_abund, 
       col="red")
## try sqrt transformation
effect_response_47$Mean_abund_sqrt = sqrt(effect_response_47$Mean_abund)
plotNormalHistogram(effect_response_47$Mean_abund_sqrt) ## bit better
qqnorm(effect_response_47$Mean_abund_sqrt,
       ylab="Sample Quantiles for Turbidity")
qqline(effect_response_47$Mean_abund_sqrt, 
       col="red")
## try log transformation
effect_response_47$Mean_abund_log = log(effect_response_47$Mean_abund)
plotNormalHistogram(effect_response_47$Mean_abund_log) ## much better
qqnorm(effect_response_47$Mean_abund_log,
       ylab="Sample Quantiles for Turbidity")
qqline(effect_response_47$Mean_abund_log, 
       col="red")

plotNormalHistogram(effect_response_47$stability) ## looks good
plotNormalHistogram(effect_response_47$FD_response) ## looks good

model7 <- lm(Mean_abund_log ~ FD_effect, data=effect_response_47)
summary(model7) ## significant (p=0.0247)
par(mfrow=c(2,2))
plot(model7) ## looks good

par(mfrow=c(1,1))
x7 = residuals(model7)
plotNormalHistogram(x7) ## residuals are normally distributed

model8 <- lm(stability ~ FD_effect, data=effect_response_47)
summary(model8) ## non-significant (0.722)
par(mfrow=c(2,2))
plot(model8) ## good model fit

par(mfrow=c(1,1))
x8 = residuals(model8)
plotNormalHistogram(x8) ## residuals are normally distributed

model9 <- lm(Mean_abund_log ~ FD_response, data=effect_response_47)
summary(model9) ## significant (p=0.0027)
par(mfrow=c(2,2))
plot(model9) ## looks good

par(mfrow=c(1,1))
x9 = residuals(model9)
plotNormalHistogram(x9) ## residuals are normally distributed

model10 <- lm(stability ~ FD_response, data=effect_response_47)
summary(model10) ## non-significant (0.353)
par(mfrow=c(2,2))
plot(model10) ## good model fit

par(mfrow=c(1,1))
x10 = residuals(model10)
plotNormalHistogram(x10) ## residuals are normally distributed


