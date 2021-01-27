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
library(tidyverse)
library(gridExtra)
library(cowplot)

## read in data
proxy_seed <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_det2.csv", header=TRUE) ## 67 spp 200 sites
## multivariate FD results
FD_multi_effect <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_effect_det2.csv", header=TRUE) 
FD_multi_response <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_response_det2.csv", header=TRUE) 
FD_multi_all <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_all_det2.csv", header=TRUE) 

##################################################################################################################################################################

##################### DATA COLLATION FOR ANALYSIS ##################### 
## remove other multivariate FD metrics that we're not using
FD_multi_effect <- FD_multi_effect[,c(1,3)]
colnames(FD_multi_effect)[2] <- "FDis_effect" ## rename FDis column
FD_multi_response <- FD_multi_response[,c(1,3)]
colnames(FD_multi_response)[2] <- "FDis_response" ## rename FDis column
FD_multi_all <- FD_multi_all[,c(1,3)]
colnames(FD_multi_all)[2] <- "FDis_all" ## rename FDis column
## merge ALL results into one file
seed_proxy_FD <- list(proxy_seed, FD_multi_effect, FD_multi_response, FD_multi_all) %>% reduce(left_join, by = "GRIDREF")
seed_proxy_FD <- na.omit(seed_proxy_FD) ## remove site which only has one species
## save file
write.csv(seed_proxy_FD, file="../Data/Analysis_data/Seed dispersal/seed_proxy_FD_det2.csv", row.names=FALSE)

####################################################################################################################################

rm(list=ls()) # clear R
## read in data
seed_proxy_FD <- read.csv("../Data/Analysis_data/Seed dispersal/seed_proxy_FD_det2.csv", header=TRUE)
spp_richness <- read.csv("../Data/Analysis_data/Seed dispersal/mean_spp_richness.csv", header=TRUE)

## check distribution of response variables (mean and stability proxy)
hist(seed_proxy_FD$mean_abund) ## slight right skew
hist(seed_proxy_FD$stability) ## slight right skew
seed_proxy_FD$sqrt_abund <- sqrt(seed_proxy_FD$mean_abund)
hist(seed_proxy_FD$sqrt_abund) ## looks better
seed_proxy_FD$sqrt_stability <- sqrt(seed_proxy_FD$stability)
hist(seed_proxy_FD$sqrt_stability) ## looks better
### use the transformed variables 

## merge species richness with FD/proxy data
seed_proxy_FD <- merge(seed_proxy_FD, spp_richness, by="GRIDREF")
## correlation between spp richness and FD metrics
cor1 <- cor.test(seed_proxy_FD$FDis_effect, seed_proxy_FD$average_spp_richness, method="pearson")
cor1 ## 0.13, p=0.06
cor2 <- cor.test(seed_proxy_FD$FDis_response, seed_proxy_FD$average_spp_richness, method="pearson")
cor2 ## 0.-0.008, p=0.9
cor3 <- cor.test(seed_proxy_FD$FDis_all, seed_proxy_FD$average_spp_richness, method="pearson")
cor3 ## 0.04, p=0.54
## low correlation between FDis and species richness
cor5 <- cor.test(seed_proxy_FD$average_spp_richness, seed_proxy_FD$sqrt_abund, method="pearson")
cor5 ## 0.52, p<0.001
cor5 <- cor.test(seed_proxy_FD$average_spp_richness, seed_proxy_FD$sqrt_stability, method="pearson")
cor5 ## 0.29, p<0.001

cor6 <- cor.test(seed_proxy_FD$sqrt_abund, seed_proxy_FD$sqrt_stability, method="pearson")
cor6 ## 0.23, p=0.0012

####################################################################
######################### MEAN PROXY ###############################
####################################################################

#### FD DENDROGRAM METRIC

## Effect traits only
mean_FD1 <- lm(sqrt_abund ~ FD_effect, data=seed_proxy_FD)
summary(mean_FD1)
par(mfrow=c(2,2))
plot(mean_FD1) ## looks good
## significant positive relationship with effect (as hypothesised)

## Effect and both traits
mean_FD2 <- lm(sqrt_abund ~ FD_effect_both, data=seed_proxy_FD)
summary(mean_FD2)
## significant positive relationship with effect and both (as hypothesised)

## Response traits only
mean_FD3 <- lm(sqrt_abund ~ FD_response, data=seed_proxy_FD)
summary(mean_FD3)
## significant positive relationship with response 

## Response and both traits
mean_FD4 <- lm(sqrt_abund ~ FD_response_both, data=seed_proxy_FD)
summary(mean_FD4)
## significant positive relationship with response and both 

## All traits
mean_FD5 <- lm(sqrt_abund ~ FD_all, data=seed_proxy_FD)
summary(mean_FD5)
## significant positive relationship with all

#### FDis METRIC

## Effect traits only
mean_FD6 <- lm(sqrt_abund ~ FDis_effect, data=seed_proxy_FD)
summary(mean_FD6)
## non-significant
results_table1 <- data.frame(summary(mean_FD6)$coefficients[,1:4])
write.csv(results_table1, file = "../Results/Seed dispersal/mean_FDis_effect.csv", row.names=TRUE)

## Effect and both traits
mean_FD7 <- lm(sqrt_abund ~ FDis_effect_both, data=seed_proxy_FD)
summary(mean_FD7)
## non-significant

## Response traits only
mean_FD8 <- lm(sqrt_abund ~ FDis_response, data=seed_proxy_FD)
summary(mean_FD8)
## significant positive relationship with response traits
results_table2 <- data.frame(summary(mean_FD8)$coefficients[,1:4])
write.csv(results_table2, file = "../Results/Seed dispersal/mean_FDis_response.csv", row.names=TRUE)

## Response and both traits
mean_FD9 <- lm(sqrt_abund ~ FDis_response_both, data=seed_proxy_FD)
summary(mean_FD9)
## significant positive relationship with response and both traits

## All traits
mean_FD10 <- lm(sqrt_abund ~ FDis_all, data=seed_proxy_FD)
summary(mean_FD10)
## significant positive relationship with all traits
results_table3 <- data.frame(summary(mean_FD10)$coefficients[,1:4])
write.csv(results_table3, file = "../Results/Seed dispersal/mean_FDis_all.csv", row.names=TRUE)

#########################################################################
######################### STABILITY PROXY ###############################
#########################################################################
#### FD DENDROGRAM METRIC

## Effect traits only
stability_FD1 <- lm(sqrt_stability ~ FD_effect, data=seed_proxy_FD)
summary(stability_FD1)
## significant positive

## Effect and both traits
stability_FD2 <- lm(sqrt_stability ~ FD_effect_both, data=seed_proxy_FD)
summary(stability_FD2)
## non-significant

## Response traits only
stability_FD3 <- lm(sqrt_stability ~ FD_response, data=seed_proxy_FD)
summary(stability_FD3)
## significant positive relationship with response (as hypothesised)

## Response and both traits
stability_FD4 <- lm(sqrt_stability ~ FD_response_both, data=seed_proxy_FD)
summary(stability_FD4)
## significant positive relationship with response and both (as hypothesised)

## all traits
stability_FD5 <- lm(sqrt_stability ~ FD_all, data=seed_proxy_FD)
summary(stability_FD5)
## significant positive relationship with all

#### FDis METRIC

## Effect traits only
stability_FD6 <- lm(sqrt_stability ~ FDis_effect, data=seed_proxy_FD)
summary(stability_FD6)
## significant negative relationship
results_table4 <- data.frame(summary(stability_FD6)$coefficients[,1:4])
write.csv(results_table4, file = "../Results/Seed dispersal/stability_FDis_effect.csv", row.names=TRUE)

## Effect and both traits
stability_FD7 <- lm(sqrt_stability ~ FDis_effect_both, data=seed_proxy_FD)
summary(stability_FD7)
## significant negative relationship with effect and both

## Response traits only
stability_FD8 <- lm(sqrt_stability ~ FDis_response, data=seed_proxy_FD)
summary(stability_FD8)
## significant positive relationship with response (as hypothesised)
results_table5 <- data.frame(summary(stability_FD8)$coefficients[,1:4])
write.csv(results_table5, file = "../Results/Seed dispersal/stability_FDis_response.csv", row.names=TRUE)

## Response and both traits
stability_FD9 <- lm(sqrt_stability ~ FDis_response_both, data=seed_proxy_FD)
summary(stability_FD9)
## significant positive relationship with response and both

## All traits
stability_FD10 <- lm(sqrt_stability ~ FDis_all, data=seed_proxy_FD)
summary(stability_FD10)
## non-significant
results_table6 <- data.frame(summary(stability_FD10)$coefficients[,1:4])
write.csv(results_table6, file = "../Results/Seed dispersal/stability_FDis_all.csv", row.names=TRUE)

## correlation between FD and FDis for each combination of traits
cor1 <- cor(seed_proxy_FD[,c(4,6)])
cor1 ## 0.6
cor2 <- cor(seed_proxy_FD[,c(5,7)])
cor2 ## 0.57
cor3 <- cor(seed_proxy_FD[,c(8,10)])
cor3 ## 0.48
cor4 <- cor(seed_proxy_FD[,c(9,11)])
cor4 ## 0.46
cor5 <- cor(seed_proxy_FD[,c(12,13)])
cor5 ## 0.47
## weak correlations between FD and FDis

## correlations between trait combinationsn for each metric
cor6 <- cor.test(seed_proxy_FD$FD_effect, seed_proxy_FD$FD_response, method="pearson")
cor6 # FD effect vs response 0.62
cor7 <- cor.test(seed_proxy_FD$FD_effect_both, seed_proxy_FD$FD_response_both, method="pearson")
cor7 # FD effect_both vs response_both 0.67
cor8 <- cor.test(seed_proxy_FD$FDis_effect, seed_proxy_FD$FDis_response, method="pearson")
cor8 # FDis effect vs response 0.068
cor9 <- cor.test(seed_proxy_FD$FDis_effect_both, seed_proxy_FD$FDis_response_both, method="pearson")
cor9 # FDis effect_both vs response_both 0.33
## species with high FD effect diversity are likely to have high FD response diveristy
## but not with FDis

## plot result stability ~ FD

r2_predict <- predict(mean_FD6,interval="confidence")
newdata_mean_effect <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_mean_effect$fit <- (newdata_mean_effect$fit)^2
newdata_mean_effect$upr <- (newdata_mean_effect$upr)^2
newdata_mean_effect$lwr <- (newdata_mean_effect$lwr)^2

mean_effect <- ggplot(newdata_mean_effect, aes(x = FDis_effect, y = mean_abund)) + 
  geom_point(colour="black", size=0.4) +
  geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Mean abundance", x=expression("F"[DIS]*" effect traits")) + 
  scale_y_continuous(breaks=seq(200,1600,200)) +
  theme_classic() +
  theme(text = element_text(size=7))
mean_effect

r2_predict <- predict(mean_FD8,interval="confidence")
newdata_mean_response <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_mean_response$fit <- (newdata_mean_response$fit)^2
newdata_mean_response$upr <- (newdata_mean_response$upr)^2
newdata_mean_response$lwr <- (newdata_mean_response$lwr)^2

mean_response <- ggplot(newdata_mean_response, aes(x = FDis_response, y = mean_abund)) + 
  geom_point(colour="black", size=0.4) +
  geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Mean abundance", x=expression("F"[DIS]*" response traits")) + 
  scale_y_continuous(breaks=seq(200,1600,200)) +
  theme_classic() +
  theme(text = element_text(size=7))
mean_response

r2_predict <- predict(mean_FD10,interval="confidence")
newdata_mean_all <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_mean_all$fit <- (newdata_mean_all$fit)^2
newdata_mean_all$upr <- (newdata_mean_all$upr)^2
newdata_mean_all$lwr <- (newdata_mean_all$lwr)^2

mean_all <- ggplot(newdata_mean_all, aes(x = FDis_all, y = mean_abund)) + 
  geom_point(colour="black", size=0.4) +
  geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Mean abundance", x=expression("F"[DIS]*" all traits")) + 
  scale_y_continuous(breaks=seq(200,1600,200)) +
  theme_classic() +
  theme(text = element_text(size=7))
mean_all

grid.arrange(mean_effect, mean_response, mean_all, nrow=1)

r2_predict <- predict(stability_FD6,interval="confidence")
newdata_stab_effect <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_stab_effect$fit <- (newdata_stab_effect$fit)^2
newdata_stab_effect$upr <- (newdata_stab_effect$upr)^2
newdata_stab_effect$lwr <- (newdata_stab_effect$lwr)^2

stability_effect <- ggplot(seed_proxy_FD, aes(x = FDis_effect, y = stability)) + 
  geom_point(colour="black", size=0.4) +
  #geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  #geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Stability of abundance", x=expression("F"[DIS]*" effect traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
stability_effect

r2_predict <- predict(stability_FD8,interval="confidence")
newdata_stab_response <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_stab_response$fit <- (newdata_stab_response$fit)^2
newdata_stab_response$upr <- (newdata_stab_response$upr)^2
newdata_stab_response$lwr <- (newdata_stab_response$lwr)^2

stability_response <- ggplot(newdata_stab_response, aes(x = FDis_response, y = stability)) + 
  geom_point(colour="black", size=0.4) +
  geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Stability of abundance", x=expression("F"[DIS]*" response traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
stability_response

r2_predict <- predict(stability_FD10,interval="confidence")
newdata_stab_all <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_stab_all$fit <- (newdata_stab_all$fit)^2
newdata_stab_all$upr <- (newdata_stab_all$upr)^2
newdata_stab_all$lwr <- (newdata_stab_all$lwr)^2

stability_all <- ggplot(seed_proxy_FD, aes(x = FDis_all, y = stability)) + 
  geom_point(colour="black", size=0.4) +
  #geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  #geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Stability of abundance", x=expression("F"[DIS]*" all traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
stability_all

grid.arrange(stability_effect, stability_response, stability_all, nrow=1)

## save all 6 plots as figure 1
fig1 <- plot_grid(mean_effect, mean_all, mean_response, stability_effect, stability_all, stability_response, 
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), label_size=8, ncol = 3, 
                  nrow = 2, vjust =c(1,1,1,1,1,1), hjust=0) 
fig1
ggsave(file="../Graphs/Figure1_seed_det2.png", fig1, width = 120, height = 80, dpi = 600, units = "mm", device='png') 

bes1 <- plot_grid(mean_effect, mean_all, mean_response, stability_effect, stability_all, stability_response, ncol = 3, nrow = 2)
ggsave(file="../Graphs/BES_seed.png", bes1, height=7, width=10) 
