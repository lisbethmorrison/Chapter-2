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

## read in data - complete case detectability
proxy_seed <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_det_comp.csv", header=TRUE) ## 67 spp 200 sites
FD_multi_effect <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_effect_det_comp.csv", header=TRUE) 
FD_multi_response <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_response_det_comp.csv", header=TRUE) 
FD_multi_all <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_all_det_comp.csv", header=TRUE) 

## read in data - interpolated detectability
proxy_seed <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_det_interpol.csv", header=TRUE) ## 67 spp 200 sites
FD_multi_effect <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_effect_det_interpol.csv", header=TRUE) 
FD_multi_response <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_response_det_interpol.csv", header=TRUE) 
FD_multi_all <- read.csv("../Data/Analysis_data/Seed dispersal/FD_multi_seed_all_det_interpol.csv", header=TRUE) 

##################################################################################################################################################################

##################### DATA COLLATION FOR ANALYSIS ##################### 

colnames(FD_multi_effect)[2] <- "FDis_effect" ## rename FDis column
colnames(FD_multi_response)[2] <- "FDis_response" ## rename FDis column
colnames(FD_multi_all)[2] <- "FDis_all" ## rename FDis column
## merge ALL results into one file
seed_proxy_FD <- list(proxy_seed, FD_multi_effect, FD_multi_response, FD_multi_all) %>% reduce(left_join, by = "GRIDREF")
## save file
write.csv(seed_proxy_FD, file="../Data/Analysis_data/Seed dispersal/seed_proxy_FD_det_interpol.csv", row.names=FALSE)

####################################################################################################################################

rm(list=ls()) # clear R
## read in data
seed_proxy_FD <- read.csv("../Data/Analysis_data/Seed dispersal/seed_proxy_FD_det_interpol.csv", header=TRUE)
spp_richness <- read.csv("../Data/Analysis_data/Seed dispersal/mean_spp_richness_det_interpol.csv", header=TRUE)

## check distribution of response variables (mean and stability proxy)
hist(seed_proxy_FD$mean_abund) ## slight right skew
hist(seed_proxy_FD$stability) ## slight right skew
seed_proxy_FD$log_abund <- log(seed_proxy_FD$mean_abund)
hist(seed_proxy_FD$log_abund) ## looks better
seed_proxy_FD$log_stability <- log(seed_proxy_FD$stability)
hist(seed_proxy_FD$log_stability) ## looks better
### use the transformed variables 

## merge species richness with FD/proxy data
seed_proxy_FD <- merge(seed_proxy_FD, spp_richness, by="GRIDREF")
## correlation between spp richness and FD metrics
cor1 <- cor.test(seed_proxy_FD$FDis_effect, seed_proxy_FD$average_spp_richness, method="pearson")
cor1 ## 0.22, p=0.01
cor2 <- cor.test(seed_proxy_FD$FDis_response, seed_proxy_FD$average_spp_richness, method="pearson")
cor2 ## -0.25 p=0.003
cor3 <- cor.test(seed_proxy_FD$FDis_all, seed_proxy_FD$average_spp_richness, method="pearson")
cor3 ## -0.12, p=0.17
## low correlation between FDis and species richness
cor5 <- cor.test(seed_proxy_FD$average_spp_richness, seed_proxy_FD$log_abund, method="pearson")
cor5 ## 0.19, p=0.03
cor5 <- cor.test(seed_proxy_FD$average_spp_richness, seed_proxy_FD$log_stability, method="pearson")
cor5 ## 0.03, p=0.73

cor6 <- cor.test(seed_proxy_FD$log_abund, seed_proxy_FD$log_stability, method="pearson")
cor6 ## 0.35, p<0.001

####################################################################
######################### MEAN PROXY ###############################
####################################################################

#### FDis METRIC

## Effect traits only
mean_effect <- lm(log_abund ~ FDis_effect, data=seed_proxy_FD)
summary(mean_effect)
## significant negative - opposite of hypothesis
results_table1 <- data.frame(summary(mean_effect)$coefficients[,1:4])
results_table1$r_sq <- summary(mean_effect)$r.squared
write.csv(results_table1, file = "../Results/Seed dispersal/mean_FDis_effect_det_interpol.csv", row.names=TRUE)

## Response traits only
mean_response <- lm(log_abund ~ FDis_response, data=seed_proxy_FD)
summary(mean_response)
## significant positive relationship with response traits
results_table2 <- data.frame(summary(mean_response)$coefficients[,1:4])
results_table2$r_sq <- summary(mean_response)$r.squared
write.csv(results_table2, file = "../Results/Seed dispersal/mean_FDis_response_det_interpol.csv", row.names=TRUE)

## All traits
mean_all <- lm(log_abund ~ FDis_all, data=seed_proxy_FD)
summary(mean_all)
## significant positive relationship with all traits
results_table3 <- data.frame(summary(mean_all)$coefficients[,1:4])
results_table3$r_sq <- summary(mean_all)$r.squared
write.csv(results_table3, file = "../Results/Seed dispersal/mean_FDis_all_det_interpol.csv", row.names=TRUE)

#########################################################################
######################### STABILITY PROXY ###############################
#########################################################################

#### FDis METRIC

## Effect traits only
stability_effect <- lm(log_stability ~ FDis_effect, data=seed_proxy_FD)
summary(stability_effect)
## non significant 
results_table4 <- data.frame(summary(stability_effect)$coefficients[,1:4])
results_table4$r_sq <- summary(stability_effect)$r.squared
write.csv(results_table4, file = "../Results/Seed dispersal/stability_FDis_effect_det_interpol.csv", row.names=TRUE)

## Response traits only
stability_response <- lm(log_stability ~ FDis_response, data=seed_proxy_FD)
summary(stability_response)
## significant positive relationship with response (as hypothesised)
results_table5 <- data.frame(summary(stability_response)$coefficients[,1:4])
results_table5$r_sq <- summary(stability_response)$r.squared
write.csv(results_table5, file = "../Results/Seed dispersal/stability_FDis_response_det_interpol.csv", row.names=TRUE)

## All traits
stability_all <- lm(log_stability ~ FDis_all, data=seed_proxy_FD)
summary(stability_all)
## significant positive
results_table6 <- data.frame(summary(stability_all)$coefficients[,1:4])
results_table6$r_sq <- summary(stability_all)$r.squared
write.csv(results_table6, file = "../Results/Seed dispersal/stability_FDis_all_det_interpol.csv", row.names=TRUE)

## correlations between effect and response FDis
cor7 <- cor.test(seed_proxy_FD$FDis_effect, seed_proxy_FD$FDis_response, method="pearson")
cor7 # FDis effect vs response r=-0.23, p=0.008

## plot result mean/stability ~ FD

## mean_effect: significant negative
r2_predict <- predict(mean_effect,interval="confidence")
newdata_mean_effect <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_mean_effect$fit <- exp(newdata_mean_effect$fit)
newdata_mean_effect$upr <- exp(newdata_mean_effect$upr)
newdata_mean_effect$lwr <- exp(newdata_mean_effect$lwr)

mean_effect <- ggplot(newdata_mean_effect, aes(x = FDis_effect, y = mean_abund)) + 
  geom_point(colour="black", size=0.4) +
  geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Mean abundance", x=expression("F"[DIS]*" effect traits")) + 
  scale_y_continuous(breaks=seq(200,1600,200)) +
  theme_classic() +
  theme(text = element_text(size=7))
mean_effect

## mean_response: significant positive
r2_predict <- predict(mean_response,interval="confidence")
newdata_mean_response <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_mean_response$fit <- exp(newdata_mean_response$fit)
newdata_mean_response$upr <- exp(newdata_mean_response$upr)
newdata_mean_response$lwr <- exp(newdata_mean_response$lwr)

mean_response <- ggplot(newdata_mean_response, aes(x = FDis_response, y = mean_abund)) + 
  geom_point(colour="black", size=0.4) +
  geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Mean abundance", x=expression("F"[DIS]*" response traits")) + 
  scale_y_continuous(breaks=seq(200,1600,200)) +
  theme_classic() +
  theme(text = element_text(size=7))
mean_response

## mean_all: significant positive
r2_predict <- predict(mean_all,interval="confidence")
newdata_mean_all <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_mean_all$fit <- exp(newdata_mean_all$fit)
newdata_mean_all$upr <- exp(newdata_mean_all$upr)
newdata_mean_all$lwr <- exp(newdata_mean_all$lwr)

mean_all <- ggplot(newdata_mean_all, aes(x = FDis_all, y = mean_abund)) + 
  geom_point(colour="black", size=0.4) +
  geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Mean abundance", x=expression("F"[DIS]*" all traits")) + 
  scale_y_continuous(breaks=seq(200,1600,200)) +
  theme_classic() +
  theme(text = element_text(size=7))
mean_all

### stability_effect: non-significant

stability_effect <- ggplot(seed_proxy_FD, aes(x = FDis_effect, y = stability)) + 
  geom_point(colour="black", size=0.4) +
  #geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  #geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Stability of abundance", x=expression("F"[DIS]*" effect traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
stability_effect

## stability_response: significant positive
r2_predict <- predict(stability_response,interval="confidence")
newdata_stab_response <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_stab_response$fit <- exp(newdata_stab_response$fit)
newdata_stab_response$upr <- exp(newdata_stab_response$upr)
newdata_stab_response$lwr <- exp(newdata_stab_response$lwr)

stability_response <- ggplot(newdata_stab_response, aes(x = FDis_response, y = stability)) + 
  geom_point(colour="black", size=0.4) +
  geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Stability of abundance", x=expression("F"[DIS]*" response traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
stability_response

## stability_all: significant positive
r2_predict <- predict(stability_all,interval="confidence")
newdata_stab_all <- cbind(data.frame(seed_proxy_FD), data.frame(r2_predict))

newdata_stab_all$fit <- exp(newdata_stab_all$fit)
newdata_stab_all$upr <- exp(newdata_stab_all$upr)
newdata_stab_all$lwr <- exp(newdata_stab_all$lwr)

stability_all <- ggplot(newdata_stab_all, aes(x = FDis_all, y = stability)) + 
  geom_point(colour="black", size=0.4) +
  geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Stability of abundance", x=expression("F"[DIS]*" all traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
stability_all


## save all 6 plots as figure 1
fig1 <- plot_grid(mean_effect, mean_all, mean_response, stability_effect, stability_all, stability_response, 
                  labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), label_size=8, ncol = 3, 
                  nrow = 2, vjust =c(1,1,1,1,1,1), hjust=0) 
fig1
ggsave(file="../Graphs/Figure1_seed_det_interpol.png", fig1, width = 120, height = 80, dpi = 600, units = "mm", device='png') 

bes1 <- plot_grid(mean_effect, mean_all, mean_response, stability_effect, stability_all, stability_response, ncol = 3, nrow = 2)
ggsave(file="../Graphs/BES_seed.png", bes1, height=7, width=10) 
