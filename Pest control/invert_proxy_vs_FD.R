###########################################################
## Title: Analyse relationship between EF proxy and FD metrics
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: August 2019
##########################################################

rm(list=ls()) # clear R
options(scipen=999)

library(rcompanion)
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(cowplot)

## read in data - complete case detectability
proxy_invert <- read.csv("../Data/Analysis_data/Pest control/BBS_proxy_invert_det_comp.csv", header=TRUE) ## 67 spp 200 sites
FD_multi_effect <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_effect_det_comp.csv", header=TRUE) 
FD_multi_response <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_response_det_comp.csv", header=TRUE) 
FD_multi_all <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_all_det_comp.csv", header=TRUE) 

## read in data - interpolated detectability
proxy_invert <- read.csv("../Data/Analysis_data/Pest control/BBS_proxy_invert_det_interpol.csv", header=TRUE) ## 67 spp 200 sites
FD_multi_effect <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_effect_det_interpol.csv", header=TRUE) 
FD_multi_response <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_response_det_interpol.csv", header=TRUE) 
FD_multi_all <- read.csv("../Data/Analysis_data/Pest control/FD_multi_invert_all_det_interpol.csv", header=TRUE) 

##################################################################################################################################################################

##################### DATA COLLATION FOR ANALYSIS ##################### 
colnames(FD_multi_effect)[2] <- "FDis_effect" ## rename FDis column
colnames(FD_multi_response)[2] <- "FDis_response" ## rename FDis column
colnames(FD_multi_all)[2] <- "FDis_all" ## rename FDis column
## merge ALL results into one file
invert_proxy_FD <- list(proxy_invert, FD_multi_effect, FD_multi_response, FD_multi_all) %>% reduce(left_join, by = "GRIDREF")
## save file
write.csv(invert_proxy_FD, file="../Data/Analysis_data/Pest control/invert_proxy_FD_det_interpol.csv", row.names=FALSE)

####################################################################################################################################

rm(list=ls()) # clear R
## read in data
invert_proxy_FD <- read.csv("../Data/Analysis_data/Pest control/invert_proxy_FD_det_interpol.csv", header=TRUE)
spp_richness <- read.csv("../Data/Analysis_data/Pest Control/mean_spp_richness_det_interpol.csv", header=TRUE)

## check distribution of response variables (mean and stability proxy)
hist(invert_proxy_FD$mean_abund) ## slight right skew
hist(invert_proxy_FD$stability) ## slight right skew
invert_proxy_FD$log_abund <- log(invert_proxy_FD$mean_abund)
hist(invert_proxy_FD$log_abund) ## looks better
invert_proxy_FD$log_stability <- log(invert_proxy_FD$stability)
hist(invert_proxy_FD$log_stability) ## looks better
### use the log transformed variables 

## merge species richness with FD/proxy data
invert_proxy_FD <- merge(invert_proxy_FD, spp_richness, by="GRIDREF")
## correlation between spp richness and FD metrics
cor1 <- cor.test(invert_proxy_FD$FDis_effect, invert_proxy_FD$average_spp_richness, method="pearson")
cor1 ## -0.0009, p=0.99
cor2 <- cor.test(invert_proxy_FD$FDis_response, invert_proxy_FD$average_spp_richness, method="pearson")
cor2 ## 0.26, p<0.001
cor3 <- cor.test(invert_proxy_FD$FDis_all, invert_proxy_FD$average_spp_richness, method="pearson")
cor3 ## 0.06, p=0.34
## low correlation between FDis and species richness
cor5 <- cor.test(invert_proxy_FD$average_spp_richness, invert_proxy_FD$log_abund, method="pearson")
cor5 ## 0.58, p<0.001
cor5 <- cor.test(invert_proxy_FD$average_spp_richness, invert_proxy_FD$log_stability, method="pearson")
cor5 ## 0.34, p<0.001

cor6 <- cor.test(invert_proxy_FD$log_abund, invert_proxy_FD$log_stability, method="pearson")
cor6 ## 0.1, p-0.16

####################################################################
######################### MEAN PROXY ###############################
####################################################################

#### FDis METRIC

## Effect traits only
mean_effect <- lm(log_abund ~ FDis_effect, data=invert_proxy_FD)
summary(mean_effect)
## non-significant
results_table1 <- data.frame(summary(mean_effect)$coefficients[,1:4])
results_table1$r_sq <- summary(mean_effect)$r.squared
write.csv(results_table1, file = "../Results/Pest control/mean_FDis_effect_det_interpol.csv", row.names=TRUE)

## Response traits only
mean_response <- lm(log_abund ~ FDis_response, data=invert_proxy_FD)
summary(mean_response)
## non-significant
results_table2 <- data.frame(summary(mean_response)$coefficients[,1:4])
results_table2$r_sq <- summary(mean_response)$r.squared
write.csv(results_table2, file = "../Results/Pest control/mean_FDis_response_det_interpol.csv", row.names=TRUE)

## All traits
mean_all <- lm(log_abund ~ FDis_all, data=invert_proxy_FD)
summary(mean_all)
## non-significant
results_table3 <- data.frame(summary(mean_all)$coefficients[,1:4])
results_table3$r_sq <- summary(mean_all)$r.squared
write.csv(results_table3, file = "../Results/Pest control/mean_FDis_all_det_interpol.csv", row.names=TRUE)

#########################################################################
######################### STABILITY PROXY ###############################
#########################################################################

#### FDis METRIC

## Effect traits only
stability_effect <- lm(log_stability ~ FDis_effect, data=invert_proxy_FD)
summary(stability_effect)
## significant negative relationship with effect
results_table4 <- data.frame(summary(stability_effect)$coefficients[,1:4])
results_table4$r_sq <- summary(stability_effect)$r.squared
write.csv(results_table4, file = "../Results/Pest control/stability_FDis_effect_det_interpol.csv", row.names=TRUE)

## Response traits only
stability_response <- lm(log_stability ~ FDis_response, data=invert_proxy_FD)
summary(stability_response)
## non-significant
results_table5 <- data.frame(summary(stability_response)$coefficients[,1:4])
results_table5$r_sq <- summary(stability_response)$r.squared
write.csv(results_table5, file = "../Results/Pest control/stability_FDis_response_det_interpol.csv", row.names=TRUE)

## All traits
stability_all <- lm(log_stability ~ FDis_all, data=invert_proxy_FD)
summary(stability_all)
## significant negative
results_table6 <- data.frame(summary(stability_all)$coefficients[,1:4])
results_table6$r_sq <- summary(stability_all)$r.squared
write.csv(results_table6, file = "../Results/Pest control/stability_FDis_all_det_interpol.csv", row.names=TRUE)

## correlations between effect and response FDis
cor7 <- cor.test(invert_proxy_FD$FDis_effect, invert_proxy_FD$FDis_response, method="pearson")
cor7 # FDis effect vs response r=-0.16, p=0.058

######################################### PLOTS ########################################

## mean_effect: non-significant
mean_effect <- ggplot(invert_proxy_FD, aes(x = FDis_effect, y = mean_abund)) + 
  geom_point(colour="black", size=0.4) +
  #geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  #geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Mean abundance", x=expression("F"[DIS]*" effect traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
mean_effect

## mean_response: non-significant
mean_response <- ggplot(invert_proxy_FD, aes(x = FDis_response, y = mean_abund)) + 
  geom_point(colour="black", size=0.4) +
  #stat_smooth(method = "lm", col = "black") +
  labs(y="Mean abundance", x=expression("F"[DIS]*" response traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
mean_response

## mean_all: non-significant
mean_all <- ggplot(invert_proxy_FD, aes(x = FDis_all, y = mean_abund)) + 
  geom_point(colour="black", size=0.4) +
  #geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  #geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Mean abundance", x=expression("F"[DIS]*" all traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
mean_all

## stability_effect: significant negative
r2_predict <- predict(stability_effect,interval="confidence")
newdata_stab_effect <- cbind(data.frame(invert_proxy_FD), data.frame(r2_predict))

newdata_stab_effect$fit <- exp(newdata_stab_effect$fit)
newdata_stab_effect$upr <- exp(newdata_stab_effect$upr)
newdata_stab_effect$lwr <- exp(newdata_stab_effect$lwr)

stability_effect <- ggplot(newdata_stab_effect, aes(x = FDis_effect, y = stability)) + 
  geom_point(colour="black", size=0.4) +
  geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Stability of abundance", x=expression("F"[DIS]*" effect traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
stability_effect

## stability_response: non-significant

stability_response <- ggplot(invert_proxy_FD, aes(x = FDis_response, y = stability)) + 
  geom_point(colour="black", size=0.4) +
  #geom_ribbon( aes(ymin = lwr, ymax = upr), fill="black", alpha = .15) +
  #geom_line(aes(y = fit), size = 0.4, color="black") +
  labs(y="Stability of abundance", x=expression("F"[DIS]*" response traits")) + 
  theme_classic() +
  theme(text = element_text(size=7))
stability_response

## stability_all: significant negative
r2_predict <- predict(stability_all,interval="confidence")
newdata_stab_all <- cbind(data.frame(invert_proxy_FD), data.frame(r2_predict))

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

## save all 6 plots as figure 2
fig2 <- plot_grid(mean_effect, mean_all, mean_response, stability_effect, stability_all, stability_response, 
          labels=c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"), label_size=8, ncol = 3, nrow = 2, 
          vjust =c(1,1,1,1,1,1), hjust=0)
fig2
ggsave(file="../Graphs/Figure2_invert_det_interpol.png", fig2, width = 120, height = 80, dpi = 600, units = "mm", device='png') 

bes2 <- plot_grid(mean_effect, mean_all, mean_response, stability_effect, stability_all, stability_response, ncol = 3, nrow = 2)
ggsave(file="../Graphs/BES_invert.png", bes2, height=7, width=10) 
