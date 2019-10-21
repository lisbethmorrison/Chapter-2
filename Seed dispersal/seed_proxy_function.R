###########################################################
## Title: Calculate mean and stability of seed dispersal function
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: June 2019
##########################################################

rm(list=ls()) # clear R
library(tidyverse)

## add data
BBS_final <- read.csv("../Data/Analysis_data/Seed dispersal/BBS_2004_2018_seed.csv", header=TRUE)
effect <- read.csv("../Data/Trait data/Seed dispersal/seed_effect_traits.csv", header=TRUE) ## 32 spp
response <- read.csv("../Data/Trait data/Seed dispersal/seed_response_traits.csv", header=TRUE) ## 38 spp
all_traits <- read.csv("../Data/Trait data/Seed dispersal/seed_all_traits.csv", header=TRUE) ## 28 spp
## read these files in to get species list

# ## remove unecessary columns from each trait data (only need name and code to merge with BBS data)
# ## only need english name and yes/no
# effect <- effect[,c(1,2)]
# response <- response[,c(1,2)]
# effect_response <- effect_response[,c(1,2)]

## merge each file with BBS_final abundance data
BBS_final_effect <- merge(BBS_final, effect, by=c("ENGLISH_NAME", "CBC_CODE"), all=FALSE)
length(unique(BBS_final_effect$ENGLISH_NAME)) ## 32 species
BBS_final_response <- merge(BBS_final, response, by=c("ENGLISH_NAME", "CBC_CODE"), all=FALSE)
length(unique(BBS_final_response$ENGLISH_NAME)) ## 38 species
BBS_final_all <- merge(BBS_final, all_traits, by=c("ENGLISH_NAME", "CBC_CODE"), all=FALSE)
length(unique(BBS_final_all$ENGLISH_NAME)) ## 28 species

#### try using partial pooling to estimate abundance for species with missing data
library(lme4)
library(dplyr)
library(tibble)
## visualise data for each species

## get subset of data on one site for now
bbs_data <- subset(BBS_final_effect, GRIDREF == "TR0642")
bbs_data <- bbs_data[,-c(2,4,7:12)]
bbs_data <- BBS_final_effect

bbs_data$ENGLISH_NAME <- as.character(bbs_data$ENGLISH_NAME)
bbs_data$YEAR <- as.numeric(bbs_data$YEAR)
bbs_data$TOT <- as.numeric(bbs_data$TOT)

ggplot(bbs_data) + 
  aes(x = YEAR, y = TOT) + 
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
  geom_point() +
  facet_wrap("ENGLISH_NAME")

#### NO POOLING #####
df_no_pooling <- lmList(TOT ~ YEAR | ENGLISH_NAME, bbs_data) %>% 
  coef() %>% 
  # Subject IDs are stored as row-names. Make them an explicit column
  rownames_to_column("ENGLISH_NAME") %>% 
  rename(Intercept = `(Intercept)`, Slope_year = YEAR) %>% 
  add_column(Model = "No pooling")

##### COMPLETE POOLING #####
m_pooled <- lm(TOT ~ YEAR, bbs_data) 

# Repeat the intercept and slope terms for each participant
df_pooled <- data_frame(
  Model = "Complete pooling",
  ENGLISH_NAME = unique(bbs_data$ENGLISH_NAME),
  Intercept = coef(m_pooled)[1], 
  Slope_year = coef(m_pooled)[2])

df_models <- bind_rows(df_pooled, df_no_pooling) %>% 
  left_join(bbs_data, by = "ENGLISH_NAME")

p_model_comparison <- ggplot(df_models) + 
  aes(x = YEAR, y = TOT) + 
  # Set the color mapping in this layer so the points don't get a color
  geom_abline(aes(intercept = Intercept, slope = Slope_year, color = Model),
              size = .75) + 
  geom_point() +
  facet_wrap("ENGLISH_NAME") +
  #labs(x = xlab, y = ylab) + 
  #scale_x_continuous(breaks = 0:4 * 2) + 
  # Fix the color palette 
  scale_color_brewer(palette = "Dark2") + 
  theme(legend.position = "top")

p_model_comparison

#### PARTIAL POOLING ####
m <- lmer(TOT ~ 1 + YEAR + (1+YEAR|ENGLISH_NAME), bbs_data) ## singular fit error
tt <- getME(m,"theta")
ll <- getME(m,"lower")
min(tt[ll==0])

m2 <- lmer(TOT ~ 1 + YEAR + (1+YEAR|ENGLISH_NAME), bbs_data, 
           control = lmerControl(optimizer = "optimx", calc.derivs = FALSE, 
                                 optCtrl = list(method = "nlminb", starttests = FALSE, kkt = FALSE)))

## try different optimizer
library(optimx)
ss <- getME(m,c("theta","fixef"))
m2 <- update(m,start=ss,control=lmerControl(optimizer="nlminb",
                                                 optCtrl=list(maxfun=2e5)))

relgrad <- with(m@optinfo$derivs,solve(Hessian,gradient))
max(abs(relgrad))

df_partial_pooling <- coef(m)[["ENGLISH_NAME"]] %>% 
  rownames_to_column("ENGLISH_NAME") %>% 
  as_tibble() %>% 
  rename(Intercept = `(Intercept)`, Slope_year = YEAR) %>% 
  add_column(Model = "Partial pooling")

df_models <- bind_rows(df_pooled, df_no_pooling, df_partial_pooling) %>% 
  left_join(bbs_data, by = "ENGLISH_NAME")

# Replace the data-set of the last plot
p_model_comparison %+% df_models

## zoom in on some key species (first 3 have missing data, last one doesn't)
df_zoom <- df_models %>% 
  filter(ENGLISH_NAME %in% c("Chaffinch", "Whimbrel", "Crossbill", "Crane"))
p_model_comparison %+% df_zoom



############################################################################################################
############################################################################################################
#### MEAN FUNCTION

# ## take the mean abundance for each species at each site over time
# BBS_mean_effect <- BBS_final_effect %>% group_by(ENGLISH_NAME, GRIDREF) %>% summarise(Mean_abund = mean(TOT)) 
# BBS_mean_response <- BBS_final_response %>% group_by(ENGLISH_NAME, GRIDREF) %>% summarise(Mean_abund = mean(TOT)) 
# BBS_mean_all <- BBS_final_all %>% group_by(ENGLISH_NAME, GRIDREF) %>% summarise(Mean_abund = mean(TOT)) 
# 
# ## calculate the sum of mean abundance at each site
# BBS_mean_effect <- BBS_mean_effect %>% group_by(GRIDREF) %>% summarise(Mean_abund = sum(Mean_abund)) ## 32 species at 199 sites
# BBS_mean_response <- BBS_mean_response %>% group_by(GRIDREF) %>% summarise(Mean_abund = sum(Mean_abund)) ## 32 species at 199 sites
# BBS_mean_all <- BBS_mean_all %>% group_by(GRIDREF) %>% summarise(Mean_abund = sum(Mean_abund)) ## 32 species at 199 sites
## this is done below anyway

############################################################################################################
############################################################################################################
#### STABILITY OF FUNCTION
## calculate mean and SD for each species and site

### effect species (=32)
## calc mean and SD abundance of each species at each site
BBS_tot_effect <- BBS_final_effect %>% group_by(GRIDREF,ENGLISH_NAME) %>% summarise(mean_abund = mean(TOT), sd=sd(TOT))
## the calculate variance of each species (squared SD)
BBS_tot_effect$variance <- (BBS_tot_effect$sd)^2
## calc total mean and SD (squareroot of sum of variance) at each site
BBS_proxy_effect <- BBS_tot_effect %>% group_by(GRIDREF) %>% summarise(sum_mean = sum(mean_abund), sum_var=sum(variance))
BBS_proxy_effect$sum_sd <- sqrt(BBS_proxy_effect$sum_var)
BBS_proxy_effect$CV <- BBS_proxy_effect$sum_sd/BBS_proxy_effect$sum_mean ## calculate CV (=SD/mean)
BBS_proxy_effect$stability <- 1/BBS_proxy_effect$CV ## calculate stability (1/CV)
## remove NA values (where there is only one species at a site - can't calc SD or stability)
BBS_proxy_effect <- na.omit(BBS_proxy_effect) ## 199 sites

## response species (=38)
BBS_tot_response <- BBS_final_response %>% group_by(GRIDREF,ENGLISH_NAME) %>% summarise(mean_abund = mean(TOT), sd=sd(TOT))
## the calculate variance of each species (squared SD)
BBS_tot_response$variance <- (BBS_tot_response$sd)^2
## calc total mean and SD (squareroot of sum of variance) at each site
BBS_proxy_response <- BBS_tot_response %>% group_by(GRIDREF) %>% summarise(sum_mean = sum(mean_abund), sum_var=sum(variance))
BBS_proxy_response$sum_sd <- sqrt(BBS_proxy_response$sum_var)
BBS_proxy_response$CV <- BBS_proxy_response$sum_sd/BBS_proxy_response$sum_mean ## calculate CV (=SD/mean)
BBS_proxy_response$stability <- 1/BBS_proxy_response$CV ## calculate stability (1/CV)
## remove NA values (where there is only one species at a site - can't calc SD or stability)
BBS_proxy_response <- na.omit(BBS_proxy_response) ## 200 sites

## all trait species (=28)
BBS_tot_all <- BBS_final_all %>% group_by(GRIDREF,ENGLISH_NAME) %>% summarise(mean_abund = mean(TOT), sd=sd(TOT))
## the calculate variance of each species (squared SD)
BBS_tot_all$variance <- (BBS_tot_all$sd)^2
## calc total mean and SD (squareroot of sum of variance) at each site
BBS_proxy_all <- BBS_tot_all %>% group_by(GRIDREF) %>% summarise(sum_mean = sum(mean_abund), sum_var=sum(variance))
BBS_proxy_all$sum_sd <- sqrt(BBS_proxy_all$sum_var)
BBS_proxy_all$CV <- BBS_proxy_all$sum_sd/BBS_proxy_all$sum_mean ## calculate CV (=SD/mean)
BBS_proxy_all$stability <- 1/BBS_proxy_all$CV ## calculate stability (1/CV)
## remove NA values (where there is only one species at a site - can't calc SD or stability)
BBS_proxy_all <- na.omit(BBS_proxy_all) ## 199 sites

## remove sum_var, sum_sd and CV columns from each
BBS_proxy_effect <- BBS_proxy_effect[,-c(3:5)]
BBS_proxy_response <- BBS_proxy_response[,-c(3:5)]
BBS_proxy_all <- BBS_proxy_all[,-c(3:5)]

## save files
write.csv(BBS_proxy_effect, file="../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_effect.csv", row.names=FALSE) ## 199 sites
write.csv(BBS_proxy_response, file="../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_response.csv", row.names=FALSE) ## 200 sites
write.csv(BBS_proxy_all, file="../Data/Analysis_data/Seed dispersal/BBS_proxy_seed_all.csv", row.names=FALSE) ## 199 sites


