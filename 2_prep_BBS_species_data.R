###########################################################
## Title: Prep BBS species data for analysis 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2019
##########################################################

rm(list=ls()) # clear R

#install.packages("dplyr", dep = TRUE)
library(dplyr)
library(tidyr)
library(ggplot2)

## read in data
BBS_2018 <- read.csv("../Data/BBS_abund_data/BBS_data_2018.csv", header=TRUE)
BBS_2004_2017 <- read.csv("../Data/BBS_abund_data/BBS_data_2004_2017.csv", header=TRUE)
seed_spp_list <- read.csv("../Data/Trait data/seed_dispersing_spp_list.csv", header=TRUE)
invert_spp_list <- read.csv("../Data/Trait data/insectivores_spp_list.csv", header=TRUE)

BBS_2018$YEAR <- as.factor(BBS_2018$YEAR)
### first calculate TOT for 2004-2017 data
## this is summed across all distance categories and transect sections for each year and visit
## first change NAs to zeros (for some reason zeros have changed to NAs)
BBS_2004_2017[is.na(BBS_2004_2017)] <- 0

str(BBS_2004_2017)
BBS_2004_2017$ï..CBC_CODE <- as.character(BBS_2004_2017$ï..CBC_CODE)
BBS_2004_2017$ENGLISH_NAME <- as.character(BBS_2004_2017$ENGLISH_NAME)
BBS_2004_2017$YEAR <- as.factor(BBS_2004_2017$YEAR)
BBS_2004_2017$BAND_1 <- as.numeric(BBS_2004_2017$BAND_1)
BBS_2004_2017$BAND_2 <- as.numeric(BBS_2004_2017$BAND_2)
BBS_2004_2017$BAND_3 <- as.numeric(BBS_2004_2017$BAND_3)
BBS_2004_2017$FLYING <- as.numeric(BBS_2004_2017$FLYING)

## calculate sum across all distance categories for each transect
BBS_2004_2017$SUM <- rowSums(BBS_2004_2017[,8:11])

## calculate TOT across all transects
BBS_2004_2017_tot <- BBS_2004_2017 %>%
  group_by(ï..CBC_CODE, ENGLISH_NAME, YEAR, OBS_DT, VISIT, GRIDREF) %>% 
  summarise(TOT = sum(SUM))
BBS_2004_2017_tot <- as.data.frame(BBS_2004_2017_tot)
head(BBS_2004_2017_tot)
## change first column to CBC_CODE
names(BBS_2004_2017_tot)[1]<-paste("CBC_CODE")
## remove some columns not needed (tenkm and region_code) from 2018 dataset
BBS_2018 <- BBS_2018[, -c(2:3)]
head(BBS_2018)

## merge together
BBS_final <- rbind(BBS_2004_2017_tot, BBS_2018)

## remove mammals from dataset (they have numbers instead of characters as code)
BBS_final <- subset(BBS_final, !grepl("[0-9]", CBC_CODE))
BBS_final <-  BBS_final[!(is.na(BBS_final$CBC_CODE) | BBS_final$CBC_CODE==""), ] ## remove blanks (sheep)

## change some common names to merge with trait & detectability species list
BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Common Crossbill" ] <- "Crossbill"
BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Feral Pigeon" ] <- "Rock Dove"
BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Cetti's Warbler" ] <- "Cetti Warbler" ## apostrophe's don't merge
BBS_final$ENGLISH_NAME[ BBS_final$ENGLISH_NAME == "Pied/White Wagtail" ] <- "Pied Wagtail" ## apostrophe's don't merge
## change rock dove CBC code to DV, not FP
BBS_final$CBC_CODE[BBS_final$ENGLISH_NAME=='Rock Dove'] <- 'DV'

## fill in data with zeros - if there is no record for a speciecs in a site-year combination, the species was not recorded = zero
BBS_final$YEAR <- as.numeric(levels(BBS_final$YEAR))[BBS_final$YEAR]
BBS_final2 <- BBS_final %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  complete(ENGLISH_NAME, CBC_CODE, VISIT, YEAR = 2004:2018)
BBS_final2$TOT[is.na(BBS_final2$TOT)] <- 0 ## change NAs to zero
## remove OBS_DATE column - not required
BBS_final2 <- BBS_final2[,-6]
## this is the complete dataset filled with zeros where a species wasn't seen at a site but the site WAS visited
length(unique(BBS_final2$ENGLISH_NAME)) ## 208 species
length(unique(BBS_final2$GRIDREF)) ## 200 sites


## merge seed dispersing species list (so we only have BBS data for seed dispersers)
seed_spp_list <- seed_spp_list[,-c(2,3,5)]
seed_spp_list <- seed_spp_list[seed_spp_list$yes_no==1,] ## 53 species which are seed dispersers
BBS_seed <- merge(BBS_final2, seed_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_seed$ENGLISH_NAME)) ## 49 species which have BBS data

## merge invertivore species list (so we only have BBS data for pest control species)
invert_spp_list <- invert_spp_list[,-c(2,4)]
invert_spp_list <- invert_spp_list[invert_spp_list$yes_no==1,] ## 100 species which are pest controllers
invert_spp_list$species <- as.character(invert_spp_list$species)
invert_spp_list$species[ invert_spp_list$species == "Cetti’s Warbler" ] <- "Cetti Warbler" ## won't merge with apostrophe
BBS_invert <- merge(BBS_final2, invert_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_invert$ENGLISH_NAME)) ## 89 species which have BBS data

BBS_final3 <- merge(BBS_seed, BBS_invert, all=TRUE)
BBS_final3 <- BBS_final3[,-7]
## save file to send to Dario
write.csv(BBS_final3, file="../Data/BBS_abund_data/BBS_2004_2018_zeros.csv", row.names=FALSE)
length(unique(BBS_final3$ENGLISH_NAME)) ## 108 unique species

seed_detect <- read.csv("../Data/BBS_abund_data/seed_dispersers_with_det_prob.csv", header=TRUE)
invert_detect <- read.csv("../Data/BBS_abund_data/insect_pred_with_det_prob.csv", header=TRUE)

bbs_detect <- rbind(seed_detect, invert_detect)
bbs_detect <- unique(bbs_detect)
length(unique(bbs_detect$ENGLISH_NAME)) ## 74 species
length(unique(bbs_detect$GRIDREF)) ## 200 sites


bbs_detect2 <- na.omit(bbs_detect)
length(unique(bbs_detect2$ENGLISH_NAME)) ## 72 species (chough and rock pipit don't have detectability data)
## some site-species combinations will also be removed (rook, stock dove, house martin, sand martin, swallow and swift)
length(unique(bbs_detect2$GRIDREF)) ## 200 sites
## Remove 2nd column (ENGLISH_NAME2)
bbs_detect2 <- bbs_detect2[,-2]
## BBS_final_detect1 = complete cases
## rename columns in bbs_detect2
colnames(bbs_detect2) <- c("GRIDREF", "ENGLISH_NAME", "CBC_CODE", "YEAR", "VISIT", "det.prob")
bbs_detect2$YEAR <- as.factor(bbs_detect2$YEAR)
BBS_presence$YEAR <- as.factor(BBS_presence$YEAR)
BBS_presence$CBC_CODE <- as.factor(BBS_presence$CBC_CODE)
BBS_presence$ENGLISH_NAME <- as.factor(BBS_presence$ENGLISH_NAME)

## merge two together (72 species, 200 sites)
BBS_detect <- merge(bbs_detect2, BBS_final3, by=c("ENGLISH_NAME", "CBC_CODE", "GRIDREF", "YEAR", "VISIT"), all=TRUE)

#### IRRELEVANT ==> ALL 208 SPECIES HAVE 15 YEARS OF DATA (BUT NOT ALL OF THAT DATA IS PRESENCE)
## only include species with at least 12 years of data (non-consecutive)
BBS_final3 <- BBS_final2 %>% group_by(GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(no_years = n())
## remove site-species combos with less than 12 years of data
BBS_final3 <- BBS_final3[BBS_final3$no_years>=12,]
length(unique(BBS_final3$ENGLISH_NAME)) ## 208 species with at least 12 years of data
length(unique(BBS_final3$GRIDREF)) ## 200 sites still
table(BBS_final3$GRIDREF)

## merge species and sites into BBS_final
BBS_final3 <- merge(BBS_final2, BBS_final, by=c("ENGLISH_NAME","CBC_CODE", "GRIDREF"), all=FALSE)
length(unique(BBS_final3$ENGLISH_NAME)) ## 110 species
length(unique(BBS_final3$GRIDREF)) ## 200 sites still
## 70,629 rows




######################################### FILL IN MISSING ABUNDANCE DATA ###############################################

## remove no.years column
BBS_final4 <- BBS_final3[,-4]

## complete data so each species/site combo has 15 years filled with NAs
BBS_final4$YEAR <- as.numeric(levels(BBS_final4$YEAR))[BBS_final4$YEAR]
BBS_final5 <- BBS_final4 %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  complete(ENGLISH_NAME, CBC_CODE, YEAR = 2004:2018)
############ WORKS!! ############ 73680 rows

## create two new columns: TOT_lag and DIFF
BBS_final5 <- BBS_final5 %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  arrange(YEAR) %>%
  mutate(TOT_lag=lag(TOT), DIFF=TOT-TOT_lag)

### take random value from DIFF to fill NAs from missing years - run multiple times? e.g. 100 for now
BBS_interpol <- NULL
for(i in 1:100){
  print(i)
  BBS_final6 <- BBS_final5 %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  mutate(DIFF=ifelse(is.na(DIFF) | is.na(TOT), sample(na.omit(DIFF)), DIFF))
  results_table_temp <- data.frame(BBS_final6$GRIDREF, BBS_final6$ENGLISH_NAME, BBS_final6$CBC_CODE, BBS_final6$YEAR, BBS_final6$TOT,
                                   BBS_final6$TOT_lag, BBS_final6$DIFF,i, stringsAsFactors = FALSE)
  BBS_interpol <- as.data.frame(rbind(BBS_interpol, results_table_temp))
}
colnames(BBS_interpol) <- c("GRIDREF", "ENGLISH_NAME", "CBC_CODE", "YEAR", "TOT", "TOT_lag", "DIFF", "i")

## yes/no for missing data
BBS_interpol$missing_data <- ifelse(is.na(BBS_interpol$TOT), "yes", "no")
BBS_interpol2 <- BBS_interpol
BBS_interpol2$i <- as.factor(BBS_interpol2$i)

## calculate new TOT using randomly sampled differences
## first descending order to fill in early years 
## run this multiple times
BBS_interpol3 <- BBS_interpol2 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(desc(YEAR)) %>%
  mutate(TOT = ifelse(is.na(TOT), (lag(TOT) - lag(DIFF)), TOT))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT <- ifelse(BBS_interpol3$TOT<0, 0, BBS_interpol3$TOT)
## run again
BBS_interpol3 <- BBS_interpol3 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(desc(YEAR)) %>%
  mutate(TOT = ifelse(is.na(TOT), (lag(TOT) - lag(DIFF)), TOT))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT <- ifelse(BBS_interpol3$TOT<0, 0, BBS_interpol3$TOT)
## run again
BBS_interpol3 <- BBS_interpol3 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(desc(YEAR)) %>%
  mutate(TOT = ifelse(is.na(TOT), (lag(TOT) - lag(DIFF)), TOT))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT <- ifelse(BBS_interpol3$TOT<0, 0, BBS_interpol3$TOT)
## then run in increasing order to fill in late years
## run multiple times
BBS_interpol3 <- BBS_interpol3 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(YEAR) %>%
  mutate(TOT = ifelse(is.na(TOT), (lag(TOT) + DIFF), TOT))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT <- ifelse(BBS_interpol3$TOT<0, 0, BBS_interpol3$TOT)
## run again
BBS_interpol3 <- BBS_interpol3 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(YEAR) %>%
  mutate(TOT = ifelse(is.na(TOT), (lag(TOT) + DIFF), TOT))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT <- ifelse(BBS_interpol3$TOT<0, 0, BBS_interpol3$TOT)
## run again
BBS_interpol3 <- BBS_interpol3 %>% 
  group_by(ENGLISH_NAME, GRIDREF, i) %>%
  arrange(YEAR) %>%
  mutate(TOT = ifelse(is.na(TOT), (lag(TOT) + DIFF), TOT))
## change any negative interpolated TOT values into zeros
BBS_interpol3$TOT <- ifelse(BBS_interpol3$TOT<0, 0, BBS_interpol3$TOT)

sum(is.na(BBS_interpol3$TOT)) ## zero NAs 

## save dataframe
write.csv(BBS_interpol3, file="../Data/BBS_abund_data/BBS_2004_2018_interpol_repeats.csv", row.names=FALSE) ## 200 sites and 110 species from 2004-2018
BBS_interpol <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_interpol_repeats.csv", header=TRUE)

## plot 6 species with missing data and highlight years where missing
BBS_final <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_interpol_repeats.csv", header=TRUE)

## plot avocet @ TG4221 - 3 missing datapoints (2004-2006)
avocet <- subset(BBS_final, ENGLISH_NAME=="Avocet" & GRIDREF=="TG4221")
avocet_plot <- ggplot(avocet, aes(x = YEAR, y = TOT, group=i)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(avocet, YEAR<=2006),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Avocet") +
  #ylim(0,25) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
avocet_plot

## plot rook @ SX8099 - 3 missing datapoints (2004, 2006, 2011)
rook <- subset(BBS_final, ENGLISH_NAME=="Rook" & GRIDREF=="SX8099")
rook_plot <- ggplot(rook, aes(x = YEAR, y = TOT, group=i)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(rook, YEAR==2004 | YEAR==2006 | YEAR==2011),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Rook") +
  #ylim(0,45) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
rook_plot

# ## plot bittern @ TG4221 - 3 missing datapoints (2005, 2007 and 2010)
# bittern <- subset(BBS_final, ENGLISH_NAME=="Bittern" & GRIDREF=="TG4221")
# bittern_plot <- ggplot(bittern, aes(x = YEAR, y = TOT)) +
#   stat_smooth(colour="black", method=loess, se=FALSE) +
#   geom_point(size=3, colour="black") + 
#   geom_point(data = filter(bittern, YEAR==2005 | YEAR==2007 | YEAR==2010),
#              aes(YEAR, TOT), colour = "red", size=3) +
#   labs(x = "Year", y = "Abundance") +
#   ggtitle("Bittern") +
#   ylim(0,4) +
#   scale_x_continuous(breaks=seq(2004,2018,2)) +
#   theme_classic() +
#   theme(text = element_text(size = 11), legend.position="none") +
#   labs(size=3)
# bittern_plot

## plot Whitethroat @ SD2707 - 3 missing datapoints (2016, 2017 and 2018)
whitethroat2 <- subset(BBS_final, subset=(YEAR>=2016 & missing_data=="yes"))

whitethroat <- subset(BBS_final, ENGLISH_NAME=="Whitethroat" & GRIDREF=="SD2707")
whitethroat_plot <- ggplot(whitethroat, aes(x = YEAR, y = TOT, group=i)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(whitethroat, YEAR>=2016),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Whitethroat") +
  #ylim(2,16) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
whitethroat_plot

## plot cetti warbler @ TG4221 - 2 missing datapoints (2013 and 2018)
cetti <- subset(BBS_final, ENGLISH_NAME=="Cetti Warbler" & GRIDREF=="TG4221")
cetti_plot <- ggplot(cetti, aes(x = YEAR, y = TOT, group=i)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(cetti, YEAR==2013 | YEAR== 2018),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Cetti's Warbler") +
  #ylim(1,4) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
cetti_plot

## plot long-tailed tit @ SK3639 - 2 missing datapoints (2006 and 2009)
long_tit <- subset(BBS_final, ENGLISH_NAME=="Long-tailed Tit" & GRIDREF=="SK3639")
long_tit_plot <- ggplot(long_tit, aes(x = YEAR, y = TOT, group=i)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(long_tit, YEAR==2006 | YEAR==2009),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Long-tailed Tit") +
  #ylim(0,5) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
long_tit_plot

## plot yellowhammer @ SJ6878 - 1 missing datapoint (2013)
yellow <- subset(BBS_final, ENGLISH_NAME=="Yellowhammer" & GRIDREF=="SJ6878")
yellow_plot <- ggplot(yellow, aes(x = YEAR, y = TOT, group=i)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(yellow, YEAR==2013),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Yellowhammer") +
  #ylim(0,6) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
yellow_plot

## plot these on one page
library(gridExtra)
plots <- arrangeGrob(avocet_plot, rook_plot, whitethroat_plot, cetti_plot, long_tit_plot, yellow_plot, nrow=2)
plots
ggsave(file="../Graphs/interpolated_abundance_average.png", plots, width=15, height=10)


## merge seed dispersing species list (so we only have BBS data for seed dispersers)
seed_spp_list <- seed_spp_list[,-c(2,3,5)]
seed_spp_list <- seed_spp_list[seed_spp_list$yes_no==1,] ## 53 species which are seed dispersers
BBS_seed <- merge(BBS_final, seed_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_seed$ENGLISH_NAME)) ## 40 species which have BBS data

## merge invertivore species list (so we only have BBS data for pest control species)
invert_spp_list <- invert_spp_list[,-c(2,4)]
invert_spp_list <- invert_spp_list[invert_spp_list$yes_no==1,] ## 100 species which are pest controllers
invert_spp_list$species <- as.character(invert_spp_list$species)
invert_spp_list$species[ invert_spp_list$species == "Cetti’s Warbler" ] <- "Cetti Warbler" ## won't merge with apostrophe
BBS_invert <- merge(BBS_final, invert_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_invert$ENGLISH_NAME)) ## 59 species which have BBS data

## save seed dispersal files
write.csv(BBS_seed, file="../Data/BBS_abund_data/BBS_2004_2018_seed_interpol_repeat.csv", row.names=FALSE) ## 200 sites and 40 species from 2004-2018
spp_list2 <- BBS_seed[,c(1,3)] 
spp_list2 <- unique(spp_list2)
## save species list (40 species)
write.csv(spp_list2, file="../Data/Trait data/Seed dispersal/BBS_seed_species_list.csv", row.names=FALSE) ## 40 species

## save insectivore files
write.csv(BBS_invert, file="../Data/BBS_abund_data/BBS_2004_2018_invert_interpol_repeat.csv", row.names=FALSE) ## 200 sites and 59 species from 2004-2018
spp_list3 <- BBS_invert[,c(1,3)] 
spp_list3 <- unique(spp_list3)
## save species list (59 species)
write.csv(spp_list3, file="../Data/Trait data/Pest control/BBS_invert_species_list.csv", row.names=FALSE) ## 59 species

BBS_seed <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_seed_interpol_repeat.csv", header=TRUE)
BBS_invert <- read.csv("../Data/BBS_abund_data/BBS_2004_2018_invert_interpol_repeat.csv", header=TRUE)
BBS <- rbind(BBS_seed, BBS_invert)
BBS <- unique(BBS)
length(unique(BBS$ENGLISH_NAME)) ## 75 species
## 244,300 rows of missing data = yes
## total dataset = 6,369,000
## % missing data = 3.8
BBS2 <- BBS[BBS$missing_data=="yes",] ## 100 species which are pest controllers
length(unique(BBS2$ENGLISH_NAME)) ## 69
