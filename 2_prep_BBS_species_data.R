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
seed_spp_list <- read.csv("../Data/Trait data/seed_dispersing_spp_list_JAT.csv", header=TRUE)
invert_spp_list <- read.csv("../Data/Trait data/insectivores_spp_list_JAT.csv", header=TRUE)

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

## calculate sum across all distance categories
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
## take max value across both visits
BBS_final <- BBS_final %>% group_by(YEAR, GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(TOT = max(TOT))
BBS_final <- as.data.frame(BBS_final)
head(BBS_final)

## remove mammals from dataset (they have numbers instead of characters as code)
BBS_final <- subset(BBS_final, !grepl("[0-9]", CBC_CODE))
## remove rows with zero TOT counts (this removes mammals and 1 house martin record and 10 rook records)
BBS_final<-BBS_final[!(BBS_final$TOT==0),]

## only include species with at least 12 years of data (non-consecutive)
BBS_final2 <- BBS_final %>% group_by(GRIDREF, ENGLISH_NAME, CBC_CODE) %>% summarise(no_years = n())
## remove site-species combos with less than 12 years of data
BBS_final2 <- BBS_final2[BBS_final2$no_years>=12,]
length(unique(BBS_final2$ENGLISH_NAME)) ## 110 species with at least 12 years of data
length(unique(BBS_final2$GRIDREF)) ## 200 sites still
table(BBS_final2$GRIDREF)

## merge species and sites into BBS_final
BBS_final3 <- merge(BBS_final2, BBS_final, by=c("ENGLISH_NAME","CBC_CODE", "GRIDREF"), all=FALSE)
length(unique(BBS_final3$ENGLISH_NAME)) ## 110 species
length(unique(BBS_final3$GRIDREF)) ## 200 sites still

## change some common names to merge with trait species list
BBS_final3$ENGLISH_NAME[ BBS_final3$ENGLISH_NAME == "Common Crossbill" ] <- "Crossbill"
BBS_final3$ENGLISH_NAME[ BBS_final3$ENGLISH_NAME == "Feral Pigeon" ] <- "Rock Dove"
BBS_final3$ENGLISH_NAME[ BBS_final3$ENGLISH_NAME == "Cetti's Warbler" ] <- "Cetti Warbler" ## apostrophe's don't merge


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

### take random value from DIFF to fill NAs from missing years
BBS_final6 <- BBS_final5 %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  mutate(DIFF=ifelse(is.na(DIFF), sample(na.omit(DIFF)), DIFF))

## yes/no for missing data
BBS_final6$missing_data <- ifelse(is.na(BBS_final6$TOT), "yes", "no")

## calculate new TOT using randomly sampled differences
## first descending order to fill in early years 
## run this multiple times
BBS_final7 <- BBS_final7 %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  arrange(desc(YEAR)) %>%
  mutate(TOT = ifelse(is.na(TOT), (lag(TOT) - lag(DIFF)), TOT))

## then run in increasing order to fill in late years
## run multiple times
BBS_final7 <- BBS_final7 %>% 
  group_by(ENGLISH_NAME, GRIDREF) %>%
  arrange(YEAR) %>%
  mutate(TOT = ifelse(is.na(TOT), (lag(TOT) + DIFF), TOT))
sum(is.na(BBS_final7$TOT))

## change any negative interpolated TOT values into zeros
BBS_final7$TOT <- ifelse(BBS_final7$TOT<0, 0, BBS_final7$TOT)

## save dataframe
write.csv(BBS_final7, file="../Data/Analysis_data/BBS_2004_2018_interpol.csv", row.names=FALSE) ## 200 sites and 110 species from 2004-2018

## plot 3 species with missing data and highlight years where missing
BBS_final <- read.csv("../Data/Analysis_data/BBS_2004_2018_interpol.csv", header=TRUE)

## calculate number of "no" i.e. years of missing data for each species at each site
BBS_missing_dat <- BBS_final %>%
  group_by(ENGLISH_NAME, GRIDREF) %>%
  summarise(missing_data = sum(missing_data == "yes"))

## plot avocet @ TG4221 - 3 missing datapoints (2004-2006)
avocet <- subset(BBS_final, ENGLISH_NAME=="Avocet" & GRIDREF=="TG4221")
avocet_plot <- ggplot(avocet, aes(x = YEAR, y = TOT)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(avocet, YEAR<=2006),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Avocet") +
  ylim(0,25) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
avocet_plot

## plot bittern @ TG4221 - 3 missing datapoints (2005, 2007 and 2010)
bittern <- subset(BBS_final, ENGLISH_NAME=="Bittern" & GRIDREF=="TG4221")
bittern_plot <- ggplot(bittern, aes(x = YEAR, y = TOT)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(bittern, YEAR==2005 | YEAR==2007 | YEAR==2010),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Bittern") +
  ylim(0,4) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
bittern_plot

## plot Whitethroat @ SD2707 - 3 missing datapoints (2016, 2017 and 2018)
whitethroat <- subset(BBS_final, ENGLISH_NAME=="Whitethroat" & GRIDREF=="SD2707")
whitethroat_plot <- ggplot(whitethroat, aes(x = YEAR, y = TOT)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(whitethroat, YEAR>=2016),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Whitethroat") +
  ylim(2,16) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
whitethroat_plot

## plot cetti warbler @ TG4221 - 2 missing datapoints (2013 and 2018)
cetti <- subset(BBS_final, ENGLISH_NAME=="Cetti Warbler" & GRIDREF=="TG4221")
cetti_plot <- ggplot(cetti, aes(x = YEAR, y = TOT)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(cetti, YEAR==2013 | YEAR== 2018),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Cetti's Warbler") +
  ylim(1,4) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
cetti_plot

## plot long-tailed tit @ SK3639 - 2 missing datapoints (2006 and 2009)
long_tit <- subset(BBS_final, ENGLISH_NAME=="Long-tailed Tit" & GRIDREF=="SK3639")
long_tit_plot <- ggplot(long_tit, aes(x = YEAR, y = TOT)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(long_tit, YEAR==2006 | YEAR==2009),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Long-tailed Tit") +
  ylim(0,5) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
long_tit_plot

## plot yellowhammer @ SJ6878 - 1 missing datapoint (2013)
yellow <- subset(BBS_final, ENGLISH_NAME=="Yellowhammer" & GRIDREF=="SJ6878")
yellow_plot <- ggplot(yellow, aes(x = YEAR, y = TOT)) +
  stat_smooth(colour="black", method=loess, se=FALSE) +
  geom_point(size=3, colour="black") + 
  geom_point(data = filter(yellow, YEAR==2013),
             aes(YEAR, TOT), colour = "red", size=3) +
  labs(x = "Year", y = "Abundance") +
  ggtitle("Yellowhammer") +
  ylim(0,6) +
  scale_x_continuous(breaks=seq(2004,2018,2)) +
  theme_classic() +
  theme(text = element_text(size = 11), legend.position="none") +
  labs(size=3)
yellow_plot

## plot these on one page
plots <- arrangeGrob(avocet_plot, bittern_plot, whitethroat_plot, cetti_plot, long_tit_plot, yellow_plot, nrow=3)
ggsave(file="../Graphs/interpolated_abundance.png", plots, width=9, height=12)


## merge seed dispersing species list (so we only have BBS data for seed dispersers)
seed_spp_list <- seed_spp_list[,-c(2,3,5)]
seed_spp_list <- seed_spp_list[seed_spp_list$yes_no==1,] ## 53 species which are seed dispersers
BBS_final4 <- merge(BBS_final3, seed_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_final4$ENGLISH_NAME)) ## 40 species which have BBS data

## merge invertivore species list (so we only have BBS data for pest control species)
invert_spp_list <- invert_spp_list[,-c(2,4)]
invert_spp_list <- invert_spp_list[invert_spp_list$yes_no==1,] ## 108 species which are pest controllers
invert_spp_list$species <- as.character(invert_spp_list$species)
invert_spp_list$species[ invert_spp_list$species == "Cetti’s Warbler" ] <- "Cetti Warbler" ## won't merge with apostrophe
BBS_final5 <- merge(BBS_final3, invert_spp_list, by.x="ENGLISH_NAME", by.y="species")
length(unique(BBS_final5$ENGLISH_NAME)) ## 66 species which have BBS data

## save seed dispersal files
write.csv(BBS_final4, file="../Data/Analysis_data/Seed dispersal/BBS_2004_2018_seed.csv", row.names=FALSE) ## 200 sites and 43 species from 2004-2018
spp_list2 <- BBS_final4[,c(1,2)] 
spp_list2 <- unique(spp_list2)
## save species list (43 species)
write.csv(spp_list2, file="../Data/Analysis_data/Seed dispersal/BBS_seed_species_list.csv", row.names=FALSE) ## 43 species

## save invertivore files
write.csv(BBS_final5, file="../Data/Analysis_data/Pest control/BBS_2004_2018_invert.csv", row.names=FALSE) ## 200 sites and 43 species from 2004-2018
spp_list3 <- BBS_final5[,c(1,2)] 
spp_list3 <- unique(spp_list3)
## save species list (73 species)
write.csv(spp_list3, file="../Data/Analysis_data/Pest control/BBS_invert_species_list.csv", row.names=FALSE) ## 73 species
