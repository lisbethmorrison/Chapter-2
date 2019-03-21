###################################################################################
## Title: Calculate species richness at each site & investigate confounding effect 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: January 2019
###################################################################################

rm(list=ls()) # clear R
options(scipen=999)

## packages
library(reshape2)
library(plyr)
library(dplyr)
library(data.table)

##############################################
## calculate species stability at each site ##
##############################################

## read in data
bbs <- read.csv("../../Chapter 1/Data/Bird_sync_data/bbs_woodland_birds.csv", header=TRUE)
bbs <- unique(bbs) # drop duplicate rows

bbs$site_code <- as.factor(bbs$site_code)
bbs$species_code <- as.factor(bbs$species_code)

## sort site code and year in ascending order
bbs <- bbs %>% arrange(site_code, ENGLISH_NAME, YEAR)

## Filter 1: only include sites with at least 8 years of continuous data
## i.e. 1994-2001 

bbs <- bbs %>%
  group_by(species_code, site_code) %>%
  mutate(Diff = YEAR - lag(YEAR))

## NAs = first value so need to change this to 1
bbs$Diff[is.na(bbs$Diff)] <- 1

#################################################################
######## code which puts consecutive years into groups ##########
#################################################################

start_time <- Sys.time()
i <- 1
site.list <- unique(bbs$site_code)
bbs_final <- NULL

for (j in site.list){
  print(j)
  bbs_1 <- bbs[bbs$site_code==j,] ## bbs_1 is the site dataframe
  
  spp.list <- unique(bbs_1$species_code) ## get unique species from site of interest
  
  for (k in spp.list){
    
    bbs_2 <- bbs_1[bbs_1$species_code==k,] ## bbs_2 is site and species combo of interest
    
      ## create findstop function which replaces values which are not 1 with NAs, and after NA's, add one to value
      ## this is so each cumulative run of 1's in Diff column are separated into blocks of different numbers
      findstop <-function(dt) {
        out <- NULL
        if(dt == 1) { out <- i}
        if(dt != 1) {
                     out <- "NA" 
                     i <<- i + 1}
        return(out)
      }
      
      i <- i+1
      
      bbs_2$cumstop  <-sapply(bbs_2$Diff, findstop) ## apply function to Diff column
      
      bbs_2$takeone <- lead(bbs_2$YEAR) ## puts the value below each year into separate column 
      bbs_2$takecumstop <- lead(bbs_2$cumstop) ## same for cumstop
 
    ## ifelse function - when cumstop does not equal NA, put cumstop value (i.e. do nothing)
    ## but when cumstop value does not equal NA AND Year == takeone -1, put takecumstop value, if not, put NA
      bbs_2<-mutate(bbs_2, fixeddiff = ifelse(cumstop != "NA",cumstop ,ifelse(YEAR == (takeone - 1), takecumstop,"NA")  )    )
    ## this part ensures that when there is an NA value which is actually consecutive, it will change it to consecutive
 
    bbs_temp <- data.frame(bbs_2$species_code, bbs_2$site_code, bbs_2$YEAR, bbs_2$ENGLISH_NAME, bbs_2$TOT, bbs_2$Diff, bbs_2$fixeddiff)
    bbs_final <- rbind(bbs_final, bbs_temp)
    
     }
} ## close loops

end_time <- Sys.time()
end_time - start_time ## 1.8 hours

## change colnames
colnames(bbs_final) <- c("species_code", "site_code", "YEAR", "ENGLISH_NAME", "TOT", "Diff", "Cum_years")
## order dataframe
bbs_final <- bbs_final %>% arrange(site_code, ENGLISH_NAME, YEAR)

## save dataframe
write.csv(bbs_final, file = "../Data/BBS sites/bbs_final_1.csv", row.names=FALSE)
## read data
#bbs_final <- read.csv("..//Data/BBS sites/bbs_final_1.csv", header=TRUE)

## now can remove NA values (as consecutive years are grouped, so NAs aren't needed)
bbs_final[bbs_final=="NA"] <- NA ## change character NAs to NA
bbs_final <- na.omit(bbs_final)

## group by Cum_years column, and use ifelse to put "yes" where length >=8, and "no" if not
bbs_final$Cum_years <- as.factor(bbs_final$Cum_years)
bbs_final <- bbs_final %>%
  group_by(Cum_years) %>%
  mutate(cons_years = ifelse(length(Cum_years)>=8, "yes", "no"))

## new dataframe with only "yes" values
bbs_final2 <- bbs_final[bbs_final$cons_years=="yes",]
length(unique(bbs_final2$site_code)) ## 2502 sites
length(unique(bbs_final2$species_code)) ## still 31 species  

## Filter 2: only include sites with a mean count of at least 9 individuals per year
## average TOT per site for each year
df <- bbs_final2 %>%
  group_by(site_code, YEAR) %>%
  summarise(Mean = mean(TOT))
## now remove sites which have any year with less than 9 individuals
df <- df %>%
  group_by(site_code) %>%
  filter(!any(Mean<9))
length(unique(df$site_code)) ## 116 sites...
## possibly too harsh a filter??

## Filter 3: only include species with data from at least 10 sites
df2 <- bbs_final2[,c(2,1)]
df2 <- unique(df2)
## calculate number of sites for each species
df2 <- df2 %>%
  group_by(species_code) %>%
  summarise(count = n())
## now remove species which have less than 10 sites
df2 <- df2[df2$count>=10,]
## removes 5 species (24 total)

## df = sites which have at least 9 individuals per year
## df2 = species which have data from at least 10 sites

## combine bbs_final2 with df2 == just removes 5 species from main dataframe
## leave df for now (takes out too many sites, maybe only suitable for UKBMS data?)

bbs_final2 <- merge(bbs_final2, df2, by="species_code", all=FALSE)
length(unique(bbs_final2$species_code)) ## 26 species

## save file 
write.csv(bbs_final2, file = "../Data/BBS sites/bbs_final_2.csv", row.names=FALSE)
#bbs_final2 <- read.csv("../Data/BBS sites/bbs_final_2.csv", header=TRUE)

###########################################################################################################
###########################################################################################################
###########################################################################################################

## Filters done - now calculate stability of each species, at each site over time
## stability = (1/CV)
## CV = SD/mean

## first calculate mean and SD
bbs_stability <- NULL
site.list <- unique(bbs_final2$site_code)
for (j in site.list){
  bbs_1 <- bbs_final2[bbs_final2$site_code==j,] ## bbs_1 is the site dataframe
  
  spp.list <- unique(bbs_1$species_code) ## get unique species from site of interest
  
  for (k in spp.list){
    
    bbs_2 <- bbs_1[bbs_1$species_code==k,] ## bbs_2 is site and species combo of interest
    
    mean <- mean(bbs_2$TOT)
    sd <- sd(bbs_2$TOT)
    species_code <- k
    site_code <- j
    results.temp<-data.frame(species_code,site_code,mean,sd)
    bbs_stability <- rbind(results.temp, bbs_stability)
  }
}

## calculate CV (mean/SD) and stability (1/CV)
bbs_stability$CV <- bbs_stability$sd/bbs_stability$mean
bbs_stability$stability <- 1/bbs_stability$CV
## 17 rows which have no varaibility (mean = 1, sd=0) so cannot calculate CV or stability
## remove these rows
bbs_stability <- bbs_stability[!bbs_stability$CV==0,] ## removes site 1677 (2501 total sites now)
## now take the average stability for each site (comparable to species richness)
av_stability <- aggregate(bbs_stability[, 6], list(bbs_stability$site_code), mean)
## change colnames
colnames(av_stability) <- c("site_code", "stability")

## save file
write.csv(av_stability, file = "../Data/BBS sites/bbs_average_stability.csv", row.names=FALSE)

#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

######### Now calculate species richness using same data (bbs_final2)
## remove zero counts
bbs_final2 <- bbs_final2[!bbs_final2$TOT=="0",]
str(bbs_final2)
bbs_final2$site_code <- as.factor(bbs_final2$site_code)
bbs_final2$species_code <- as.factor(bbs_final2$species_code)
bbs_final2$YEAR <- as.factor(bbs_final2$YEAR)

## add column which assigns each species a number (1)
## therefore sum of that column == species richness
bbs_final2$sppnumber <- 1

richness <- tapply(bbs_final2$sppnumber, bbs_final2[, c("site_code", "YEAR")], sum)
richness <- na.omit(melt(richness)) ## melt into 3 columns (site, year, richness)
## rename columns
colnames(richness) <- c("site_code", "year", "spp_richness")

## take the average richness over time
## gives average richness at each site - comparable to stability
av_richness <- aggregate(richness[, 3], list(richness$site_code), mean)
## change colnames
colnames(av_richness) <- c("site_code", "spp_richness")

#######################################################################################

## merge average stability and average richness together
richness_stability <- merge(av_richness, av_stability, all=FALSE)

model <- lm(spp_richness ~ stability, data=richness_stability)
summary(model)
## very significant (positive)
## as species richness increases, stability increases

plot(richness_stability$spp_richness, richness_stability$stability)
## positive relationship

