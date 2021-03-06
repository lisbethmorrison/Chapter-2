###################################################################################
## Title: Prep BBS site data for Chapter 2 analysis
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: March 2019
###################################################################################

rm(list=ls()) # clear R

library(dplyr)
library(Hmisc)
library(ggplot2)

## read in site data file
site_data <- read.csv("../Data/BBS_site_data/BBS_species_richness_site_year.csv", header=TRUE)
load("../../Chapter 1/Chapter1_script/Bird_scripts/gr_let2num.Rdata") # load function to convert gridrefs to easting and northing

## remove year and species richness columns
site_data <- site_data[-c(1,3)]
## remove duplicated sites
site_data <- unique(site_data)
site_data$gridref <- as.character(site_data$gridref)

## remove rows with gridrefs starting with 'I' (these are Irish sites)
site_data <- site_data[!grepl("I", site_data$gridref),]
site_data <- as.data.frame(site_data)
colnames(site_data) <- "gridref"

## convert gridrefs to Easting and Northing
easting_northing_values <- gr_let2num(site_data$gridref)
site_data <- cbind(site_data, easting_northing_values)

## export csv file
write.csv(site_data, file="../Data/BBS_site_data/BBS_East_North.csv", row.names=FALSE)

######################################################################################################################
######################################################################################################################
#################### Filter 1: Only include sites within one Bioclimatic Zone (Atlantic Central) ##################### 
######################################################################################################################
######################################################################################################################

rm(list=ls()) # clear R

## once ArcGIS analysis is done, read in site data again 
## these are the sites which are within one bioclimatic zone
bbs_sites_final <- read.csv("../Data/BBS_site_data/BBS_sites_final.csv", header=TRUE, stringsAsFactors=FALSE) # 3608 sites
## read in site_data file again (so we have species richness data)
site_data <- read.csv("../Data/BBS_site_data/BBS_species_richness_site_year.csv", header=TRUE, stringsAsFactors=FALSE)

## merge these site_data files together
site_data_final <- merge(site_data, bbs_sites_final, by="gridref", all=FALSE)
length(unique(site_data_final$gridref)) # 3608
## remove FID column (automatically generated by arcmap)
site_data_final$FID <- NULL

## group by gridref to find out how many years of data for each site
site_data_final_summ <- site_data_final %>% 
  group_by(gridref) %>% 
  summarise(freq = n()) ## nrow of each site 
## remove sites with less than 15 years of data
site_data_final_summ <- site_data_final_summ[site_data_final_summ$freq>=15,] # 1353 sites
## merge back into site_data_final
site_data_final <- merge(site_data_final, site_data_final_summ, by="gridref", all=FALSE)
length(unique(site_data_final$gridref)) # 1353 sites

## Filter 2: Only include sites with at least 15 years of consecutive data
site_data_final <- site_data_final %>% arrange(gridref, year)
site_data_final <- site_data_final %>%
  group_by(gridref) %>%
  mutate(Diff = year - lag(year))

## NAs = first value so need to change this to 1
site_data_final$Diff[is.na(site_data_final$Diff)] <- 1
str(site_data_final)

######### find sites which have at least 15 years of consecutive data
start_time <- Sys.time()
i <- 1
site.list <- unique(site_data_final$gridref)
site_data_final2 <- NULL

for (j in site.list){
  print(j)
  site_data_1 <- site_data_final[site_data_final$gridref==j,] ## site_data_1 is the site dataframe
  
  # spp.list <- unique(bbs_1$species_code) ## get unique species from site of interest
  # 
  # for (k in spp.list){
  #   
  #   bbs_2 <- bbs_1[bbs_1$species_code==k,] ## bbs_2 is site and species combo of interest
    
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
    
    site_data_1$cumstop  <-sapply(site_data_1$Diff, findstop) ## apply function to Diff column
    
    site_data_1$takeone <- lead(site_data_1$year) ## puts the value below each year into separate column 
    site_data_1$takecumstop <- lead(site_data_1$cumstop) ## same for cumstop
    
    ## ifelse function - when cumstop does not equal NA, put cumstop value (i.e. do nothing)
    ## but when cumstop value does not equal NA AND Year == takeone -1, put takecumstop value, if not, put NA
    site_data_1<-mutate(site_data_1, fixeddiff = ifelse(cumstop != "NA",cumstop ,ifelse(year == (takeone - 1), takecumstop,"NA")  )    )
    ## this part ensures that when there is an NA value which is actually consecutive, it will change it to consecutive
    
    site_temp <- data.frame(site_data_1$gridref, site_data_1$year, site_data_1$spp_richness, site_data_1$EASTING, 
                            site_data_1$NORTHING, site_data_1$freq, site_data_1$Diff, site_data_1$fixeddiff, stringsAsFactors = FALSE)
    site_data_final2 <- rbind(site_data_final2, site_temp)
    
    i <- i+1
  
  } ## close loops

end_time <- Sys.time()
end_time - start_time ## 8.4 seconds

## change colnames
colnames(site_data_final2) <- c("gridref", "year", "spp_richness", "EASTING", "NORTHING", "freq", "Diff", "Cum_years")
## order dataframe
site_data_final2 <- site_data_final2 %>% arrange(gridref, year)
## now can remove NA values (as consecutive years are grouped, so NAs aren't needed)
site_data_final2[site_data_final2=="NA"] <- NA ## change character NAs to NA
site_data_final2 <- na.omit(site_data_final2)
## group by Cum_years column, and use ifelse to put "yes" where length >=15, and "no" if not
site_data_final2$Cum_years <- as.factor(site_data_final2$Cum_years)
site_data_final2 <- site_data_final2 %>%
  group_by(Cum_years) %>%
  mutate(cons_years = ifelse(length(Cum_years)>=15, "yes", "no"))
## new dataframe with only "yes" values
site_data_final2 <- site_data_final2[site_data_final2$cons_years=="yes",]
length(unique(site_data_final2$gridref)) ## 672 sites

## loop through each 15 year period from 1994-2018, and calculate how many sites have data for each 15 year period

year.list<-1994:2004
num_sites_final <- NULL

for (i in year.list){ # loop through years
  start.year<-i
  end.year<-i+14

  num_sites <- site_data_final2 %>% group_by(gridref) %>% filter(min(year) <= start.year & max(year) >= end.year)
  
  if (dim(num_sites)[1] == 0) {
    next
  }
  

  num_sites2 <- as.data.frame(unique(num_sites$gridref))
  num_sites2$time_period <- i
  colnames(num_sites2)[1] <- "gridref"
  num_sites_final <- rbind(num_sites_final, num_sites2)
}

length(unique(num_sites_final$time_period)) ## 11 time periods were run 
length(unique(num_sites_final$gridref)) ## 672 sites in total
## how many rows per 15-year time period
num_sites_summ <- num_sites_final %>% group_by(time_period) %>% summarise(freq = n())
## time period with most sites is 2004-2018 (516 sites)
## filter these sites
num_sites_final <- num_sites_final[num_sites_final$time_period==2004,] # 516 rows
## merge this in with site_data_final2
site_data_final2 <- merge(site_data_final2, num_sites_final, by="gridref", all=FALSE)
length(unique(site_data_final2$gridref)) ## 516 sites
## now remove years <2004 (as we are only interested in time period 2004-2018)
site_data_final2 <- site_data_final2[!site_data_final2$year<2004,]
## check each year has 15 years of data
site_data_final_summ2 <- site_data_final2 %>% group_by(gridref) %>% summarise(freq = n()) ## Yes!

## Filter 3: Choose a random 20 sites from 10 categories of species richness

## take average richness over time for each site
av_richness <- aggregate(site_data_final2[, 3], list(site_data_final2$gridref), mean)
## change colnames
colnames(av_richness) <- c("gridref", "av_spp_richness")

## split species richness in 10 equal sized categories
av_richness$spp_richness_cat <- as.numeric(cut_number(av_richness$av_spp_richness,10))
## number of rows per group (we want same number per group)
nrow_group <- av_richness %>% 
  group_by(spp_richness_cat) %>% 
  summarise(freq = n()) ## nrow of each group 
## nearly equal number of observations (varies between 51 and 52 sites per group)
av_richness2 <- av_richness %>% group_by(spp_richness_cat) %>% sample_n(20) ## take random sample of 20 sites per group
## check that this produces 20 sites per group
nrow_group2 <- av_richness2 %>% 
  group_by(spp_richness_cat) %>% 
  summarise(freq = n()) ## nrow of each group (20 in each group)
## 200 sites in total

## merge these sites back into site_data_final2
site_data_final3 <- merge(av_richness2, site_data_final2, by="gridref", all=FALSE)
length(unique(site_data_final3$gridref)) ## 200 sites!!

## remove uneccessary columns
site_data_final3 <- site_data_final3[-c(8:12)]

## save file
write.csv(site_data_final3, file="../Data/BBS_site_data/BBS_filtered_sites_final.csv", row.names=FALSE)
  
  
  
  


