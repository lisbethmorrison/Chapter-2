###########################################################
## Title: Calculate landscape heterogeneity for BBS sites 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: January 2019
##########################################################

## packages
library(vegan)

### add data
bbs_hab <- read.table("../Data/BBS sites/BBS.CBC.UKBMS.complete.landcover.data.all.scales.soil.dem.configuration.txt", header=TRUE)

## subset BBS data
bbs_hab <- bbs_hab[bbs_hab$Surv=="BBS",]

## remove columns not needed
bbs_hab <- bbs_hab[-c(3,6:16,30:85)]

## create total_count column
bbs_hab$total_count <- rowSums(bbs_hab[5:17])

# divide all columns by total land area
habs <- c("A","BgRo","Br","BW","C","CW","F","G","H","M","S","R","UG")

for (i in habs){
  bbs_hab[,i] <- bbs_hab[,i]/bbs_hab[,"total_count"]
}
### now we have the proportion of each land cover type at each site
### we can calculate the Shannon Diversity at each site (how much diversity of habitats there is at each site)
### sites with higher Shannon diversity are more heterogeneous

### calculate heterogeneity using Shannon diversity index
### for each site under 4 buffers: 500m, 2km, 5km, and 10km 
bbs_hab$shannon_div <- diversity(bbs_hab[,5:17], index="shannon")


