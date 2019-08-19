###################################################################################
## Title: Map of Chapter 2 sites 
## User: Lisbeth Morrison 
## email: l.morrison@pgr.reading.ac.uk
## Date: May 2019
###################################################################################

rm(list=ls()) # clear R


library(ggmap)
library(maps)
library(ggplot2)
library(gganimate)
library(tidyverse)
library(sf)
library(gifski)
library(mapproj)

## add data
BBS_sites <- read.csv("../Data/BBS_site_data/BBS_filtered_sites_final.csv", header=TRUE)

sites <- unique(subset(BBS_sites[c(1,6,7)]))

lat_long <- sites %>%
  st_as_sf(coords = c("EASTING", "NORTHING"), crs = 27700) %>%
  st_transform(4326) %>%
  st_coordinates() %>%
  as_tibble()
site_list <- cbind(sites, lat_long)
site_list <- transform( site_list, gridref = sample(gridref))

UK <- map_data("world") %>% filter(region=="UK")

p <- 
  ggplot() +
  geom_polygon(data = UK, aes(x=long, y = lat, group = group), fill="grey35", colour="black", alpha=0.3) +
  geom_point(data=site_list, aes(x=X, y=Y), size = 2, shape = 21, fill="black", alpha = 0.7) +
  theme_void() + coord_map() +
  
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(hjust = 0, size = 16, vjust=0))
p  
ggsave(filename="../Graphs/map_of_sites.png", width=5, height=7)
