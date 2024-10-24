rm(list=ls())
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(partitions)
library(Rfast)
library(evd)
library(matrixStats)
library(Matrix)
library(sf)
library(ggplot2)
library(ggmap)
# load the data##
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
load("data/application_florida_list.RData")
grid.sf <- read_sf("data/Florida/DOPGrid/DOPGrid.shp")
intb.sf <- read_sf("data/Florida/INTB_Basins/INTB_Basins.shp")
fill_basin <- as.factor(rep(1:8,length.out=172))
all_IDs <- names(read.csv("data/Florida/PixelRain15min_1995.csv", header = TRUE, nrows = 1))[-1]
all_IDs_num <- as.numeric(stringi::stri_extract_first(all_IDs, regex = "[0-9]+"))
fill_grid = grid.sf$PIXEL %in% all_IDs_num

# p <- ggplot() + geom_sf(data=intb.sf, aes(fill=fill_basin)) + 
#     scale_fill_brewer(palette = "RdBu",name="Basins")  + 
#     geom_sf(data=grid.sf, aes(colour=as.factor(fill_grid)),alpha=0.1) + 
#     scale_color_brewer(palette = "Set1",name="Grid") + xlim(c(317734,433386)) + ylim(3047561,3191794)
# p

register_google(key="AIzaSyA4GHwem8rdsXo1mgSNx-8uUqEDOjr39Ds")
map <- get_googlemap(center=c(-82.273528,28.209394),zoom=9,maptype = "terrain",style = "feature:all|element:all|saturation:-100|lightness:50")
intb.transform <- st_transform(intb.sf, crs = 4326)
grid.transform <- st_transform(grid.sf, crs = 4326)

save(intb.transform,grid.transform,map,grid.sf,intb.sf,file="data/maps_florida.RData")

ggmap(map) + theme_void() + 
ggtitle("Tampa Bay") + theme(plot.title = element_text(hjust = 0.5))  + 
geom_sf(data=grid.transform, aes(colour=as.factor(fill_grid)),alpha=0.1,inherit.aes = FALSE) + scale_color_brewer(palette = "Set1",name="Grid")

