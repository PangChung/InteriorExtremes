rm(list=ls())
args <- commandArgs(TRUE)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
computer = "local"
for (arg in args) eval(parse(text = arg))
switch(computer,
    "ws" = {DataPath<-"~/Desktop/InteriorExtremes/"},
    "hpc" = {DataPath<-"/srv/scratch/z3536974/";.libPaths("../src")},
    "local" = {DataPath<-"~/Documents/Github/InteriorExtremes/"}
)
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(partitions)
library(Rfast)
library(matrixStats)
library(sf)
library(ggplot2)
set.seed(12342)
## load the data##

data.basin <- read.csv2("data/Florida/BasinRain15min.csv", header = TRUE, sep = ",")
grid.sf <- read_sf("data/Florida/DOPGrid/DOPGrid.shp")
intb.sf <- read_sf("data/Florida/INTB_Basins/INTB_Basins.shp")

fill_factor <- as.factor(rep(1:8,length.out=172))

ggplot() + geom_sf(data=grid.sf,fill="red",color="black") + coord_sf()

data <- read.csv2("data/Florida/PixelRain15min_1995.csv", header = TRUE, sep = ",")
IDs <- as.numeric(stringi::stri_extract_first(names(data), regex = "[0-9]+"))

fill_grid = grid.sf$PIXEL %in% IDs

ggplot() + geom_sf(data=intb.sf, aes(fill=fill_factor)) + scale_fill_brewer(palette = "RdBu",name="Basins")  + geom_sf(data=grid.sf, aes(colour=as.factor(fill_grid)),alpha=0.1) + scale_color_brewer(palette = "RdBu",name="Grid")

