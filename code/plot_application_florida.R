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
library(Matrix)
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

load("data/maps_florida.RData")
register_google(key=system("echo $g_key",intern = TRUE))
map <- get_googlemap(center=c(-82.273528,28.209394),zoom=9,maptype = "terrain",style = "feature:all|element:all|saturation:-100|lightness:50")
intb.transform <- st_transform(intb.sf, crs = 4326)
grid.transform <- st_transform(grid.sf, crs = 4326)

save(intb.transform,grid.transform,map,grid.sf,intb.sf,file="data/maps_florida.RData")

ggmap(map) + theme_void() + 
ggtitle("Tampa Bay") + theme(plot.title = element_text(hjust = 0.5))  + 
geom_sf(data=grid.transform, aes(colour=as.factor(fill_grid)),alpha=0.5,inherit.aes = FALSE) + scale_color_brewer(palette = "Set1",name="Grid")

coord.geo <- as.data.frame(st_coordinates(grid.transform))
idx.pixel <- unlist(lapply(all_IDs_num,function(x){which(grid.transform$PIXEL==x)}))
coord.geo <- matrix(unlist(lapply(idx.pixel, function(i){apply(coord.geo[coord.geo$L2==i,1:2][1:4,],2,mean)})),ncol=2,byrow=TRUE)
coord.grid = coord.grid/60000
coord.grid.new <- apply(coord.grid,2,function(x){x-median(x)})
load("data/application_florida/application_florida_results_1_Nelder-Mead.RData")
basis <- lapply(1:nrow(basis.centers),function(i){
    y=dnorm(sqrt((coord.grid[,1]-basis.centers[i,1])^2 + (coord.grid[,2]-basis.centers[i,2])^2),mean=0,sd=ncol(coord.grid)*2)
    y=y-mean(y)
    y/sqrt(sum(y^2))
})
basis <- matrix(unlist(basis),nrow=nrow(coord.grid),byrow=FALSE)
alpha <- alpha.func(fit.result$par[-c(1:4)],basis)
cov.mat <- vario.func2(coord.grid.new,fit.result$par[1:4],ncores=5)
par.list <- alpha2delta(list(cov.mat,alpha))
pairs <- comb_n(1:nrow(coord.grid),2)
idx = pairs[,1]
V_bi_logskew(c(1,1),list(par.list[[1]][idx,idx],par.list[[2]][idx]),alpha.para=FALSE)

fitted.extcoef <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x];V_bi_logskew(c(1,1),list(par.list[[1]][x,x],par.list[[2]][x]),alpha.para=FALSE)},mc.cores=5,mc.set.seed = FALSE))

range(fitted.extcoef)

fitted.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=fitted.extcoef,symmetric = TRUE,dimnames=NULL) 