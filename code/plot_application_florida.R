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
load("data/application_florida/application_florida_results_3_L-BFGS-B.RData")
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

fitted.extcoef <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x];V_bi_logskew(c(1,1),list(par.list[[1]][x,x],par.list[[2]][x]),alpha.para=FALSE)},mc.cores=5,mc.set.seed = FALSE))

range(fitted.extcoef)

fitted.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=fitted.extcoef,symmetric = TRUE,dimnames=NULL) 


coord.ratio <- 0.01796407/0.01912047
basis.centers.geo <- expand.grid(quantile(coord.geo[,1],seq(0.1,0.9,length.out=10)),quantile(coord.geo[,2],seq(0.1,0.9,length.out=10)))
p <- list()
brks = round(quantile(fitted.extcoef.mat@x,probs=c(0.001,0.005,0.01,0.05,0.1,0.2,0.5,0.8),na.rm=TRUE),4)
for(i in 1:nrow(basis.centers.geo)){
    center.coord <- basis.centers.geo[i,]
    idx.center = which.min(apply(coord.geo,1,function(x){sum(x-center.coord)^2}))
    data.df <- data.frame(lon=round(coord.geo[,1],5),lat=round(coord.geo[,2],5),z=fitted.extcoef.mat[,idx.center])
    data.df$z[idx.center] = NA
    p[[i]]<-ggmap(map) + theme_void() + 
    ggtitle("Tampa Bay") + theme(plot.title = element_text(hjust = 0.5)) + geom_tile(data=data.df,aes(x=lon,y=lat,fill=z),alpha=0.5) + scale_fill_distiller(name="Extremal Coefficient",palette = "Spectral",trans="reverse") + coord_fixed(ratio=1/coord.ratio) + stat_contour(data=data.df,aes(x=lon,y=lat,z=z),breaks = brks,colour = "black")
}

for(i in 1:length(p)){
    png(paste0("figures/application/florida_extcoef_",sprintf(i,fmt="%.3d"),".png"),width=800,height=800)
    print(p[[i]])
    dev.off()
}

load("data/application_florida/application_florida_results_ext_1.RData",e<-new.env())
load("data/application_florida/application_florida_fitted_extcoef.RData",e1<-new.env())
load("data/application_florida/application_florida_fitted_extcoef_3.RData",e2<-new.env())
emp.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=e$emp.extcoef,symmetric = TRUE,dimnames=NULL)
p <- list()
brks = round(quantile(c(e1$fitted.extcoef.mat@x,e2$fitted.extcoef.mat@x),probs=c(0.001,0.005,0.01,0.05,0.1,0.2,0.5,0.8),na.rm=TRUE),4)
for(i in 1:nrow(basis.centers.geo)){
    center.coord <- basis.centers.geo[i,]
    idx.center = which.min(apply(coord.geo,1,function(x){sum(x-center.coord)^2}))
    data.df <- data.frame(lon=round(coord.geo[,1],5),lat=round(coord.geo[,2],5),
                emp=2-emp.extcoef.mat[,idx.center],br=e1$fitted.extcoef.mat[,idx.center],
                sbr=e2$fitted.extcoef.mat[,idx.center])
    data.df[idx.center,-c(1:2)] = NA
    p[[i]]<-ggmap(map) + theme_void() + 
    ggtitle("Tampa Bay") + theme(plot.title = element_text(hjust = 0.5)) + geom_tile(data=data.df,aes(x=lon,y=lat,fill=emp),alpha=0.5) + scale_fill_distiller(name="Extremal Coefficient",palette = "Spectral",trans="reverse") + coord_fixed(ratio=1/coord.ratio) + stat_contour(data=data.df,aes(x=lon,y=lat,z=br),breaks = brks,colour = "black",linetype="dashed") + stat_contour(data=data.df,aes(x=lon,y=lat,z=sbr),breaks = brks,colour = "black")
}

for(i in 1:length(p)){
    png(paste0("figures/application/florida_extcoef_",sprintf(i,fmt="%.3d"),".png"),width=800,height=800)
    print(p[[i]])
    dev.off()
}

save(fitted.extcoef.mat, par.list, file="data/application_florida/application_florida_fitted_extcoef_3.RData")

# data.pareto <- mclapply(data,function(x){list(x[[1]],qgpd(x[[2]],1,1,1))},mc.cores=4)
# len.row <- unlist(lapply(1:length(data.pareto),function(i){length(data.pareto[[i]][[1]])}))
# data.pareto.mat <- sparseMatrix(i=rep(1:length(data.pareto),times=len.row),j=unlist(lapply(1:length(data.pareto),function(i){data.pareto[[i]][[1]]})),x=unlist(lapply(1:length(data.pareto),function(i){data.pareto[[i]][[2]]})),dimnames=NULL,symmetric = FALSE)
# rm(data.pareto,data);gc()
# empirical.extcoef <- function(data){
#     u=10
#     x = data[,1]
#     y= data[,2]
#     return( sum(x>u & y>u)/(sum(x>u)+sum(y>u))*2)
# }
# data.pareto.mat.nonsparse <- as.matrix(data.pareto.mat)
# system.time({emp.extcoef <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x]; empirical.extcoef(data.pareto.mat.nonsparse[,x])},mc.cores=5,mc.set.seed = FALSE))})

# save(data.pareto.mat,emp.extcoef,fitted.extcoef.mat,file="data/application_florida/application_florida_results_ext_1.RData")

