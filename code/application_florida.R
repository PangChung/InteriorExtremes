rm(list=ls())
args <- commandArgs(TRUE)
computer = "local"
id = 1
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
library(evd)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
set.seed(12342)

## load the data##
# library(matrixStats)
# library(Matrix)
# library(sf)
# library(ggplot2)
# grid.sf <- read_sf("data/Florida/DOPGrid/DOPGrid.shp")
# intb.sf <- read_sf("data/Florida/INTB_Basins/INTB_Basins.shp")

# fill_basin <- as.factor(rep(1:8,length.out=172))
# all_IDs <- names(read.csv("data/Florida/PixelRain15min_1995.csv", header = TRUE, nrows = 1))[-1]
# all_IDs_num <- as.numeric(stringi::stri_extract_first(all_IDs, regex = "[0-9]+"))
# fill_grid = grid.sf$PIXEL %in% all_IDs_num

# p <- ggplot() + geom_sf(data=intb.sf, aes(fill=fill_basin)) + scale_fill_brewer(palette = "RdBu",name="Basins")  + geom_sf(data=grid.sf, aes(colour=as.factor(fill_grid)),alpha=0.1) + scale_color_brewer(palette = "Set1",name="Grid") + xlim(c(317734,433386)) + ylim(3047561,3191794)

## extract the data for the year 1995--2019 ##
# format_string <- "%d-%b-%Y %H:%M:%S"
# column_classes <- rep("NULL", length(all_IDs));column_classes[1] <- "character"
# for(y in 1995:2019){
#     column_classes[i] <- "numeric"
#     pipe_command <- paste0("cut -f1 -d, ","data/Florida/PixelRain15min_",y,".csv")
#     date.i <- scan(pipe(pipe_command),what="character")[-1]
#     date.i <- paste(date.i[1:length(date.i) %% 2 == 1], date.i[1:length(date.i) %% 2 == 0], sep=" ")
#     date.i <- as.POSIXct(strptime(date.i, format_string),tz="UTC")
#     ind.date <- date.i >= as.POSIXct(paste0(y,"-06-20 00:00:00"),tz="UTC") & date.i <= as.POSIXct(paste0(y,"-09-22 23:45:00"),tz="UTC")
#     rm(date.i,data)
#     data = read.csv(paste0("data/Florida/PixelRain15min_",y,".csv"), header = TRUE,sep=",")[ind.date,]
#     if(y==1995){
#         write.table(data, file="data/Florida/PixelRain15min_1995_2019.csv",sep=",",row.names=FALSE,col.names=TRUE)
#     }else{
#         write.table(data, file="data/Florida/PixelRain15min_1995_2019.csv",sep=",",row.names=FALSE,col.names=FALSE,append=TRUE)
#     }
#     print(y)
# }

# date <- scan(pipe("cut -f1 -d, data/Florida/PixelRain15min_1995_2019.csv"),what="character")[-1]
# data <- scan(pipe("cut -f2 -d, data/Florida/PixelRain15min_1995_2019.csv"),what="numeric")[-1]
# data <- as.numeric(data)

# tmp <- read.csv("data/Florida/PixelRain15min_1995_2019.csv", header = TRUE,sep=",")
# save(tmp,file="data/application_florida.RData")

# load("data/application_florida.RData")
# num.nonzeros=rep(NA,length(fill_grid))
# num.nonzeros[fill_grid] <- apply(tmp[,-1],2,function(x) sum(x>0))

# p <- ggplot() + geom_sf(data=intb.sf,color="black") + geom_sf(data=grid.sf, color="grey50",aes(fill=num.nonzeros),alpha=0.75) +
# scale_fill_distiller(type="seq",name="Intensity",na.value="grey50",trans="reverse") + xlim(c(317734,433386)) + ylim(3047561,3191794) 
# p

# func <- function(i){
#     x = tmp[,i]
#     ind.x = x>0
#     x[ind.x] = rank(x[ind.x])/(sum(ind.x)+1)
#     return(x)
# }

# data.uniform <- matrix(unlist(mclapply(2:ncol(tmp),func,mc.cores=4)),nrow=nrow(tmp),byrow=FALSE)
# rm(tmp);gc()

# func <- function(i){
#     x = data.uniform[i,]
#     ind.x = which(x>0)
#     return(list(ind.x,x[ind.x]))
# }

# nonzeros_row <- apply(data.uniform,1,function(x) any(x>0))

# num.nonzeros_row <- apply(data.uniform,1,function(x){sum(x>0)})

# system.time({data <- mclapply(which(nonzeros_row),func,mc.cores=4)})

# save(data, num.nonzeros_row, nonzeros_row, file="data/application_florida_list.RData")

# data.pareto <- mclapply(data,function(x){list(x[[1]],qgpd(x[[2]],1,1,1))},mc.cores=4)

# data.sum <- unlist(mclapply(data.pareto,function(x){mean(x[[2]])},mc.cores=4))

# data.max <- unlist(mclapply(data.pareto,function(x){max(x[[2]])},mc.cores=4))

# thres <- quantile(data.sum, seq(0.999,0.9999,length.out=30))
# tstab.gpd(data.sum,thresh=thres,plot=TRUE)

# thres <- quantile(data.max, seq(0.999,0.9999,length.out=30))
# tstab.gpd(data.max,thresh=thres,plot=TRUE)

# data.fit.sum <- data.pareto[data.sum>quantile(data.sum,0.9995)]
# data.fit.max <- data.pareto[data.max>quantile(data.max,0.9995)]

# coord <- as.data.frame(st_coordinates(grid.sf))
# idx.pixel <- unlist(lapply(all_IDs_num,function(x){which(grid.sf$PIXEL==x)}))
# coord.grid <- matrix(unlist(lapply(idx.pixel, function(i){apply(coord[coord$L2==i,1:2][1:4,],2,mean)})),ncol=2,byrow=TRUE)

# ggplot() + geom_point(data=coord,aes(X,Y),color="black") + 
# geom_point(data=subset(coord,L2 %in% idx.pixel),aes(x=X,y=Y),color="red") + geom_point(data=as.data.frame(coord.grid),aes(V1,V2),color="yellow")


# system.time({cov.mat <- vario.func2(coord.grid,c(60000,1,1,1))})
# a = as.matrix(dist(coord.grid))

### fit the model ###
load("data/application_florida_list.RData")
coord.grid = coord.grid/60000
### fit the skewd-BR model ###
basis.centers <- expand.grid(quantile(coord.grid[,1],c(0.2,0.5,0.8)),quantile(coord.grid[,2],c(0.2,0.5,0.8)))
basis <- lapply(1:nrow(basis.centers),function(i){
    y=dnorm(sqrt((coord.grid[,1]-basis.centers[i,1])^2 + (coord.grid[,2]-basis.centers[i,2])^2),mean=0,sd=ncol(coord.grid)*2)
    y=y-mean(y)
    y/sqrt(sum(y^2))
})
unit = 2/0.03128403
basis <- matrix(unlist(basis),nrow=nrow(coord.grid),byrow=FALSE)
idx.para = c(1:4)
ncores = detectCores()/2
init = c(0.1,1,0,1,rep(0,nrow(basis.centers)))
ub = c(Inf,1.99,pi/4,Inf,rep(Inf,nrow(basis.centers)))
lb = c(0.01,0.01,-pi/4,0.01,rep(-Inf,nrow(basis.centers)))

switch(id,
    {fit.result <- fit.model(data=data.fit.sum,loc=coord.grid,init=init,fixed=c(F,F,F,F,rep(T,nrow(basis.centers))),model="logskew",maxit=1000,FUN=vario.func2,basis=basis,alpha.func=alpha.func,ncores=ncores,method="Nelder-Mead",lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)},
    {fit.result <- fit.model(data=data.fit.max,loc=coord.grid,init=init,fixed=c(F,F,F,F,rep(T,nrow(basis.centers))),model="logskew",maxit=1000,FUN=vario.func2,basis=basis,alpha.func=alpha.func,ncores=ncores,method="Nelder-Mead",lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)},
    {fit.result <- fit.model(data=data.fit.sum,loc=coord.grid,init=init,fixed=c(F,F,F,F,rep(F,nrow(basis.centers))),model="logskew",maxit=1000,FUN=vario.func2,basis=basis,alpha.func=alpha.func,ncores=ncores,method="Nelder-Mead",lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)},
    {fit.result <- fit.model(data=data.fit.max,loc=coord.grid,init=init,fixed=c(F,F,F,F,rep(F,nrow(basis.centers))),model="logskew",maxit=1000,FUN=vario.func2,basis=basis,alpha.func=alpha.func,ncores=ncores,method="Nelder-Mead",lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)})

save(fit.result,basis.centers,file=paste0(DataPath,"/data/application_florida_results_",id,".RData"))