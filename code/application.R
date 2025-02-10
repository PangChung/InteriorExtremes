rm(list=ls())
args <- commandArgs(TRUE)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
computer = "local"
id = 1
idx.jack = 1
method = "Nelder-Mead"

for (arg in args) eval(parse(text = arg))
switch(computer,
    "ws" = {DataPath<-"~/Desktop/InteriorExtremes/"},
    "hpc" = {DataPath<-"/srv/scratch/z3536974/";.libPaths("../src")},
    "local" = {DataPath<-"~/Documents/Github/InteriorExtremes/"}
)
idx.jack=idx.jack-1
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(partitions)
library(Rfast)
library(matrixStats)
init.seed = as.integer((as.integer(Sys.time()) + sample.int(10^5,1))%%10^5)
set.seed(init.seed)
## load the data##
load("data/data_application.RData")
load("data/Trends_fits.RData")

#     data <- residuals
#     data[data<0] = 0
#     data <- apply(data,2,function(x) {ind = x>0;x[ind] = qgpd(rank(x[ind])/(sum(ind)+1),1,1,1); x })
#     data = apply(data,1,function(x){list(which(x>0),x[x>0])})
#     idx.data = which(sapply(data,function(x) length(x[[2]]))>0)
#     data = data[idx.data]
#     data.sum = sapply(data,function(x) mean(x[[2]]))
#     data.max = sapply(data,function(x) max(x[[2]]))

#     # thres <- quantile(data.sum, seq(0.9,0.9999,length.out=30))
#     # mev::tstab.gpd(data.sum,thres)

#     # thres <- quantile(data.max, seq(0.9,0.9999,length.out=30))
#     # mev::tstab.gpd(data.max,thres)

#     data.fit.sum = data[data.sum>quantile(data.sum,0.995)]
#     data.fit.max = data[data.max>quantile(data.max,0.995)]

data <- maxima.frechet

if(idx.jack != 0){
    boot.ind <- sample(1:nrow(data),nrow(data),replace = TRUE)
    data <- data[boot.ind,]
}
data.sum = apply(data,1,sum)
ind.data = which(data.sum > 10*1043)
length(ind.data)

data.fit = data[ind.data,]
data.fit = apply(data.fit,1,function(x){list(which(x>0),x[x>0])})

D = nrow(loc.sub.trans)
ncores = detectCores()
idx.centers = unlist(lapply(quantile(loc.sub.trans[,1],seq(0.1,0.9,length.out=5)),function(x){ idx = abs(loc.sub.trans[,1] - x) < 5; which(idx)[which.min(abs(loc.sub.trans[idx,2] - median(loc.sub.trans[idx,2])))]}))

basis <- sapply(idx.centers,function(x){y=dnorm(distmat[x,],mean=0,sd=ncol(distmat)*2);y=y-mean(y);y/sqrt(sum(y^2))})

n.alpha = ncol(basis)
init = c(100,1,0,1,rep(0,n.alpha))
idx.para = 1:4

ub = c(Inf,1.99,pi/4,Inf,rep(Inf,n.alpha))
lb = c(0.01,0.01,-pi/4,0.01,rep(-Inf,n.alpha))

file.save = paste0(DataPath,"/data/application_RedSea_results_",id,"_",method,"_",idx.jack,".RData")
file.origin = paste0(DataPath,"/data/application_RedSea_results_",id,"_",method,".RData")

if(file.exists(file.save)){
    init = e$fit.result$par
    # stop("job already done")
}

if(file.exists(file.origin)){
    load(file.origin,e<-new.env())
    init = e$fit.result$par
}

all.index <- combn(x=D,m=2)
dist.max = distmat[t(all.index)]
dist.ind<-which(rank(dist.max,ties.method = "first") <= D*2)
all.index = all.index[,dist.ind]

vecchia.seq <- order(distmat[which.min(apply(distmat,2,mean))[1],]) ## Vecchia sequence based on middle-out ordering
neighbours.mat <- sapply(1:D,FUN=neighbours,vecchia.seq=vecchia.seq,q=4,loc=distmat)

switch(id,
        {fit.result <- fit.model(data=data.fit,loc=loc.sub.trans,init=init,fixed=c(F,F,F,F,rep(T,ncol(basis))),model="logskew",maxit=1e6,FUN=vario.func2,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)},
        {fit.result <- fit.model(data=data.fit,loc=loc.sub.trans,init=init,fixed=c(F,F,F,F,rep(F,ncol(basis))),model="logskew",maxit=1e6,FUN=vario.func2,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)},
        {fit.result <- MCLE(data=data,init=init[1:4],fixed=c(F,F,F,F),loc=loc.sub.trans,FUN=vario.func2,index=all.index,model="BR",lb=lb[1:4],ub=ub[1:4],ncores=ncores,maxit=1e6,trace=TRUE,idx.para=idx.para)},
        {fit.result <- MVLE(data=data,init=init[1:4],fixed=c(F,F,F,F),loc=loc.sub.trans,FUN=vario.func2,vecchia.seq=vecchia.seq,neighbours=neighbours.mat,model="BR",lb=lb[1:4],ub=ub[1:4],ncores=ncores,maxit=1e6,trace=TRUE,idx.para=idx.para)}
    )
    
if(idx.jack!=0){
    save(fit.result,idx.centers,basis,file=file.save)
}else{
    save(fit.result,idx.centers,basis,file=file.origin)
}
