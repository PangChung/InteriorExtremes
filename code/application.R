rm(list=ls())
args <- commandArgs(TRUE)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
computer = "hpc";id=1
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
set.seed(12342)
## load the data##
load("data/data_application.RData")
D = ncol(maxima.frechet)
ncores = floor(detectCores()/2)
idx.centers = unlist(lapply(quantile(loc.sub.trans[,1],c(0.2,0.5,0.8)),function(x){ idx = abs(loc.sub.trans[,1] - x) < 5; which(idx)[which.min(abs(loc.sub.trans[idx,2] - median(loc.sub.trans[idx,2])))]}))

#plot(x=loc.sub.trans[,1],y=loc.sub.trans[,2],xlab="x",ylab="y",main="Data",pch=20)
#points(x=loc.sub.trans[idx.centers,1],y=loc.sub.trans[idx.centers,2],col="red",pch=20)

basis <- matrix(0,nrow=D,ncol=4)
basis[1:floor(D/2),1] = 0.1;basis[(D-floor(D/2)+1):D,1] = -0.1
basis[,-1] <- sapply(idx.centers,function(x){y=dnorm(distmat[,x],mean=0,sd=500);y=y-mean(y);y/sqrt(sum(y^2))})

switch(id,
    results1 <- fit.model(data=maxima.frechet,loc=distmat,init=c(60,1.2,0,0,0),fixed=c(F,F,F,F,F),basis=basis,thres=0.9,model="logskew",maxit=1000,FUN=cov.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,3)),ub=c(Inf,1.99,rep(Inf,3)),hessian=FALSE,opt=TRUE,trace=FALSE), 
    results2 <- fit.model(data=maxima.frechet,loc=loc.sub.trans,init=c(100,1,0,0,0),fixed=c(F,F,F,F,F),basis=basis,thres=0.9,model="logskew",maxit=1000,FUN=vario.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,3)),ub=c(Inf,1.99,rep(Inf,3)),hessian=FALSE,opt=TRUE,trace=FALSE),
    results3 <- fit.model(data=maxima.frechet,loc=distmat,init=c(60,1.2,0,0,0),fixed=c(F,F,T,T,T),basis=basis,thres=0.9,model="logskew",maxit=1000,FUN=cov.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,3)),ub=c(Inf,1.99,rep(Inf,3)),hessian=FALSE,opt=TRUE,trace=FALSE),
    results4 <- fit.model(data=maxima.frechet,loc=loc.sub.trans,init=c(60,1.2,0,0,0),fixed=c(F,F,T,T,T),basis=matrix(0,ncol=4,nrow=D),thres=0.9,model="logskew",maxit=1000,FUN=vario.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,3)),ub=c(Inf,1.99,rep(Inf,3)),hessian=FALSE,opt=TRUE,trace=FALSE)
)

save.image(file=paste0(DataPath,"data/application_results",id,".RData"))
