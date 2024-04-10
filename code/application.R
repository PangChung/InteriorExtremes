rm(list=ls())
args <- commandArgs(TRUE)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
computer = "ws";id=2
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
idx.centers = unlist(lapply(quantile(loc.sub.trans[,1],seq(0.15,0.85,length.out=3)),function(x){ idx = abs(loc.sub.trans[,1] - x) < 5; which(idx)[which.min(abs(loc.sub.trans[idx,2] - median(loc.sub.trans[idx,2])))]}))
#idx.centers = floor(matrix(seq(1,nrow(distmat),length.out=6),ncol=2,3))

#plot(x=loc.sub.trans[,1],y=loc.sub.trans[,2],xlab="x",ylab="y",main="Data",pch=20)
#points(x=loc.sub.trans[idx.centers,1],y=loc.sub.trans[idx.centers,2],col="red",pch=20)

loc.sub.trans = apply(loc.sub.trans,2,function(x) x-mean(x))
basis <- matrix(0,nrow=D,ncol=length(idx.centers)+1)
basis[1:floor(D/2),1] = 0.1;basis[(D-floor(D/2)+1):D,1] = -0.1
basis[,-1] <- sapply(idx.centers,function(x){y=-distmat[,x]/max(distmat[,x]);y=y-mean(y)})


# basis[1:floor(D/2),1] = 0.1;basis[(D-floor(D/2)+1):D,1] = -0.1
# basis[,-1] <- apply(idx.centers,1,function(x){y <- rep(0,nrow(distmat));y[x] <- c(-2,2);y})

init = c(109,1.16,rep(0,length(idx.centers)))
n.alpha = ncol(basis)-1
switch(id,
    results1 <- fit.model(data=maxima.frechet,loc=distmat,init=init,fixed=c(F,F,rep(F,n.alpha)),basis=basis,thres=2,model="logskew",maxit=1000,FUN=cov.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,n.alpha)),ub=c(Inf,1.99,rep(Inf,n.alpha)),hessian=FALSE,opt=TRUE,trace=TRUE), 
    
    {results4 <- fit.model(data=maxima.frechet,loc=loc.sub.trans,init=init,fixed=c(F,F,rep(T,n.alpha)),basis=matrix(0,ncol=n.alpha+1,nrow=D),thres=15,model="logskew",maxit=1000,FUN=vario.func,alpha.func=alpha.func,ncores=ncores,method="Nelder-Mead",lb=c(0.01,0.0,rep(-Inf,n.alpha)),ub=c(Inf,1.99,rep(Inf,n.alpha)),hessian=FALSE,opt=TRUE,trace=TRUE,idx.para=1:2)

    results2 <- fit.model(data=maxima.frechet,loc=loc.sub.trans,init=results4$par,fixed=c(T,T,rep(F,n.alpha)),basis=basis,thres=15,model="logskew",maxit=1000,FUN=vario.func,alpha.func=alpha.func,ncores=ncores,method="Nelder-Mead",lb=c(0.01,0.0,rep(-Inf,n.alpha)),ub=c(Inf,1.99,rep(Inf,n.alpha)),hessian=FALSE,opt=TRUE,trace=TRUE,step2=FALSE,idx.para=1:2)},

    results3 <- fit.model(data=maxima.frechet,loc=distmat,init=init,fixed=c(F,F,rep(T,n.alpha)),basis=basis,thres=0.9,model="logskew",maxit=1000,FUN=cov.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,n.alpha)),ub=c(Inf,1.99,rep(Inf,n.alpha)),hessian=FALSE,opt=TRUE,trace=TRUE)
)

save.image(file=paste0(DataPath,"data/application_results_new",id,".RData"))
