rm(list=ls())
args <- commandArgs(TRUE)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
computer = "local";id=2
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
idx.centers = unlist(lapply(quantile(loc.sub.trans[,1],seq(0.1,0.9,length.out=5)),function(x){ idx = abs(loc.sub.trans[,1] - x) < 5; which(idx)[which.min(abs(loc.sub.trans[idx,2] - median(loc.sub.trans[idx,2])))]}))
#idx.centers = floor(matrix(seq(1,nrow(distmat),length.out=6),ncol=2,3))

#plot(x=loc.sub.trans[,1],y=loc.sub.trans[,2],xlab="x",ylab="y",main="Data",pch=20)
#points(x=loc.sub.trans[idx.centers,1],y=loc.sub.trans[idx.centers,2],col="red",pch=20)

loc.sub.trans = apply(loc.sub.trans,2,function(x) x-min(x)+1)

# basis <- matrix(0,nrow=D,ncol=length(idx.centers))
# basis[1:floor(D/2),1] = 0.1;basis[(D-floor(D/2)+1):D,1] = -0.1
basis <- sapply(idx.centers,function(x){y=-distmat[,x]/max(distmat[,x]);y=y-mean(y)})


# basis[1:floor(D/2),1] = 0.1;basis[(D-floor(D/2)+1):D,1] = -0.1
# basis[,-1] <- apply(idx.centers,1,function(x){y <- rep(0,nrow(distmat));y[x] <- c(-2,2);y})

init = c(100,1,rep(0,ncol(basis)))
n.alpha = ncol(basis)
idx.para = 1:2
switch(id,
    results1 <- fit.model(data=maxima.frechet,loc=distmat,init=init,fixed=c(F,F,rep(F,n.alpha)),basis=basis,thres=2,model="logskew",maxit=1000,FUN=cov.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,n.alpha)),ub=c(Inf,1.99,rep(Inf,n.alpha)),hessian=FALSE,opt=TRUE,trace=TRUE,idx.para=idx.para), 
    
    {results4 <- fit.model(data=maxima.frechet,loc=loc.sub.trans,init=init[idx.para],fixed=c(F,F),thres=14,model="BR",maxit=1000,FUN=vario.func,ncores=ncores,method="Nelder-Mead",lb=c(0.01,0.0),ub=c(Inf,1.99),hessian=FALSE,opt=TRUE,trace=TRUE,idx.para=1:2)
    init[1:2] = results4$par
    results2 <- fit.model(data=maxima.frechet,loc=loc.sub.trans,init=init,fixed=c(F,F,rep(F,n.alpha)),basis=basis,thres=14,model="logskew",maxit=1000,FUN=vario.func,alpha.func=alpha.func,ncores=ncores,method="Nelder-Mead",lb=c(0.01,0.0,rep(-Inf,n.alpha)),ub=c(Inf,1.99,rep(Inf,n.alpha)),hessian=FALSE,opt=TRUE,trace=TRUE,step2=TRUE,idx.para=1:2)},

    results3 <- fit.model(data=maxima.frechet,loc=distmat,init=init,fixed=c(F,F,rep(T,n.alpha)),basis=basis,thres=0.9,model="logskew",maxit=1000,FUN=cov.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,n.alpha)),ub=c(Inf,1.99,rep(Inf,n.alpha)),hessian=FALSE,opt=TRUE,trace=TRUE,idx.para=1:2)
)

save.image(file=paste0(DataPath,"data/application_results_new_",id,".RData"))
