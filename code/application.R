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
set.seed(12342)
## load the data##
load("data/data_application.RData")
D = ncol(maxima.frechet)
ncores = floor(detectCores()/2)
idx.centers = unlist(lapply(quantile(loc.sub.trans[,1],seq(0.1,0.9,length.out=5)),function(x){ idx = abs(loc.sub.trans[,1] - x) < 5; which(idx)[which.min(abs(loc.sub.trans[idx,2] - median(loc.sub.trans[idx,2])))]}))
#idx.centers = floor(matrix(seq(1,nrow(distmat),length.out=6),ncol=2,3))

#plot(x=loc.sub.trans[,1],y=loc.sub.trans[,2],xlab="x",ylab="y",main="Data",pch=20)
#points(x=loc.sub.trans[idx.centers,1],y=loc.sub.trans[idx.centers,2],col="red",pch=20)

#loc.sub.trans = apply(loc.sub.trans,2,function(x) x-min(x)+1)

# basis <- matrix(0,nrow=D,ncol=length(idx.centers))
# basis[1:floor(D/2),1] = 0.1;basis[(D-floor(D/2)+1):D,1] = -0.1
# basis <- sapply(idx.centers,function(x){y=-distmat[,x]/max(distmat[,x]);y=y-mean(y)})
basis <- sapply(idx.centers,function(x){y=dnorm(distmat[x,],mean=0,sd=ncol(distmat)*2);y=y-mean(y);y/sqrt(sum(y^2))})

# basis[1:floor(D/2),1] = 0.1;basis[(D-floor(D/2)+1):D,1] = -0.1
# basis[,-1] <- apply(idx.centers,1,function(x){y <- rep(0,nrow(distmat));y[x] <- c(-2,2);y})

init = c(100,1,rep(0,ncol(basis)))
n.alpha = ncol(basis)
idx.para = 1:2

thres = 10

results2 <- fit.model(data=maxima.frechet,loc=loc.sub.trans,init=init,fixed=c(F,F,rep(F,n.alpha)),basis=basis,thres=thres,model="logskew",maxit=1000,FUN=vario.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,n.alpha)),ub=c(Inf,1.99,rep(Inf,n.alpha)),hessian=FALSE,opt=TRUE,trace=TRUE,step2=TRUE,idx.para=1:2)

results21 <- fit.model(data=maxima.frechet,loc=loc.sub.trans,init=results2$par,fixed=c(F,F,rep(F,n.alpha)),basis=basis,thres=thres,model="logskew",maxit=1000,FUN=vario.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,n.alpha)),ub=c(Inf,1.99,rep(Inf,n.alpha)),hessian=TRUE,opt=FALSE,trace=TRUE,step2=FALSE,idx.para=1:2)

U = eigen(results21$K)$values
V = eigen(results21$K)$vectors
K.inv = V %*% diag(1/U) %*% t(V)
sqrt(diag(K.inv %*% results21$hessian %*% K.inv))

results4 <- fit.model(data=maxima.frechet,loc=loc.sub.trans,init=init[idx.para],fixed=c(F,F),thres=thres,model="BR",maxit=1000,FUN=vario.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0),ub=c(Inf,1.99),hessian=FALSE,opt=TRUE,trace=TRUE,idx.para=1:2)

results41 <- fit.model(data=maxima.frechet,loc=loc.sub.trans,init=results4$par[idx.para],fixed=c(F,F),thres=thres,model="BR",maxit=1000,FUN=vario.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0),ub=c(Inf,1.99),hessian=TRUE,opt=FALSE,trace=TRUE,idx.para=1:2)

sqrt(diag(solve(results41$K) %*% results41$hessian %*% solve(results41$K)))

### jackknife results ###
data.avg = rowMeans(maxima.frechet)
data = maxima.frechet[data.avg>thres,]

results22 <- mcmapply(1:nrow(data),function(i){fit.model(data=data[-i,],loc=loc.sub.trans,init=results2$par,fixed=c(F,F,rep(F,n.alpha)),basis=basis,thres=thres,model="logskew",maxit=1000,FUN=vario.func,alpha.func=alpha.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0,rep(-Inf,n.alpha)),ub=c(Inf,1.99,rep(Inf,n.alpha)),hessian=FALSE,opt=TRUE,trace=TRUE,step2=FALSE,idx.para=1:2)},mc.cores=ncores,mc.set.seed = TRUE) 

results42 <- mcmapply(1:nrow(data),function(i){fit.model(data=data[-i,],loc=loc.sub.trans,init=results4$par[idx.para],fixed=c(F,F),thres=thres,model="BR",maxit=1000,FUN=vario.func,ncores=ncores,method="L-BFGS-B",lb=c(0.01,0.0),ub=c(Inf,1.99),hessian=FALSE,opt=TRUE,trace=TRUE,idx.para=1:2)},mc.cores=ncores,mc.set.seed = TRUE) 

save.image(file=paste0(DataPath,"data/application_results_new_.RData"))


