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
load("data/Trends_fits.RData")
data <- residuals2
data[data<0] = 0
data <- apply(data,2,function(x) {ind = x>0;x[ind] = qgpd(rank(x[ind])/(sum(ind)+1),1,1,1); x })
data = apply(data,1,function(x){list(which(x>0),x[x>0])})
idx.data = which(sapply(data,function(x) length(x[[2]]))>0)
data = data[idx.data]
data.sum = sapply(data,function(x) mean(x[[2]]))
data.max = sapply(data,function(x) max(x[[2]]))

# thres <- quantile(data.sum, seq(0.9,0.9999,length.out=30))
# mev::tstab.gpd(data.sum,thres)

# thres <- quantile(data.max, seq(0.9,0.9999,length.out=30))
# mev::tstab.gpd(data.max,thres)

data.fit.sum = data[data.sum>quantile(data.sum,0.99)]
data.fit.max = data[data.sum>quantile(data.max,0.99)]

D = nrow(loc.sub.trans)
ncores = floor(detectCores()/2)
idx.centers = unlist(lapply(quantile(loc.sub.trans[,1],seq(0.1,0.9,length.out=5)),function(x){ idx = abs(loc.sub.trans[,1] - x) < 5; which(idx)[which.min(abs(loc.sub.trans[idx,2] - median(loc.sub.trans[idx,2])))]}))

basis <- sapply(idx.centers,function(x){y=dnorm(distmat[x,],mean=0,sd=ncol(distmat)*2);y=y-mean(y);y/sqrt(sum(y^2))})

init = c(100,1,rep(0,ncol(basis)))
n.alpha = ncol(basis)
idx.para = 1:2
ub = c(Inf,1.99,rep(Inf,n.alpha))
lb = c(0.01,0.01,rep(-Inf,n.alpha))

file.save = paste0(DataPath,"/data/application_RedSea_results_",id,"_",method,"_",idx.jack,".RData")
file.origin = paste0(DataPath,"/data/application_RedSea_results_",id,"_",method,".RData")

if(file.exists(file.save)){
    stop("job already done")
}

if(idx.jack!=0){
    switch(id,
        {fit.result <- fit.model(data=data.fit.sum[-idx.jack],loc=loc.sub.trans,init=init,fixed=c(F,F,F,F,rep(T,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)},
        {fit.result <- fit.model(data=data.fit.max[-idx.jack],loc=loc.sub.trans,init=init,fixed=c(F,F,F,F,rep(T,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)},
        {t0 <- proc.time();
        fit.result <- fit.model(data=data.fit.sum[-idx.jack],loc=loc.sub.trans,init=init,fixed=c(F,F,F,F,rep(F,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)
        fit.result$time = proc.time() - t0},
        {fit.result <- fit.model(data=data.fit.max[-idx.jack],loc=loc.sub.trans,init=init,fixed=c(F,F,F,F,rep(F,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE);
        # fit.result <- fit.model(data=data.fit.max[-idx.jack],loc=loc.sub.trans,init=fit.result$par,fixed=c(T,T,T,T,rep(F,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE);
        # fit.result <- fit.model(data=data.fit.max[-idx.jack],loc=loc.sub.trans,init=fit.result$par,fixed=c(F,F,F,F,rep(T,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE);fit.result$time = proc.time() - t0
        })
        save(fit.result,idx.centers,basis,file=file.save)
}else{
     switch(id,
        {fit.result <- fit.model(data=data.fit.sum,loc=loc.sub.trans,init=init,fixed=c(F,F,F,F,rep(T,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)},
        {fit.result <- fit.model(data=data.fit.max,loc=loc.sub.trans,init=init,fixed=c(F,F,F,F,rep(T,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE)},
        {#t0 <- proc.time();
        fit.result <- fit.model(data=data.fit.sum,loc=loc.sub.trans,init=init,fixed=c(F,F,F,F,rep(F,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE);
        # fit.result <- fit.model(data=data.fit.sum,loc=loc.sub.trans,init=fit.result$par,fixed=c(T,T,T,T,rep(F,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE);
        # fit.result <- fit.model(data=data.fit.sum,loc=loc.sub.trans,init=fit.result$par,fixed=c(F,F,F,F,rep(T,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE);fit.result$time = proc.time() - t0
        },
        {#t0 <- proc.time();
        fit.result <- fit.model(data=data.fit.max,loc=loc.sub.trans,init=init,fixed=c(F,F,F,F,rep(F,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE);
        # fit.result <- fit.model(data=data.fit.max,loc=loc.sub.trans,init=fit.result$par,fixed=c(T,T,T,T,rep(F,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE);
        # fit.result <- fit.model(data=data.fit.max,loc=loc.sub.trans,init=fit.result$par,fixed=c(F,F,F,F,rep(T,ncol(basis))),model="logskew",maxit=1000,FUN=vario.func,basis=basis,alpha.func=alpha.func,ncores=ncores,method=method,lb=lb,ub=ub,opt=TRUE,idx.para=idx.para,pareto=TRUE,partial=TRUE,step2=FALSE,trace=TRUE);fit.result$time = proc.time() - t0
        })
        save(fit.result,idx.centers,basis,file=file.origin)
}
