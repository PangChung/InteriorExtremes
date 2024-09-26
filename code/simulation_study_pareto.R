rm(list=ls())
args <- commandArgs(TRUE)
computer = "local"
id = 1
d <- 15 ## 10 * 10 grid on [0,1]^2
m <- 2000 ## number of samples
basis.idx = 1 # 1 for Gaussian Kernel and 2 for binary basis
model = "logskew"; # "logskew" or "truncT"
xi=3
for (arg in args) eval(parse(text = arg))
switch(computer,
    "ws" = {DataPath<-"~/Desktop/InteriorExtremes/"},
    "hpc" = {DataPath<-"/srv/scratch/z3536974/";.libPaths("../src")},
    "local" = {DataPath<-"~/Documents/Github/InteriorExtremes/"}
)
# settings 
coord = as.matrix(expand.grid(1:d,1:d))
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
para.range = c(5,10) # range for the covariance function ##      
para.shape = c(1,1.5) #c(1,1.5) ## smoothness parameter for the covariance function ##
idx.para = 1:2 # variogram parameters; otherwise 1:3 for cov.func
para.alpha = rbind(c(1,0,0),c(1,-1,-2),c(1,-1,1)) ## slant parameter for skewed norm model ##
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
# loading library and setting path
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(partitions)
library(Rfast)
library(matrixStats)
library(splines)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
source("code/pareto_inference.R")
ncores=detectCores()
file2save = paste0(DataPath,"data/simulation_pareto_",model,"_",id,"_",xi,".RData")
file.samples = paste0(DataPath,"data/samples/samples_pareto_",model,"_",id,"_",xi,".RData")
if(file.exists(file2save)){stop("job already finished")}
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
set.seed(init.seed)

## compute the basis ###
if(basis.idx == 1){
    centers <- rbind(c(0.25,0.25),c(0.25,0.25),c(0.75,0.75))*d
    idx.centers <- apply(centers,1,function(x){which.min(apply(coord,1,function(y){sum((x-y)^2)}))})
    basis <- sapply(idx.centers,function(x){y=dnorm(diff.mat[x,],mean=0,sd=d*2);y=y-mean(y);y/sqrt(sum(y^2))})
}else{
    idx = floor(matrix(seq(1,nrow(coord),length.out=6),ncol=2,3))
    basis <- sapply(1:(ncol(para.alpha)+1),function(x){y <- rep(0,nrow(coord));y[idx[x,]] <- c(-2,2);y})
}

basis[,1] = rep(0,d^2);basis[1:floor(d^2/2),1] = 0.1; basis[(d^2-floor(d^2/2)+1):d^2,1] = -0.1

########################################################################
### simulation study for the log-skew normal based max-stable process ##
########################################################################


rFun <- function(x){
    val = sum((x)^xi)^{1/xi}
    return(val)
}

lb=c(0.01,0.01,rep(-Inf,ncol(para.alpha)))
ub=c(Inf,1.99,rep(Inf,ncol(para.alpha)))
init = c(1,1,1,0,0)
# par.skew.normal <- as.matrix(expand.grid(para.range,para.nu,para.shape,1:nrow(para.alpha)))
# par.skew.normal <- cbind(par.skew.normal[,idx.para],para.alpha[par.skew.normal[,-idx.para],]);colnames(par.skew.normal) <- NULL
par.skew.normal <- as.matrix(expand.grid(para.range,para.shape,1:nrow(para.alpha)))
par.skew.normal <- cbind(par.skew.normal[,idx.para],para.alpha[par.skew.normal[,-idx.para],]);colnames(par.skew.normal) <- NULL

simu <- function(i){
    par.skew.list <- list(sigma=vario.func(coord,par.skew.normal[i,idx.para]))
    par.skew.list$alpha <- alpha.func(par=par.skew.normal[i,-idx.para],b.mat=basis)
    samples.skew.normal <- simu_Pareto_logskew(m=m,par=alpha2delta(par.skew.list),rFun,ncores=NULL)
    return(samples.skew.normal)
}

model.fit <- function(i){
    data = samples.skew.normal[[i]]
    data.sum = apply(data,1,rFun)
    u = quantile(data.sum,0.95)
    data = data[data.sum>u,]/u
    t0 <-  proc.time()
    fit.result1 <- fit.scoreMatching(init=init, obs=data, loc=coord, fixed=c(F,F,T,T,T), model="logskew", vario.func=vario.func, idx.para=idx.para, alpha.func=alpha.func, basis=basis, weightFun = weightFun, dWeightFun = dWeightFun, method="Nelder-Mead", maxit=1000, ncores = NULL)
    fit.result1 <- fit.scoreMatching(init=fit.result1$par, obs=data, loc=coord, fixed=c(F,F,T,F,F), model="logskew", vario.func=vario.func, idx.para=idx.para, alpha.func=alpha.func, basis=basis, weightFun = weightFun, dWeightFun = dWeightFun, method="Nelder-Mead", maxit=1000, ncores = NULL)
    t0 <- proc.time() - t0
    fit.result1$time <- t0

    fit.result2 <- fit.model(data=samples.skew.normal[[i]],loc=coord,init=init,fixed=c(F,F,T,F,F),basis=basis,thres=u/ncol(data),model="logskew",FUN=vario.func,alpha.func=alpha.func,ncores=NULL,maxit=1000,method="Nelder-Mead",lb=lb,ub=ub,hessian=FALSE,opt=TRUE,trace=FALSE,step2=TRUE,idx.para=idx.para,pareto=TRUE)

    return(list(fit.result1,fit.result2))
}

if(file.exists(file.samples)){load(file.samples,e<-new.env());samples.skew.normal<-e$samples.skew.normal}else{ 
    samples.skew.normal <- mclapply(1:nrow(par.skew.normal),simu,mc.cores=ncores,mc.set.seed = TRUE)
    save(samples.skew.normal,basis,par.skew.normal,xi,file=file.samples)
}
fit.logskew <- mclapply(1:nrow(par.skew.normal),model.fit,mc.cores=ncores,mc.set.seed = TRUE)
save(fit.logskew,par.skew.normal,basis,xi,file=file2save)

# d = 2
# loc = matrix(rnorm(d*2),ncol=2)*10

# alpha = rnorm(d)
# data = exp(rnorm(d))

# cov.mat = matrix(runif(d^2),ncol=d);cov.mat = (cov.mat + t(cov.mat))/2
# cov.mat = cov.mat + diag(d)*d

# delta=alpha2delta(list(cov.mat,alpha))[[2]]
# par = list(cov.mat,delta)
# intensity_logskew(data,list(cov.mat=cov.mat,delta=rep(0,d)),alpha.para=FALSE,log=FALSE)
# for(i in 1:d){
#     print(intensity_HR(data,par,i))
# }

# intensity_logskew(data,par,alpha.para=FALSE,log=FALSE)
# for(i in 1:d){
#     print(intensity_skewedHR(data,par,i))
# }

# V_logskew(rep(1,d),list(cov.mat,alpha=rep(0,d)),alpha.para=TRUE)
# for(i in 1:d){
#     print(V_HR(rep(1,d),par,i))
# }

# V_logskew(rep(1,d),par,alpha.para=FALSE)
# for(i in 1:d){
#     print(V_skewedHR(rep(1,d),par,i))
# }






