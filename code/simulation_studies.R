rm(list=ls())
args <- commandArgs(TRUE)
computer = "local"
id = 1
d <- 15 ## 10 * 10 grid on [0,1]^2
m <- 2000 ## number of samples
basis.idx = 1 # 1 for Gaussian Kernel and 2 for binary basis
model = "logskew"; # "logskew" or "truncT"
#model = "truncT"; # "logskew" or "truncT"
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
para.nu = 10 # ## variance parameter for the covariance function ##
para.shape = c(1,1.5) #c(1,1.5) ## smoothness parameter for the covariance function ##
idx.para = 1:2 # variogram parameters; otherwise 1:3 for cov.func
para.alpha = rbind(c(1,0,0),c(1,-1,-2),c(1,-1,1)) ## slant parameter for skewed norm model ##
para.deg = c(2,3) ## degree of the freedom for the truncated t model ##
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
thres = c(50,100)
if(model=="truncT"){thres=c(10,50,100);para.range = c(3,5)}
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
ncores=detectCores()
file2save = paste0(DataPath,"data/simulation_study2_",model,"_",id,"_",m,"_",basis.idx,".RData")
file.samples = paste0(DataPath,"data/samples/simulation_","model","_",id,"_",m,"_",basis.idx,".RData")
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

t0 <- proc.time()
if(model == "logskew"){
    # lb=c(0.01,0.01,0.01,rep(-Inf,ncol(para.alpha)))
    # ub=c(Inf,Inf,1.99,rep(Inf,ncol(para.alpha)))
    # init = c(1,para.nu,1,0,0)
    lb=c(0.01,0.01,rep(-Inf,ncol(para.alpha)))
    ub=c(Inf,1.99,rep(Inf,ncol(para.alpha)))
    init = c(1,1,0,0,0)
    # par.skew.normal <- as.matrix(expand.grid(para.range,para.nu,para.shape,1:nrow(para.alpha)))
    # par.skew.normal <- cbind(par.skew.normal[,idx.para],para.alpha[par.skew.normal[,-idx.para],]);colnames(par.skew.normal) <- NULL
    par.skew.normal <- as.matrix(expand.grid(para.range,para.shape,1:nrow(para.alpha)))
    par.skew.normal <- cbind(par.skew.normal[,idx.para],para.alpha[par.skew.normal[,-idx.para],]);colnames(par.skew.normal) <- NULL
    par.skew.list <- list()
    ec.logskew <- list()
    tc.logskew <- list()
    fit.logskew.angular <- list()
    if(file.exists(file.samples)){load(file.samples,e<-new.env());samples.skew.normal<-e$samples.skew.normal} else samples.skew.normal <- list()
    for(i in 1:nrow(par.skew.normal)){
        par.skew.list[[i]] <- list(sigma=vario.func(coord,par.skew.normal[i,idx.para]))
        par.skew.list[[i]]$alpha <- alpha.func(par=par.skew.normal[i,-idx.para],b.mat=basis)
        # par.skew.list[[i]] <- list(sigma=cov.func(diff.mat,par.skew.normal[i,idx.para]))
        # par.skew.list[[i]]$alpha <- alpha.func(par=par.skew.normal[i,-idx.para],b.mat=basis)
        if(!file.exists(file.samples)){
            samples.skew.normal[[i]] <- simu_logskew(m=m,par=alpha2delta(par.skew.list[[i]]),ncores=ncores)
        }
        # true.ext.coef <- unlist(mclapply(all.pairs.list,true_extcoef,par=alpha2delta(par.skew.list[[i]]),model="logskew1",mc.cores=ncores,mc.set.seed = FALSE))
        # print(range(true.ext.coef))
        # samples.skew.normal[[i]] <- simu_logskew(m=m,par=alpha2delta(par.skew.list[[i]]),ncores=ncores)
        # init = par.skew.normal[i,]
        # fit.result1 <- fit.model(data=samples.skew.normal[[i]],loc=diff.mat,init=init,fixed=c(F,T,F,F,F),basis=basis,thres=30,model="logskew",FUN=cov.func,alpha.func=alpha.func,ncores=ncores,maxit=1000,method="Nelder-Mead",lb=lb,ub=ub,hessian=FALSE,opt=TRUE,trace=FALSE,step2=TRUE,idx.para=idx.para)
        fit.result1 <- fit.model(data=samples.skew.normal[[i]],loc=coord,init=init,fixed=c(F,F,F,F,F),basis=basis,thres=30,model="logskew",FUN=vario.func,alpha.func=alpha.func,ncores=ncores,maxit=1000,method="Nelder-Mead",lb=lb,ub=ub,hessian=FALSE,opt=TRUE,trace=FALSE,step2=TRUE,idx.para=idx.para)
        print(fit.result1$par)
        print(par.skew.normal[i,])
        fit.logskew.angular[[i]] <- fit.result1
        print(i)
    }
    save(fit.logskew.angular,par.skew.normal,m,d,basis,file=file2save)
    if(!file.exists(file.samples)) save(samples.skew.normal,basis,coord,par.skew.normal,cov.func,alpha.func,file=file.samples)
}


print(t0 <- proc.time() - t0)

if(model == "truncT"){
    lb=c(0.01,0.01,0.01,0)
    ub=c(Inf,Inf,1.99,Inf)
    fixed = c(F,T,F,T)
    init = c(1,1,2,2)
    par.truncT <- as.matrix(expand.grid(para.range,para.nu,para.shape,para.deg))
    samples.truncT <- par.truncT.list <- ec.truncT  <- tc.truncT <- fit.truncT.angular <-  list()
    for(i in 1:nrow(par.truncT)){
        fit.truncT <- list()
        par.truncT.list[[i]] <- list(sigma=cov.func(coord,par.truncT[i,idx.para]),nu=par.truncT[i,-idx.para])
        set.seed(init.seed)
        samples.truncT[[i]] <- simu_truncT(m=m,par=par.truncT.list[[i]],ncores=ncores)
        init[!fixed] = par.truncT[i,!fixed]
        for(j in 1:length(thres)){
           fit.result1 <- fit.model(data=samples.skew.normal[[i]],loc=diff.mat,init=init,fixed=c(F,T,F,T),thres=30,model="truncT",FUN=cov.func,ncores=ncores,maxit=1000,method="Nelder-Mead",lb=lb,ub=ub,hessian=FALSE,opt=TRUE,trace=FALSE,step2=TRUE,idx.para=idx.para)
        }
        fit.truncT.angular[[i]] <- fit.truncT
        print(i)
    }
    save(fit.truncT.angular,par.truncT,file=file2save)
}




