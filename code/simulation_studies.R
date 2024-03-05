rm(list=ls())
args <- commandArgs(TRUE)
computer = "hpc"
id = 1
d <- 15 ## 10 * 10 grid on [0,1]^2
m <- 1000 ## number of samples
model = "logskew"; # "logskew" or "truncT"
#model = "truncT"; # "logskew" or "truncT"
for (arg in args) eval(parse(text = arg))
switch(computer,
    "ws" = {DataPath<-"~/Desktop/InteriorExtremes/"},
    "hpc" = {DataPath<-"/srv/scratch/z3536974/";.libPaths("../src")},
    "local" = {DataPath<-"~/Documents/Github/InteriorExtremes/"}
)
# settings 
coord = as.matrix(expand.grid(0:(d-1),0:(d-1)))
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
para.range = c(2,4) #c(0.5,1,2) ## range for the correlation function ##      
para.nu = c(1,1.5) #c(0.5,1,1.5) ## smoothness parameter for the correlation function ##
para.alpha = rbind(c(0,0),c(-1,-2),c(2,1)) ## slant parameter for skewed norm model ##
para.deg = 2 ## degree of the freedom for the truncated t model ##
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
thres = c(0.95,0.9)
# loading library and setting path
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
#library(gridExtra)
library(partitions)
#library(ggplot2)
library(Rfast)
library(matrixStats)
library(splines)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
ncores=detectCores()
file2save = paste0(DataPath,"data/simulation_study_",model,"_",id,"_",m,".RData")
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
set.seed(init.seed)

#compute the basis ###
# centers <- rbind(c(0.5,0.5),c(0.25,0.25),c(0.75,0.75))
# idx.centers <- apply(centers,1,function(x){which.min(apply(coord,1,function(y){sum((x-y)^2)}))})
# basis <- sapply(idx.centers,function(x){ y=dnorm(diff.mat[x,],mean=0,sd=1);y=y/max(abs(y));y=y-mean(y)})

basis <- sapply(1:(ncol(para.alpha)+1),function(x){y <- rep(0,nrow(coord));y[sample(1:nrow(coord),2)] <- c(-2,2);y})

########################################################################
### simulation study for the log-skew normal based max-stable process ##
########################################################################
t0 <- proc.time()
if(model == "logskew"){
    lb=c(0.01,0.01,rep(-Inf,ncol(para.alpha)))
    ub=c(10,1.99,rep(Inf,ncol(para.alpha)))
    init = c(1,1,0,0)
    par.skew.normal <- as.matrix(expand.grid(para.range,para.nu,1:3))
    par.skew.normal <- cbind(par.skew.normal[,-3],para.alpha[par.skew.normal[,3],]);colnames(par.skew.normal) <- NULL
    samples.skew.normal <- list()
    par.skew.list <- list()
    ec.logskew <- list()
    tc.logskew <- list()
    fit.logskew.angular <- list()
    fit.logskew.angular2 <- list()
    for(i in 1:nrow(par.skew.normal)){
        fit.logskew <- list()
        fit.logskew2 <- list()
        par.skew.list[[i]] <- list(sigma=cov.func(coord,par.skew.normal[i,1:2]),alpha=alpha.func(par=par.skew.normal[i,-c(1:2)]))
        set.seed(init.seed)
        samples.skew.normal[[i]] <- simu_logskew(m=m,par=alpha2delta(par.skew.list[[i]]),ncores=ncores)
        # ec.logskew[[i]] <- unlist(lapply(all.pairs.list,empirical_extcoef,data=samples.skew.normal[[i]]))
        # tc.logskew[[i]] <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=alpha2delta(par.skew.list[[i]]),model="logskew1"),mc.cores=ncores,mc.set.seed=FALSE)
        for(j in 1:length(thres)){
            fit.result1 <- fit.model(data=samples.skew.normal[[i]],loc=coord,init=init,fixed=c(F,F,F,F),thres=thres[j],model="logskew",ncores=ncores,maxit=1000,method="Nelder-Mead",lb=lb,ub=ub,bootstrap=FALSE,hessian=FALSE,opt=TRUE,trace=FALSE)
            a = matrix(rnorm(ncol(para.alpha)*ncores),nrow=ncores)*2
            a <- sweep(a,1,sqrt(rowSums(a^2)),FUN="/")
            init.mat = cbind(fit.result1$par[1],fit.result1$par[2],a)
            init.list = split(init.mat,row(init.mat))
            fit.result = mcmapply(FUN=fit.model,init=init.list,MoreArgs=list(data=samples.skew.normal[[i]],loc=coord,fixed=c(T,T,F,F),thres=thres[j],model="logskew",ncores=NULL,maxit=1000,method="Nelder-Mead",lb=lb,ub=ub,bootstrap=FALSE,hessian=FALSE,opt=TRUE,trace=FALSE),mc.set.seed = FALSE,mc.cores=ncores,SIMPLIFY = FALSE)
            results.mat <- matrix(unlist(lapply(fit.result,function(x){c(x$value)})),nrow=ncores,byrow=TRUE)
            dist.mat <- as.matrix(dist(results.mat))
            fit.logskew[[j]] = fit.result[[which.min(colSums(dist.mat))]]
            fit.logskew2[[j]] = fit.result 
            print(c(i,j))
        }
        fit.logskew.angular[[i]] <- fit.logskew
        fit.logskew.angular2[[i]] <- fit.logskew2
        print(i)
    }
    save(fit.logskew.angular,fit.logskew.angular2,par.skew.normal,file=file2save)
}

print(t0 <- proc.time() - t0)

if(model == "truncT"){
    lb=c(0.01,0.01,-Inf)
    ub=c(10,1.99,Inf)
    par.truncT <- as.matrix(expand.grid(para.range,para.nu,para.deg))
    samples.truncT <- par.truncT.list <- ec.truncT  <- tc.truncT <- fit.truncT.angular <-  list()
    for(i in 1:nrow(par.truncT)){
        fit.truncT <- list()
        par.truncT.list[[i]] <- list(sigma=cov.func(coord,par.truncT[i,1:2]),nu=par.truncT[i,3])
        set.seed(init.seed)
        samples.truncT[[i]] <- simu_truncT(m=m,par=par.truncT.list[[i]],ncores=ncores)
        # ec.truncT[[i]] <- unlist(lapply(all.pairs.list,empirical_extcoef,data=samples.truncT[[i]]))
        # tc.truncT[[i]] <- true_extcoef(all.pairs,par=par.truncT.list[[i]],model="truncT2")
        for(j in 1:length(thres)){
            fit.truncT[[j]] <- fit.model(data=samples.truncT[[i]],loc=coord,init=c(1,1,par.truncT.list[[i]]$nu),fixed=c(F,F,T),thres=thres[j],model="truncT",ncores=ncores,maxit=1000,lb=lb,ub=ub,method="Nelder-Mead",bootstrap=FALSE,hessian=FALSE,opt=TRUE,trace=FALSE)
        }
        fit.truncT.angular[[i]] <- fit.truncT
        print(i)
    }
    save(fit.truncT.angular,par.truncT,file=file2save)
}


# for(idx in 1:nrow(par.skew.normal)){
#     print(fit.logskew.angular[[idx]][[1]]$par - par.skew.normal[idx,])
# }