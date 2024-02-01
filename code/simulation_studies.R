library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(gridExtra)
library(partitions)
library(ggplot2)
library(Rfast)
library(matrixStats)
library(splines)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")

ncores=detectCores()
d <- 10 ## 10 * 10 grid on [0,1]^2
m <- 10000 ## number of samples
set.seed(1342342)
coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
para.range = c(0.5,1,2) ## range for the correlation function ##
para.nu = c(0.5,1,1.5) ## smoothness parameter for the correlation function ##
para.alpha = rbind(c(0,0,0),c(-1,2,3),c(-2,-1,4)) ## slant parameter for skewed norm model ##
para.deg = 2 ## degree of the freedom for the truncated t model ##
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))

########################################################################
### simulation study for the log-skew normal based max-stable process ##
########################################################################
par.skew.normal <- as.matrix(expand.grid(para.range,para.nu,1:3))
par.skew.normal <- cbind(par.skew.normal[,-3],para.alpha[par.skew.normal[,3],]);colnames(par.skew.normal) <- NULL

samples.skew.normal <- list()
par.skew.list <- list()
ec.logskew <- list()
tc.logskew <- list()
for(i in 1:nrow(par.skew.normal)){
    par.skew.list[[i]] <- list(sigma=cov.func(coord,par.skew.normal[i,1:2]),alpha=alpha.func(coord,par.skew.normal[i,-c(1:2)]))
    samples.skew.normal[[i]] <- simu_logskew(m=m,par=alpha2delta(par.skew.list[[i]]),ncores=ncores)
    ec.logskew[[i]] <- lapply(all.pairs,empirical_extcoef,data=samples.skew.normal[[i]])
    tc.logskew[[i]] <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=alpha2delta(par.skew.list[[i]]),model="logskew1"),mc.cores=ncores,mc.set.seed=FALSE)
}



## trainning the model with the log-skew normal based max-stable process ##
idx.case = 11
fit.logskew.comp <- MCLE(data=samples.skew.normal[[idx.case]],init=c(0.5,1,0,0,0),fixed=c(F,F,F,F,F),loc=coord,FUN=cov.func,index=all.pairs[,sample(1:ncol(all.pairs),1000,replace=FALSE)],ncores=ncores,maxit=200,model="logskew",lb=c(0.1,0.1,-Inf,-Inf,-Inf),ub=c(10,2.5,Inf,Inf,Inf), alpha.func=alpha.func,hessian=TRUE)

vecchia.seq <- sample(1:nrow(coord),size=nrow(coord),replace=FALSE)
neighbours.mat <- sapply(1:nrow(coord),FUN=neighbours,vecchia.seq=vecchia.seq,
					q=3,loc=diff.mat)
fit.logskew.vecchia <- MVLE(data=samples.skew.normal[[idx.case]],init=c(0.5,1,0,0,0),fixed=c(F,F,F,F,F),loc=coord,FUN=cov.func,vecchia.seq=vecchia.seq,neighbours = neighbours.mat,alpha.func=alpha.func,maxit=200,model="logskew",lb=c(0.1,0.1,-Inf,-Inf,-Inf),ub=c(10,2,Inf,Inf,Inf),ncores=ncores)

fit.logskew.angular <- fit.model(data=samples.skew.normal[[idx.case]],loc=coord,init=c(0.4,0.8,0,0,0),fixed=c(F,F,F,F,F),thres=0.9,model="logskew",ncores=ncores,maxit=100,lb=c(0.01,0.01,-Inf,-Inf,-Inf),ub=c(10,2.0,Inf,Inf,Inf),bootstrap=FALSE,hessian=TRUE,opt=TRUE)
