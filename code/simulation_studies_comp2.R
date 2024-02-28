args <- commandArgs(TRUE)
id = 1
computer = "hpc"
d <- 15 ## 10 * 10 grid on [0,1]^2
m <- 100 ## number of samples
# loading library and setting path
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
library(splines)
library(SpatialExtremes)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")

ncores=detectCores()
coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
para.range = c(0.5,1) ## range for the correlation function ##
para.nu = c(0.5,1) ## smoothness parameter for the correlation function ##
#para.alpha = rbind(c(0,0),c(-1,-2),c(-2,-1),c(2,1)) 
para.alpha = matrix(c(0,0),nrow=1)
para.deg = c(2,3) ## degree of the freedom for the truncated t model ##
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
file2save = paste0(DataPath,"data/simulation_study_comp2_",id,".RData")
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
set.seed(init.seed)
vecchia.seq <- 1:nrow(coord)#sample(1:nrow(coord),size=nrow(coord),replace=FALSE)
neighbours.mat <- sapply(1:nrow(coord),FUN=neighbours,vecchia.seq=vecchia.seq,
					q=3,loc=diff.mat)
lb=c(0.01,0.01,rep(-Inf,ncol(para.alpha)))
ub=c(10,2.0,rep(Inf,ncol(para.alpha)))
init = c(1,1,0,0)

##compute the basis ###
centers <- rbind(c(0.5,0.5),c(0.25,0.25),c(0.75,0.75))
idx.centers <- apply(centers,1,function(x){which.min(apply(coord,1,function(y){sum((x-y)^2)}))})
basis <- sapply(idx.centers,function(x){ y=dnorm(diff.mat[x,],mean=0,sd=1);y=y-mean(y) })
basis[,1] <- rep(0,nrow(coord))

pairs.idx = rank(diff.mat[t(all.pairs)]) < 4000
########################################################################
### simulation study for the log-skew normal based max-stable process ##
########################################################################
par.skew.normal <- as.matrix(expand.grid(para.range,para.nu,1:nrow(para.alpha)))
par.skew.normal <- cbind(par.skew.normal[,-3],para.alpha[par.skew.normal[,3],]);colnames(par.skew.normal) <- NULL
samples.skew.normal <- list()
par.skew.list <- list()
ec.logskew <- list()
tc.logskew <- list()
fit.logskew.vecchia <- list()
fit.logskew.comp <- list()
for(i in 1:nrow(par.skew.normal)){
    par.skew.list[[i]] <- list(sigma=cov.func(coord,par.skew.normal[i,1:2]),alpha=alpha.func(par=par.skew.normal[i,-c(1:2)]))
    samples.skew.normal[[i]] <- simu_logskew(m=m,par=alpha2delta(par.skew.list[[i]]),ncores=ncores)
    #fit.logskew.vecchia[[i]] <- MVLE(data=samples.skew.normal[[i]],init=par.skew.normal[i,],fixed=c(F,F,F,F,F),loc=coord,FUN=cov.func,vecchia.seq=vecchia.seq,neighbours = neighbours.mat,alpha.func=alpha.func,maxit=500,model="logskew",lb=lb,ub=ub,ncores=ncores)
    print(par.skew.normal[i,])
    
    # init = c(par.skew.normal[i,1:2],0,0)
    # fit.logskew.comp[[i]] <- MCLE(data=samples.skew.normal[[i]],init=init,fixed=c(F,F,T,T),loc=coord,FUN=cov.func,index=all.pairs[,pairs.idx],alpha.func=alpha.func,model="logskew",lb=lb,ub=ub,ncores=ncores,maxit=500,trace=FALSE)
    
    fit.logskew.comp[[i]] <- MCLE(data=samples.skew.normal[[i]],init=par.skew.normal[i,1:2],fixed=c(F,F),loc=coord,FUN=cov.func,index=all.pairs[,pairs.idx],model="BR",lb=lb[1:2],ub=ub[1:2],ncores=ncores,maxit=1000,trace=FALSE)

    fit.logskew.vecchia[[i]] <- MVLE(data=samples.skew.normal[[i]],init=par.skew.normal[i,1:2],fixed=c(F,F),loc=coord,FUN=cov.func,vecchia.seq=vecchia.seq,neighbours = neighbours.mat,maxit=1000,model="BR",lb=lb[1:2],ub=ub[1:2],ncores=ncores,trace=FALSE)
    fit.logskew.comp2[[i]] <- SpatialExtremes::fitmaxstab(data=samples.skew.normal[[i]],coord,cov.mod="brown",fit.marge = FALSE)
}

save(fit.logskew.comp,fit.logskew.vecchia,par.skew.normal,init.seed,file=file2save)

