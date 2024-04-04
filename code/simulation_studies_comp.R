args <- commandArgs(TRUE)
id = 1
computer = "ws"
d <- 10 ## 10 * 10 grid on [0,1]^2
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
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")

ncores=detectCores()
coord = as.matrix(expand.grid(1:d,1:d))
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
para.range = c(4,8) ## range for the correlation function ##
para.nu = c(0.5,1,1.5)*4 ## smoothness parameter for the correlation function ##
para.alpha = rbind(c(0,0)) 
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
file2save = paste0(DataPath,"data/simulation_study_comp_",id,".RData")
file.samples = paste0(DataPath,"data/samples/simulation_logskew_comp_",id,"_",m,".RData")
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
set.seed(init.seed)
vecchia.seq <- 1:nrow(coord) #sample(1:nrow(coord),size=nrow(coord),replace=FALSE)
neighbours.mat <- sapply(1:nrow(coord),FUN=neighbours,vecchia.seq=vecchia.seq,
					q=2,loc=diff.mat)
lb=c(0.01,0.01,rep(-Inf,ncol(para.alpha)))
ub=c(Inf,Inf,rep(Inf,ncol(para.alpha)))
init = c(1,1,0,0)
pairs.idx = rank(diff.mat[t(all.pairs)]) < 2000

##compute the basis ###
basis = matrix(0,nrow=nrow(coord),ncol=3)

########################################################################
### simulation study for the log-skew normal based max-stable process ##
########################################################################
par.skew.normal <- as.matrix(expand.grid(para.range,para.nu,1:nrow(para.alpha)))
par.skew.normal <- cbind(par.skew.normal[,-3],para.alpha[par.skew.normal[,3],]);colnames(par.skew.normal) <- NULL
if(file.exists(file.samples)){load(file.samples)} else samples.skew.normal <- list()
par.skew.list <- list()
ec.logskew <- list()
tc.logskew <- list()
fit.logskew.vecchia <- list()
fit.logskew.angular <- list()
fit.logskew.comp <- list()

#if(file.exists(file2save)){
#    load(file2save)
for(i in 1:nrow(par.skew.normal)){
    par.skew.list[[i]] <- list(sigma=cov.func(diff.mat,par.skew.normal[i,1:2]))
    par.skew.list[[i]]$alpha <- alpha.func(par=par.skew.normal[i,-c(1:2)],b.mat=basis / sqrt(diag(par.skew.list[[i]]$sigma)))
    if(!file.exists(file.samples)){
        samples.skew.normal[[i]] <- simu_logskew(m=m,par=alpha2delta(par.skew.list[[i]]),ncores=ncores)
    }
    fit.logskew.angular[[i]] <- fit.model(data=samples.skew.normal[[i]],init=init,fixed=c(F,F,T,T),loc=diff.mat,FUN=cov.func,alpha.func=alpha.func,model="logskew",lb=lb,ub=ub,ncores=ncores,maxit=1000,trace=FALSE,method="Nelder-Mead",opt=TRUE,hessian=FALSE,basis=basis)
    fit.logskew.comp[[i]] <- MCLE(data=samples.skew.normal[[i]],init=init,fixed=c(F,F,T,T),loc=diff.mat,FUN=cov.func,index=all.pairs[,pairs.idx],alpha.func=alpha.func,model="logskew",lb=lb,ub=ub,ncores=ncores,maxit=1000,trace=TRUE,basis=basis)
}
save(fit.logskew.comp,fit.logskew.angular,basis,par.skew.normal,init.seed,m,d,file=file2save)
#}
if(!file.exists(file.samples)) save(samples.skew.normal,basis,coord,par.skew.normal,cov.func,alpha.func,file=file.samples)
