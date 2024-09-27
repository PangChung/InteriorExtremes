rm(list=ls())
args <- commandArgs(TRUE)
computer = "local"
id = 1
d <- 10 ## 10 * 10 grid on [0,1]^2
m <- 200 ## number of samples
model = "truncT"; # "logskew" or "truncT"
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
para.range = c(3,5) # range for the covariance function ##      
para.shape = c(1,1.5) #c(1,1.5) ## smoothness parameter for the covariance function ##
idx.para = 1:2 # variogram parameters; otherwise 1:3 for cov.func
para.deg = c(2,3) ## degree of the freedom for the truncated t model ##
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

########################################################################
### simulation study for the log-skew normal based max-stable process ##
########################################################################
rFun <- function(x){
    val = sum((x)^xi)^{1/xi}
    return(val)
}

lb=c(0.01,0.01,0)
ub=c(Inf,1.99,Inf)
fixed = c(F,F,T)
init = c(1,0.5,2)
par.truncT <- as.matrix(expand.grid(para.range,para.shape,para.deg))
set.seed(init.seed)

simu <- function(i){
    par.truncT.list <- list(sigma=cov.func(diff.mat,c(par.truncT[i,idx.para])),nu=par.truncT[i,-idx.para])    
    samples.truncT <- simu_Pareto_truncT(m=m,par=par.truncT.list,rFun,ncores=NULL)
    return(samples.truncT)
}

model.fit <- function(i){
    set.seed(init.seed)
    init[fixed] = par.truncT[i,fixed]
    data = samples.truncT[[i]]
    data.sum = apply(data,1,rFun)
    u = quantile(data.sum,0.95)
    data = data[data.sum>u,]/u
    fit.result2 <- fit.scoreMatching(init=init[-3],obs=data,loc=diff.mat,fixed=c(F,F), model="truncT",cov.func=cov.func,idx.para=idx.para,dof=par.truncT[i,3],weightFun = weightFun , dWeightFun = dWeightFun , method="Nelder-Mead", maxit=1000,lb=lb[-3],ub=ub[-3])
    fit.result1 <- fit.model(data=samples.truncT[[i]],loc=diff.mat,init=init,fixed=c(F,F,T),thres=u/ncol(data),model="truncT",FUN=cov.func,ncores=ncores,maxit=1000,method="Nelder-Mead",lb=lb,ub=ub,hessian=FALSE,opt=TRUE,trace=TRUE,idx.para=idx.para,pareto=TRUE)
    return(list(fit.result1,fit.result2))
}

if(file.exists(file.samples)){
    load(file.samples)
}else{
    samples.truncT <- mclapply(1:nrow(par.truncT),simu,mc.cores=ncores,mc.set.seed = TRUE)
    save(par.truncT,xi,samples.truncT,file=file.samples)
}

fit.truncT <- mclapply(1:nrow(par.truncT),model.fit,mc.cores=ncores,mc.set.seed = TRUE)
save(fit.truncT,xi,par.truncT,file=file2save)








