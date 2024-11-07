args <- commandArgs(TRUE)
id = 1
computer = "local"
d <- 15
m <- 500 ## number of samples
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
para.range = c(3,6) ## range for the correlation function ##
para.nu = 8
para.shape = c(0.5,1,1.5) ## smoothness parameter for the correlation function ##
para.alpha = rbind(c(0,0,0)) 
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
file2save = paste0(DataPath,"data/simulation_study_comp_",id,".RData")
file.samples = paste0(DataPath,"data/samples/simulation_logskew_comp_",id,"_",m,".RData")
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
if(file.exists(file2save)){stop("job already finished")}
idx.para=1:2 # variogram parameters; otherwise 1:3 for cov.func
set.seed(init.seed)
vecchia.seq <- 1:nrow(coord) #sample(1:nrow(coord),size=nrow(coord),replace=FALSE)
neighbours.mat <- sapply(1:nrow(coord),FUN=neighbours,vecchia.seq=vecchia.seq,
					q=2,loc=diff.mat)
lb=c(0.01,0.01,rep(-Inf,ncol(para.alpha)))
ub=c(Inf,1.99,rep(Inf,ncol(para.alpha)))
init = c(2,1,0,0,0)

# lb=c(0.01,0.01,0.01,rep(-Inf,ncol(para.alpha)))
# ub=c(Inf,Inf,1.99,rep(Inf,ncol(para.alpha)))
# init = c(1,1,1,0,0)

pairs.idx = rank(diff.mat[t(all.pairs)]) < nrow(coord)*10

##compute the basis ###
basis = matrix(0,nrow=nrow(coord),ncol=3)

########################################################################
### simulation study for the log-skew normal based max-stable process ##
########################################################################
# par.skew.normal <- as.matrix(expand.grid(para.range,para.nu,para.shape,1:nrow(para.alpha)))
par.skew.normal <- as.matrix(expand.grid(para.range,para.shape,1:nrow(para.alpha)))
par.skew.normal <- cbind(par.skew.normal[,idx.para],para.alpha[par.skew.normal[,-idx.para],]);colnames(par.skew.normal) <- NULL

simu <- function(i){
    par.skew.list <- list(sigma=vario.func(coord,par.skew.normal[i,idx.para]))
    par.skew.list$alpha <- alpha.func(par=par.skew.normal[i,-idx.para],b.mat=basis)
    samples.skew.normal <- simu_Pareto_logskew(m=m,par=alpha2delta(par.skew.list),rFun,ncores=NULL)
    return(samples.skew.normal)
}

model.fit <- function(i){
    set.seed(init.seed)

    data = samples.skew.normal[[i]]
    data.sum = apply(data,1,sum)
    u = max(quantile(data.sum,0.95),nrow(coord)^(1-1/xi))
    data = data[data.sum>u,]/u

    fit.result <- fit.model(data=data,loc=coord,init=init,fixed=c(F,F,T,T,T),basis=basis,thres=0,model="logskew",FUN=vario.func,alpha.func=alpha.func,maxit=1000,method="L-BFGS-B",lb=lb,ub=ub,hessian=FALSE,opt=TRUE,trace=FALSE,idx.para=idx.para,pareto=TRUE,step2=FALSE,ncores=3)
    
    return(list(fit.result))
}

if(file.exists(file.samples)){load(file.samples,e<-new.env());samples.skew.normal<-e$samples.skew.normal}else{ 
    samples.skew.normal <- mclapply(1:nrow(par.skew.normal),simu,mc.cores=ncores,mc.set.seed = TRUE)
    save(samples.skew.normal,basis,par.skew.normal,xi,file=file.samples)
}

fit.logskew <- mclapply(1:nrow(par.skew.normal),model.fit,mc.cores=ncores,mc.set.seed = TRUE)

save(fit.logskew,par.skew.normal,basis,xi,file=file2save)

save(fit.logskew.comp,fit.logskew.angular,basis,par.skew.normal,init.seed,m,d,file=file2save)

if(!file.exists(file.samples)) save(samples.skew.normal,basis,coord,par.skew.normal,cov.func,alpha.func,file=file.samples)


