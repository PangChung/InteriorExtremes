args <- commandArgs(TRUE)
id = 1
computer = "local"
d <- 15## 10 * 10 grid on [0,1]^2
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
para.alpha = rbind(c(0,0)) 
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
file2save = paste0(DataPath,"data/simulation_study_comp_",id,".RData")
file.samples = paste0(DataPath,"data/samples/simulation_logskew_comp_",id,"_",m,".RData")
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
idx.para=1:2 # variogram parameters; otherwise 1:3 for cov.func
set.seed(init.seed)
vecchia.seq <- 1:nrow(coord) #sample(1:nrow(coord),size=nrow(coord),replace=FALSE)
neighbours.mat <- sapply(1:nrow(coord),FUN=neighbours,vecchia.seq=vecchia.seq,
					q=2,loc=diff.mat)
lb=c(0.01,0.01,rep(-Inf,ncol(para.alpha)))
ub=c(Inf,1.99,rep(Inf,ncol(para.alpha)))
init = c(2,1,0,0)
lb=c(0.01,0.01,0.01,rep(-Inf,ncol(para.alpha)))
ub=c(Inf,Inf,1.99,rep(Inf,ncol(para.alpha)))
init = c(1,1,1,0,0)
pairs.idx = rank(diff.mat[t(all.pairs)]) < nrow(coord)*10

##compute the basis ###
basis = matrix(0,nrow=nrow(coord),ncol=3)
basis[,2] = rep(0,d^2);basis[1:floor(d^2/2),2] = 0.1; basis[(d^2-floor(d^2/2)+1):d^2,2] = -0.1

########################################################################
### simulation study for the log-skew normal based max-stable process ##
########################################################################
# par.skew.normal <- as.matrix(expand.grid(para.range,para.nu,para.shape,1:nrow(para.alpha)))
par.skew.normal <- as.matrix(expand.grid(para.range,para.shape,1:nrow(para.alpha)))
par.skew.normal <- cbind(par.skew.normal[,idx.para],para.alpha[par.skew.normal[,-idx.para],]);colnames(par.skew.normal) <- NULL
if(file.exists(file.samples)){load(file.samples,e<-new.env());samples.skew.normal<-e$samples.skew.normal} else samples.skew.normal <- list()
par.skew.list <- list()
ec.logskew <- list()
tc.logskew <- list()
fit.logskew.vecchia <- list()
fit.logskew.angular <- list()
fit.logskew.comp <- list()

for(i in 1:nrow(par.skew.normal)){
    # par.skew.list[[i]] <- list(sigma=cov.func(diff.mat,par.skew.normal[i,idx.para]))
    par.skew.list[[i]] <- list(sigma=vario.func(coord,par.skew.normal[i,idx.para])) #semivarigoram
    par.skew.list[[i]]$alpha <- alpha.func(par=par.skew.normal[i,-idx.para],b.mat=basis)
    # true.ext.coef <- unlist(mclapply(all.pairs.list,true_extcoef,par=alpha2delta(par.skew.list[[i]]),model="logskew1",mc.cores=ncores,mc.set.seed = FALSE))
    #print(range(true.ext.coef))#}
    if(!file.exists(file.samples)){
        samples.skew.normal[[i]] <- simu_logskew(m=m,par=alpha2delta(par.skew.list[[i]]),ncores=ncores)
    }
    # system.time(samples.skew.normal[[i]] <- simu_logskew(m=m,par=alpha2delta(par.skew.list[[i]]),ncores=ncores))
    # init = par.skew.normal[i,]
    # fit.logskew.angular[[i]] <- fit.model(data=samples.skew.normal[[i]],init=init,fixed=c(F,T,F,T,T),loc=diff.mat,FUN=cov.func,alpha.func=alpha.func,thres=50,model="logskew",lb=lb,ub=ub,ncores=ncores,maxit=1000,trace=FALSE,method="Nelder-Mead",opt=TRUE,hessian=FALSE,basis=basis,idx.para=idx.para,step2=FALSE)
    fit.logskew.angular[[i]] <- fit.model(data=samples.skew.normal[[i]],init=init,fixed=c(F,F,T,T),loc=coord,FUN=vario.func,alpha.func=alpha.func,thres=30,model="logskew",lb=lb,ub=ub,ncores=ncores,maxit=1000,trace=FALSE,method="Nelder-Mead",opt=TRUE,hessian=FALSE,basis=basis,idx.para=idx.para,step2=FALSE)
    #fit.logskew.comp[[i]] <- MCLE(data=samples.skew.normal[[i]][1:100,],init=init,fixed=c(F,F,F,T,T),loc=diff.mat,FUN=vario.func,index=all.pairs[,pairs.idx],alpha.func=alpha.func,model="logskew",lb=lb,ub=ub,ncores=ncores,maxit=1000,trace=TRUE,basis=basis,idx.para=idx.para)
    print(fit.logskew.angular[[i]]$par-par.skew.normal[i,])

    fit.logskew.comp[[i]] <- MCLE(data=samples.skew.normal[[i]],init=init,fixed=c(F,F,T,T),loc=coord,FUN=vario.func,index=all.pairs[,pairs.idx],alpha.func=alpha.func,model="logskew",lb=lb,ub=ub,ncores=ncores,maxit=1000,trace=FALSE,basis=basis,idx.para=idx.para)

    print(fit.logskew.comp[[i]]$par - par.skew.normal[i,])
}

    

save(fit.logskew.comp,fit.logskew.angular,basis,par.skew.normal,init.seed,m,d,file=file2save)

if(!file.exists(file.samples)) save(samples.skew.normal,basis,coord,par.skew.normal,cov.func,alpha.func,file=file.samples)

# nu = 4
# data <- simu_logskew(m=1000,par=alpha2delta(list(cov.func(diff.mat,c(2,nu,1)),alpha.func(c(1,0),basis))),ncores=ncores)
# #init = cbind(par.skew.normal[i,1],seq(1,30,length.out=100),par.skew.normal[i,3])
# init = cbind(seq(1,20,length.out=100),nu,1,1,0)
# result2 <- apply(init, 1, function(x) fit.model(data=data,init=x,fixed=c(F,F,F),loc=diff.mat,FUN=cov.func,thres=100,model="BR",lb=lb[idx.para],ub=ub[idx.para],ncores=ncores,maxit=1000,trace=FALSE,method="Nelder-Mead",opt=FALSE,hessian=FALSE,idx.para=idx.para))
# result2 <- apply(init, 1, function(x) fit.model(data=data,init=x,fixed=c(F,F,F,F,F),loc=diff.mat,FUN=cov.func,alpha.func=alpha.func,basis=basis,thres=100,model="logskew",lb=lb,ub=ub,ncores=ncores,maxit=1000,trace=FALSE,method="Nelder-Mead",opt=FALSE,hessian=FALSE,idx.para=idx.para))
# #result1 <- fit.model(data=samples.skew.normal[[i]],init=init,fixed=c(F,T,F,T,T),loc=diff.mat,FUN=cov.func,alpha.func=alpha.func,thres=50,model="logskew",lb=lb,ub=ub,ncores=ncores,maxit=1000,trace=FALSE,method="Nelder-Mead",opt=TRUE,hessian=FALSE,basis=basis,idx.para=idx.para,step2=FALSE)

# plot(init[,1],result2,xlab="range",ylab="values",main=paste("Log-likelihood (BR):",paste(c(4,nu,1),collapse = ",")),type="l")
# init[which.min(result2),1]

