library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(gridExtra)
library(partitions)
library(ggplot2)
library(Rfast)
library(matrixStats)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")

### setting ###
d <- 10
coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
cov.mat <- cov.func(coord,c(0.5,1))
chol(cov.mat)
nu = 2
m = 1e+4
ncores=detectCores()
all.pairs <- combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))

### simulate the truncated extremal-t model ###
par1 <- list(sigma=cov.mat,nu=nu)
system.time(Z.trunc <- simu_truncT(m=m,par=par1,ncores=10))

png("figures/marginal_qqplot_truncT.png",width=d*600,height=d*600,res=300)
par(mfrow=c(d,d),mgp=c(2,1,0),mar=c(2,2,3,1),cex=0.5)
y=(1:m)/(1+m)
for(idx in 1:ncol(Z.trunc)){
    print(idx)
    x = pgev(Z.trunc[,idx],1,1,1)
    #p.value = ks.test(x,y)$p.value
    qqplot(x,y,main=paste0("Station ",idx),xlab="Data",ylab="Theoretical")
    abline(0,1,col="red")
}
dev.off()

ec.trunc <- apply(all.pairs,2,empirical_extcoef,data=Z.trunc)
tc.truncT1 <- true_extcoef(all.pairs,par=par1,model="truncT1")

png("figures/extcoef_truncT.png",width=6*300,height=4*300,res=300)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.trunc,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
        main="Truncated extremal t processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=tc.truncT1,type="p",cex=0.5,col="red",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Theoretical"),col=c("black","red"),
    bty="n",pch=20,cex=1)
dev.off()

# Simulate a log-skew normal based max-stable process
alpha = alpha.func(coord,-3)
#par2 <- list(sigma=sd.mat %*% cov.mat %*% sd.mat,alpha=alpha)
par2 <- list(sigma=cov.mat,alpha=alpha)
system.time(Z.logskew <- simu_logskew(m=m,par=alpha2delta(par2),ncores=ncores))


png("figures/marginal_qqplot_logskew.png",width=d*600,height=d*600,res=300)
par(mfrow=c(d,d),mgp=c(2,1,0),mar=c(2,2,3,1),cex=0.5)
for(idx in 1:ncol(Z.logskew)){
z.order = order(Z.logskew[,idx],decreasing=FALSE)
x = pgev(Z.logskew[z.order,idx],1,1,1);y=sort(runif(length(x)))
qqplot(x,y,main=paste0("Station ",idx),xlab="Data",ylab="Theoretical")
abline(0,1,col="red")
}
dev.off()

ec.logskew <- apply(all.pairs,2,empirical_extcoef,data=Z.logskew)
tc.logskew1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=alpha2delta(par2),model="logskew1"),mc.cores=ncores,mc.set.seed=FALSE)
tc.logskew2 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par2,model="BR"),mc.cores=ncores,mc.set.seed=FALSE)

max(abs(tc.logskew1 - tc.logskew2))

png("figures/extcoef_logskew.png",width=6*300,height=4*300,res=300)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=1,col="red",pch=20)
points(x=diff.mat[t(all.pairs)],y=tc.logskew2,type="p",cex=1,col="blue",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Theoretical","BR"),col=c("black","red","blue"),
    bty="n", pch=20,cex=1)
dev.off()

######################
## fit the model #####
#####################
# fit the truncated extremal t model: the angular density approach
system.time( fit.truncT <- fit.model(data=Z.trunc,loc=coord,init=c(0.4,0.8,2),fixed=c(F,F,T),thres=0.9,model="truncT",method="Nelder-Mead",ncores=ncores,maxit=100,lb=c(0.01,0.01,-Inf),ub=c(10,2.0,Inf),bootstrap=FALSE,hessian=FALSE,opt=TRUE) )

## the composite likelihood approach
#fit.truncT.comp <- MCLE(data=Z.trunc,init=c(0.5,1,2),fixed=c(F,F,T),loc=coord,FUN=cov.func,index=all.pairs[,diff.mat[t(all.pairs)] < 0.3],maxit=200,model="truncT", lb=c(0.1,0.1,-Inf),ub=c(10,2.5,Inf),ncores=ncores)

## the vecchia approximation approach
#fit.trunct.vecchia <- MVLE(data=Z.trunc,init=c(0.5,1,2),fixed=c(F,T,T),loc=coord,FUN=cov.func,index=all.pairs[,diff.mat[t(all.pairs)]<0.5],maxit=200,model="truncT",lb=c(0.1,0.1,-Inf),ub=c(10,2.5,Inf),ncores=ncores)

# fit the log-skew based model
set.seed(422242)
alpha.para = 2
alpha = alpha.func(coord,alpha.para) 
cov.mat = cov.func(coord,c(0.5,1))
#cov.mat = diag(seq(1,2,length.out=d^2)) %*% cov.mat %*% diag(seq(1,2,length.out=d^2))
par2 <- list(sigma=cov.mat,alpha=alpha)
range(alpha2delta(par2)[[2]])
system.time(Z.logskew <- simu_logskew(m=m,par=alpha2delta(par2),ncores=ncores))
system.time(fit.logskew <- fit.model(data=Z.logskew,loc=coord,init=c(0.3,0.5,-4),fixed=c(F,F,F),thres=0.98,model="logskew",method="Nelder-Mead",lb=c(0.1,0.1,-Inf),ub=c(10,1.9,Inf),bootstrap=FALSE,ncores=ncores,maxit=10000,hessian=TRUE,opt=TRUE) )

## plot the intensity function ##
alpha.seq = seq(-5,5,0.1)
paras.list = as.matrix(expand.grid(fit.logskew$par[1],fit.logskew$par[2],alpha.seq))
paras.list = split(paras.list,row(paras.list))
paras.list.2 = lapply(paras.list,FUN=function(x){list(sigma=cov.func(coord,x[1:2]),alpha=alpha.func(coord,x[3]))})
logskew.vals <- lapply(paras.list,FUN=fit.model,data=Z.logskew,loc=coord,thres=0.98,model="logskew",ncores=ncores,opt=FALSE)
logskew.vals = sapply(logskew.vals,mean)
plot(alpha.seq,logskew.vals,cex=0.5,pch=20,,main="Log-skew normal based max-stable processes",xlab="alpha",ylab="Approxi. likelihood")
idx.min = which.min(logskew.vals)
points(alpha.seq[idx.min],logskew.vals[idx.min],cex=1,pch=20,col="red")
idx.est = which.min(abs(alpha.seq - fit.logskew$par[3]))
points(alpha.seq[idx.est],logskew.vals[idx.est],cex=1,pch=2,col="blue")
idx.true = which.min(abs(alpha.seq - alpha.para))
points(alpha.seq[idx.true],logskew.vals[idx.true],cex=1,pch=2,col="green")

system.time(logskew.comp.vals <- mclapply(paras.list.2,FUN=nlogcomplik,data=Z.logskew[1:1000,],index=all.pairs[,diff.mat[t(all.pairs)]<0.2],ncores=NULL,model="logskew",mc.cores=ncores,mc.set.seed = FALSE))
logskew.comp.vals <- sapply(logskew.comp.vals,mean)
plot(alpha.seq,logskew.comp.vals,cex=0.5,pch=20,,main="Log-skew normal based max-stable processes",xlab="alpha",ylab="Approxi. likelihood")
idx.min = which.min(logskew.comp.vals)
points(alpha.seq[idx.min],logskew.comp.vals[idx.min],cex=1,pch=20,col="red")
idx.true = which.min(abs(alpha.seq - alpha.para))
points(alpha.seq[idx.true],logskew.comp.vals[idx.true],cex=1,pch=2,col="green")

## fit the composite likelihood 
fit.logskew.comp <- MCLE(data=Z.logskew[1:1000,],init=c(0.5,1,-3),fixed=c(F,F,F),loc=coord,FUN=cov.func,index=all.pairs[,diff.mat[t(all.pairs)] < 0.5],ncores=ncores,maxit=200,model="logskew",lb=c(0.1,0.1,-Inf),ub=c(10,2.5,Inf),alpha.func=alpha.func,hessian=TRUE)

fit.logskew.comp.2 <- MCLE(data=Z.logskew[1:1000,],init=c(0.5,1),fixed=c(F,F),loc=coord,FUN=cov.func,index=all.pairs,ncores=ncores,maxit=200,model="BR",lb=c(0.1,0.1),ub=c(10,2.5),alpha.func=alpha.func,hessian=TRUE)

fit.logskew.comp3 <- MCLE(data=Z.logskew[1:100,],init=c(0.5,1,0.5),fixed=c(T,T,F),loc=coord,FUN=cov.func,index=all.pairs[,diff.mat[t(all.pairs)] < 0.3],ncores=ncores,maxit=200,model="logskew",lb=c(0.1,0.1,-5),ub=c(10,2.5,5),alpha.func=alpha.func,hessian=TRUE)

save.image("data/simulation_test.RData")

library(cubature)
library(SimplicialCubature)
n = 2
loc = cbind(seq(1,0,length.out=n),0)
sigma = cov.func(loc,c(0.5,1))
par = alpha2delta(list(sigma=sigma,alpha=rep(-5,n)))

func <- function(dat){
    val = exp(-nloglik(par=par,data=dat,model="logskew"))
}

func <- function(dat){
    val = exp(-nloglik(par=list(sigma),data=dat,model="BR"))
}

func <- function(dat){
    val = exp(-nloglik(par=list(sigma=sigma,nu=2),data=dat,model="truncT"))
}

print(res <- adaptIntegrate(func,rep(0,n),rep(Inf,n),tol=1e-5,maxEval=1e+5))

n = 10
loc = cbind(seq(1,0,length.out=n),0)
sigma = cov.func(loc,c(0.5,1))
par = alpha2delta(list(sigma=sigma,alpha=rep(-10,n)))

func <- function(dat){    
    dat = t(dat)
    if(!is.matrix(dat)) {data  <- c(dat,1-sum(dat))} else {data <- cbind(dat,1-rowSums(dat))}
    val <- intensity_logskew(data,list(sigma,rep(0,n)),alpha.para=TRUE,log=FALSE)
    return(val)
}

func <- function(dat){    
    dat = t(dat)
    if(!is.matrix(dat)) {data  <- c(dat,1-sum(dat))} else {data <- cbind(dat,1-rowSums(dat))}
    val <- nVI(data,sigma=sigma,I=1:n)
    return(val)
}

func <- function(dat){    
    dat = t(dat)
    if(!is.matrix(dat)) {data  <- c(dat,1-sum(dat))} else { data <- cbind(dat,1-rowSums(dat)) }
    val <- intensity_truncT(data,list(sigma,nu=2),log=FALSE)
    return(val)
}

S = CanonicalSimplex(n-1)
print(res <- adaptIntegrateSimplex(func,S,maxEvals=1e+6,absError=1e-4))

func <- function(x,idx,par){
    x = rbind(x,x) 
    x[1,idx] = x[1,idx] + 1e-3
    x[2,idx] = x[2,idx] - 1e-3
    val = V_truncT(x,par=par)
    val = (val[2]-val[1])/(2*1e-3)
    return(val)
}

x = c(1,2,3)
func(x,par=list(sigma=sigma,nu=2),idx=3)
partialV_truncT(x,par=list(sigma=sigma,nu=2),idx=3,log=FALSE)

func <- function(x,par,idx){
    x = rbind(x,x) 
    x[1,idx] = x[1,idx] + 1e-3
    x[2,idx] = x[2,idx] - 1e-3
    val = V_logskew(x,par=par,alpha.para=FALSE)
    val =(val[2]-val[1])/(2*1e-3)
    return(val)
}

func(x,idx=1,par=par)
partialV_logskew(x,par,idx=1,alpha.para=FALSE)

func <- function(x,par,idx){
    x = rbind(x,x) 
    x[1,idx] = x[1,idx] + 1e-8
    x[2,idx] = x[2,idx] - 1e-8
    val = partialV_logskew(x,par=par,alpha.para=FALSE,idx=(1:ncol(x))[-idx])
    val =(val[2]-val[1])/(2*1e-8)
    return(val)
}

func(x,idx=1,par=par)
intensity_logskew(x,par=par,alpha.para=FALSE,log=FALSE)
