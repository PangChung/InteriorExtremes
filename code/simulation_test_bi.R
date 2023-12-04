rm(list=ls())
library(TruncatedNormal)
library(parallel)
library(evd)
library(gridExtra)
library(ggplot2)
library(partitions)
library(matrixStats)
library(mvtnorm)
source("code/likelihood_inference.R")
source("code/simulation.R")
source("code/exponent_functions.R")

## the spatial setting ##
d <- 1000
coord = as.matrix(expand.grid(0,0:(d-1))/d)  
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),
                         as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
cov.mat <- cov.func(coord,c(0.5,1))
pairs <- rbind(1,2:d);pairs.list = split(pairs,col(pairs))
chol(cov.mat)
nu = 2
m = 10000
random.seed = 34234
set.seed(random.seed)
ncores=10

## the truncatd extremal-t max-stable processes ##
par1 <- list(sigma=cov.mat,nu=nu)
system.time(Z.trunc <- bi.simu(m=m,par=par1,ncores=ncores,model="truncT"))
Z.trunc.val <- do.call(rbind,Z.trunc$val)
which(apply(Z.trunc.val,1,anyDuplicated)>0)

ec.trunc <- unlist(lapply(Z.trunc$val,empirical_extcoef,idx=c(1,2)))
tc.truncT2 <- mapply(true_extcoef,pairs.list,MoreArgs=list(par=par1,model="truncT2"))


pdf("figures/extcoef_truncT_bi.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(pairs)],y=ec.trunc,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main="Truncated extremal t processes",col="black",pch=20)
lines(x=diff.mat[t(pairs)],y=tc.truncT2,cex=0.5,col="red")
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Theoretical"),col=c("black","red"),
    bty="n",lwd=1,cex=1)
dev.off()

## log-skew normal based max-stable processes ##

alpha.range = 10
alpha = - 1 - coord[,2] + exp(sin(5*coord[,2]))
alpha.1 = (alpha - min(alpha))/(max(alpha)-min(alpha))*alpha.range-alpha.range/2
cov.mat.logskew = matrix(c(1,0.5,0.5,1),2,2)
par2.1 <- list(cov.mat.logskew,alpha.1)

alpha = 1 + 1.5*coord[,2] - exp(2*sin(10*coord[,2]))
alpha.2 = (alpha - min(alpha))/(max(alpha)-min(alpha))*alpha.range-alpha.range/2
par2.2 <- list(cov.mat.logskew,alpha.2)

alpha = 2.25 * sin(9*coord[,2])*cos(9*coord[,2])
alpha.3 = (alpha - min(alpha))/(max(alpha)-min(alpha))*alpha.range-alpha.range/2
par2.3 <- list(cov.mat.logskew,alpha.3)

alpha.4 = rep(0,d)
par2.4 <- list(cov.mat.logskew,alpha.4)

pdf("figures/alpha_bi.pdf",width=10,height=8)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.lab=1.5,cex=1.5,mgp=c(2,1,0))
plot(x=(1:d)/(d+1),y=alpha.1,type="l",ylim=c(-alpha.range+1,alpha.range+1),xlab="s",ylab="alpha",col="black",lwd=2)
lines(x=(1:d)/(d+1),y=alpha.2,type="l",col="red",lwd=2)
lines(x=(1:d)/(d+1),y=alpha.3,type="l",col="blue",lwd=2)
lines(x=(1:d)/(d+1),y=alpha.4,type="l",col="purple",lwd=2)
legend("topleft",legend=c("alpha 1","alpha 2","alpha 3","alpha 4"),col=c("black","red","blue","purple"),
    bty="n",lwd=2,cex=1)
dev.off()


set.seed(random.seed)
system.time(Z.logskew.1 <- bi.simu(m=m ,par=par2.1,ncores=ncores, model="logskew2"))

set.seed(random.seed)
system.time(Z.logskew.2 <- bi.simu(m=m,par=par2.2,ncores=ncores, model="logskew2"))

set.seed(random.seed)
system.time(Z.logskew.3 <- bi.simu(m=m,par=par2.3,ncores=ncores, model="logskew2"))

set.seed(random.seed)
system.time(Z.logskew.4 <- bi.simu(m=m,par=par2.4,ncores=ncores, model="logskew2"))

cov.mat.logskew = cov.func(cbind(1,seq(0,1,length.out=2)),c(1,1))
cov.mat.logskew = diag(c(2,2)) %*% cov.mat.logskew %*% diag(c(2,2))
alpha.vec = seq(-10,10,0.1)
val.1 = unlist(mclapply(alpha.vec,function(x) V_logskew(alpha.para=TRUE,par=list(cov.mat.logskew,c(x,x)),x=c(1,1)),mc.cores=ncores,mc.set.seed=FALSE))
plot(alpha.vec,val.1,pch=20,cex=0.5,ylim=c(1,2),xlab="x",ylab="Bivariate extremal coefficient",main="alpha=(x,x)")

val.2 = unlist(mclapply(alpha.vec,function(x) V_logskew(alpha.para=TRUE,par=list(cov.mat.logskew,c(-x,x)),x=c(1,1)),mc.cores=ncores,mc.set.seed=FALSE))
plot(alpha.vec,val.2,pch=20,cex=0.5,ylim=c(1,2),xlab="x",ylab="Bivariate extremal coefficient",main="alpha=(-x,x)",col=2)

val.3 = unlist(mclapply(alpha.vec,function(x) V_logskew(alpha.para=TRUE,par=list(cov.mat.logskew,c(-sin(x),sin(x))),x=c(1,1)),mc.cores=ncores,mc.set.seed=FALSE))
points(alpha.vec,val.3,pch=20,cex=0.5,ylim=c(1,2),xlab="x",ylab="Bivariate extremal coefficient",main="alpha=(-sin(x),sin(x))",col=3)

val.4 = unlist(mclapply(alpha.vec,function(x) V_logskew(alpha.para=TRUE,par=list(cov.mat.logskew,c(-sin(x),cos(x))),x=c(1,1)),mc.cores=ncores,mc.set.seed=FALSE))
points(alpha.vec,val.4,pch=20,cex=0.5,ylim=c(1,2),xlab="x",ylab="Bivariate extremal coefficient",main="alpha=(-sin(x),cos(x))",col=4)

val.5 = unlist(mclapply(alpha.vec,function(x) V_logskew(alpha.para=TRUE,par=list(cov.mat.logskew,c(sin(x),cos(x))),x=c(1,1)),mc.cores=ncores,mc.set.seed=FALSE))
points(alpha.vec,val.4,pch=20,cex=0.5,ylim=c(1,2),xlab="x",ylab="Bivariate extremal coefficient",main="alpha=(sin(x),cos(x))",col=5)

val.6 = unlist(mclapply(alpha.vec,function(x) V_logskew(alpha.para=TRUE,par=list(cov.mat.logskew,c(0,0)),x=c(1,1)),mc.cores=ncores,mc.set.seed=FALSE))

points(alpha.vec,val.6,pch=20,cex=0.5,ylim=c(1,2),xlab="x",ylab="Bivariate extremal coefficient",main="alpha=(0,0)",col=6)

val.7 = unlist(mclapply(alpha.vec,function(x) V_logskew(alpha.para=TRUE,par=list(cov.mat.logskew,c(x,10)),x=c(1,1)),mc.cores=ncores,mc.set.seed=FALSE))

points(alpha.vec,val.7,pch=20,cex=0.5,ylim=c(1,2),xlab="x",ylab="Bivariate extremal coefficient",main="alpha=(x,10)",col=7)

val.8 = unlist(mclapply(alpha.vec,function(x) V_logskew(alpha.para=TRUE,par=list(cov.mat.logskew,c(-10,x)*2),x=c(1,1)),mc.cores=ncores,mc.set.seed=FALSE))

points(alpha.vec,val.8,pch=20,cex=0.5,ylim=c(1,2),xlab="x",ylab="Bivariate extremal coefficient",main="alpha=(-10,x)",col=8)

val.9 = unlist(mclapply(alpha.vec,function(x) V_logskew(alpha.para=TRUE,par=list(cov.mat.logskew,c(-10,x)*2),x=c(1,1)),mc.cores=ncores,mc.set.seed=FALSE))

points(alpha.vec,val.8,pch=20,cex=0.5,ylim=c(1,2),xlab="x",ylab="Bivariate extremal coefficient",main="alpha=(-10,x)*2",col=9)


ec.logskew.1 <- unlist(lapply(Z.logskew.1$val,empirical_extcoef,idx=1:2))
ec.logskew.2 <- unlist(lapply(Z.logskew.2$val,empirical_extcoef,idx=1:2))
ec.logskew.3 <- unlist(lapply(Z.logskew.3$val,empirical_extcoef,idx=1:2))
ec.logskew.4 <- unlist(lapply(Z.logskew.4$val,empirical_extcoef,idx=1:2))

tc.logskew.1 <- mcmapply(true_extcoef,Z.logskew.1$par,MoreArgs=list(idx=1:2,model="logskew1"),mc.cores=ncores)
tc.logskew.2 <- mcmapply(true_extcoef,Z.logskew.2$par,MoreArgs=list(idx=1:2,model="logskew1"),mc.cores=ncores)
tc.logskew.3 <- mcmapply(true_extcoef,Z.logskew.3$par,MoreArgs=list(idx=1:2,model="logskew1"),mc.cores=ncores)
tc.logskew.4 <- mcmapply(true_extcoef,Z.logskew.4$par,MoreArgs=list(idx=1:2,model="logskew1"),mc.cores=ncores)

pdf("figures/extcoef_logskew_bi.pdf",width=10*4,height=8)
par(mfrow=c(1,4),mar=c(4,4,1,1),cex.main=1.5,cex.lab=2,cex=2,mgp=c(2.2,1,0))

plot(x=alpha.1,y=ec.logskew.1,type="p",cex=0.5,ylim=c(1,2),xlab="Alpha",ylab="Extremal coefficient",
    main="",col="black",pch=20)
points(x=alpha.1,y=tc.logskew.1,cex=0.5,col="red",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)

plot(x=alpha.2,y=ec.logskew.2,type="p",cex=0.5,ylim=c(1,2),xlab="Alpha",ylab="Extremal coefficient",
    main="",col="black",pch=20)
points(x=alpha.2,y=tc.logskew.2,cex=0.5,col="red",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)

plot(x=alpha.3,y=ec.logskew.3,type="p",cex=0.5,ylim=c(1,2),xlab="Alpha",ylab="Extremal coefficient",
    main="",col="black",pch=20)
points(x=alpha.3,y=tc.logskew.3,cex=0.5,col="red",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)

plot(x=alpha.4,y=ec.logskew.4,type="p",cex=0.5,ylim=c(1,2),xlab="Alpha",ylab="Extremal coefficient",
    main="",col="black",pch=20)
points(x=alpha.4,y=tc.logskew.4,cex=0.5,col="red",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
dev.off()

pdf("figures/delta_bi.pdf",width=10,height=8)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.lab=1.5,cex=1.5,mgp=c(2,1,0))
delta.1 = unlist(lapply(Z.logskew.1$par,function(x){x[[2]][1]}))

plot(x=alpha.1,y=delta.1,type="l",ylim=c(-1,1),xlab="alpha",ylab="delta",col="black",lwd=2)

delta.2 = unlist(lapply(Z.logskew.2$par,function(x){x[[2]][1]}))

lines(x=alpha.2,y=delta.2,type="l",col="red",lwd=2)

delta.3 = unlist(lapply(Z.logskew.3$par,function(x){x[[2]][1]}))

lines(x=alpha.3,y=delta.3,type="l",col="blue",lwd=2)

delta.4 = unlist(lapply(Z.logskew.4$par,function(x){x[[2]][1]}))

lines(x=alpha.4,y=delta.4,type="l",col="purple",lwd=2)

legend("topleft",legend=c("delta 1","delta 2","delta 3","delta 4"),col=c("black","red","blue","purple"),
    bty="n",lwd=2,cex=1)
dev.off()

## simulate the data for the d locations jointly ##
set.seed(random.seed)
system.time(Z.logskew.1.1 <- bi.simu(m=m ,par=alpha2delta(list(cov.mat,alpha.1)),ncores=ncores, model="logskew",alpha.para=FALSE))

set.seed(random.seed)
system.time(Z.logskew.2.1 <- bi.simu(m=m ,par=alpha2delta(list(cov.mat,alpha.2)),ncores=ncores, model="logskew",alpha.para=FALSE))

set.seed(random.seed)
system.time(Z.logskew.3.1 <- bi.simu(m=m ,par=alpha2delta(list(cov.mat,alpha.3)),ncores=ncores, model="logskew",alpha.para=FALSE))

set.seed(random.seed)
system.time(Z.logskew.4.1 <- bi.simu(m=m ,par=alpha2delta(list(cov.mat,alpha.4)),ncores=ncores, model="logskew",alpha.para=FALSE))

ec.logskew.1.1 <- unlist(lapply(Z.logskew.1.1$val,empirical_extcoef,idx=1:2))
ec.logskew.2.1 <- unlist(lapply(Z.logskew.2.1$val,empirical_extcoef,idx=1:2))
ec.logskew.3.1 <- unlist(lapply(Z.logskew.3.1$val,empirical_extcoef,idx=1:2))
ec.logskew.4.1 <- unlist(lapply(Z.logskew.4.1$val,empirical_extcoef,idx=1:2))

tc.logskew.1.1 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=alpha2delta(list(cov.mat,alpha.1)),model="logskew1"),mc.cores=ncores)
tc.logskew.2.1 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=alpha2delta(list(cov.mat,alpha.2)),model="logskew1"),mc.cores=ncores)
tc.logskew.3.1 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=alpha2delta(list(cov.mat,alpha.3)),model="logskew1"),mc.cores=ncores)
tc.logskew.4.1 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=alpha2delta(list(cov.mat,alpha.4)),model="logskew1"),mc.cores=ncores)

pdf("figures/extcoef_logskew_joint.pdf",width=10*4,height=8)
par(mfrow=c(1,4),mar=c(4,4,1,1),cex.main=1.5,cex.lab=2,cex=2,mgp=c(2.2,1,0))

plot(x=diff.mat[t(pairs)],y=ec.logskew.1.1,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal coefficient",
    main="",col="black",pch=20)
lines(x=diff.mat[t(pairs)],y=tc.logskew.1.1,cex=1,lwd=2,col="red",lty=1)
abline(h=c(1,2),col="grey",lty=2,cex=2)

plot(x=diff.mat[t(pairs)],y=ec.logskew.2.1,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal coefficient",
    main="",col="black",pch=20)
lines(x=diff.mat[t(pairs)],y=tc.logskew.2.1,cex=1,lwd=2,col="red",lty=1)
abline(h=c(1,2),col="grey",lty=2,cex=2)

plot(x=diff.mat[t(pairs)],y=ec.logskew.3.1,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal coefficient",
    main="",col="black",pch=20)
lines(x=diff.mat[t(pairs)],y=tc.logskew.3.1,cex=1,lwd=2,col="red",lty=1)
abline(h=c(1,2),col="grey",lty=2,cex=2)

plot(x=diff.mat[t(pairs)],y=ec.logskew.4.1,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal coefficient",
    main="",col="black",pch=20)
lines(x=diff.mat[t(pairs)],y=tc.logskew.4.1,cex=1,lwd=2,col="red",lty=1)
abline(h=c(1,2),col="grey",lty=2,cex=2)
dev.off()

pdf("figures/delta_joint.pdf",width=10,height=8)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.lab=1.5,cex=1.5,mgp=c(2,1,0))
plot(x=alpha.1,y=alpha2delta(list(cov.mat,alpha.1))[[2]],type="l",ylim=c(-1,1),xlab="alpha",ylab="delta",col="black",lwd=2)
lines(x=alpha.2,y=alpha2delta(list(cov.mat,alpha.2))[[2]],type="l",col="red",lwd=2)
lines(x=alpha.3,y=alpha2delta(list(cov.mat,alpha.3))[[2]],type="l",col="blue",lwd=2)
lines(x=alpha.4,y=alpha2delta(list(cov.mat,alpha.4))[[2]],type="l",col="purple",lwd=2)
legend("topleft",legend=c("delta 1","delta 2","delta 3","delta 4"),col=c("black","red","blue","purple"),
    bty="n",lwd=2,cex=1)
dev.off()

save.image(file="data/bi_simulation.RData")

## compute the likelihood for $alpha=0$ ## 
## it seems all the values agree with each other ##
par = Z.logskew.4$par[[1]]
data = Z.logskew.4$val[[1]]

val.1 = V(data,par[[1]])
val.1.1 = nVI(data,par[[1]],I = c(1,2))
val.1.2 = nVI(data,par[[1]],I = 1)
val.1.3 = nloglik(par,data,model="BR")
val.1.4 = nVI(data,par[[1]],I = 2)

val.2 = V_logskew(data,par,alpha.para=FALSE)
val.2.1 = intensity_logskew(data,par,alpha.para=FALSE,log=FALSE)
val.2.2 = partialV_logskew(data,idx = 1,par,alpha.para=FALSE)
val.2.3 = nloglik(par,data,model="logskew")
val.2.4 = partialV_logskew(data,idx = 2,par,alpha.para=FALSE)

max(abs(val.2.1 - val.1.1))
max(abs(val.2 - val.1))
max(abs(val.2.2 - val.1.2))
max(abs(val.2.3 - val.1.3))
max(abs(val.2.4 - val.1.4)^2)






