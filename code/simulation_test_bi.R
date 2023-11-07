rm(list=ls())
library(tmvnsim)
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
#Z.trunc.gev <- apply(Z.trunc.val,1,function(x){pgev(x,1,1,1)})
#ks.test.result <- unlist(mclapply(split(Z.trunc.gev,col(Z.trunc.gev)),ks.test.new,mc.cores=8,mc.set.seed = TRUE))
#sum(!ks.test.result)

# for(i in 1:nrow(Z.trunc.val)){
# if(!ks.test.result[i]){
#     qqplot(Z.trunc.gev[,i],(1:m)/(1+m),cex=0.5,main=paste("p-value:",ks.test.val,"stat: ",i))
#     abline(0,1,col="red")
#     hist(Z.trunc.gev[,i],50)
#     print(i)
#     Sys.sleep(10)
# }
# }
ec.trunc <- unlist(lapply(Z.trunc$val,empirical_extcoef,idx=c(1,2)))
tc.truncT2 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=par1,model="truncT2"),mc.cores=ncores)

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
alpha = (alpha - min(alpha))/(max(alpha)-min(alpha))*alpha.range-alpha.range/2
par2.1 <- list(sigma=cov.mat,alpha=alpha)

alpha = 1 + 1.5*coord[,2] - exp(2*sin(10*coord[,2]))
alpha = (alpha - min(alpha))/(max(alpha)-min(alpha))*alpha.range-alpha.range/2
par2.2 <- list(sigma=cov.mat,alpha=alpha)

alpha = 2.25 * sin(9*coord[,2])*cos(9*coord[,2])
alpha = (alpha - min(alpha))/(max(alpha)-min(alpha))*alpha.range-alpha.range/2
par2.3 <- list(sigma=cov.mat,alpha=alpha)

alpha = rep(0,d)
par2.4 <- list(sigma=cov.mat,alpha=alpha)

pdf("figures/alpha_bi.pdf",width=10,height=8)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.lab=1.5,cex=1.5,mgp=c(2,1,0))
plot(x=(1:d)/(d+1),y=par2.1$alpha,type="l",ylim=c(-alpha.range+1,alpha.range+1),xlab="s",ylab="Alpha",col="black",lwd=2)
lines(x=(1:d)/(d+1),y=par2.2$alpha,type="l",col="red",lwd=2)
lines(x=(1:d)/(d+1),y=par2.3$alpha,type="l",col="blue",lwd=2)
legend("topleft",legend=c("Alpha 1","Alpha 2","Alpha 3"),col=c("black","red","blue"),
    bty="n",lwd=2,cex=1)
dev.off()

set.seed(random.seed)
system.time(Z.logskew.1 <- bi.simu(m=m ,par=par2.1,ncores=ncores, model="logskew"))
which(apply(do.call(rbind,Z.logskew.1$val),1,anyDuplicated)>0)

set.seed(random.seed)
system.time(Z.logskew.2 <- bi.simu(m=m,par=par2.2,ncores=ncores, model="logskew"))
which(apply(do.call(rbind,Z.logskew.2$val),1,anyDuplicated)>0)

set.seed(random.seed)
system.time(Z.logskew.3 <- bi.simu(m=m,par=par2.3,ncores=ncores, model="logskew"))
which(apply(do.call(rbind,Z.logskew.3$val),1,anyDuplicated)>0)

set.seed(random.seed)
system.time(Z.logskew.4 <- bi.simu(m=m,par=par2.4,ncores=ncores, model="logskew"))
which(apply(do.call(rbind,Z.logskew.3$val),1,anyDuplicated)>0)

ec.logskew.1 <- unlist(lapply(Z.logskew.1$val,empirical_extcoef,idx=1:2))
ec.logskew.2 <- unlist(lapply(Z.logskew.2$val,empirical_extcoef,idx=1:2))
ec.logskew.3 <- unlist(lapply(Z.logskew.3$val,empirical_extcoef,idx=1:2))
ec.logskew.4 <- unlist(lapply(Z.logskew.4$val,empirical_extcoef,idx=1:2))
tc.logskew.1 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=par2.1,model="logskew2"),mc.cores=ncores)
tc.logskew.2 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=par2.2,model="logskew2"),mc.cores=ncores)
tc.logskew.3 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=par2.3,model="logskew2"),mc.cores=ncores)
tc.logskew.4 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=par2.4,model="logskew2"),mc.cores=ncores)

pdf("figures/extcoef_logskew_bi.pdf",width=10*4,height=8)

par(mfrow=c(1,4),mar=c(4,4,1,1),cex.main=1.5,cex.lab=2,cex=2,mgp=c(2.2,1,0))
plot(x=diff.mat[t(pairs)],y=ec.logskew.1,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal coefficient",
    main="",col="black",pch=20)
lines(x=diff.mat[t(pairs)],y=tc.logskew.1,cex=1,lwd=2,col="red",lty=1)
abline(h=c(1,2),col="grey",lty=2,cex=2)

plot(x=diff.mat[t(pairs)],y=ec.logskew.2,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal coefficient",
    main="",col="black",pch=20)
lines(x=diff.mat[t(pairs)],y=tc.logskew.2,cex=1,lwd=2,col="red",lty=1)
abline(h=c(1,2),col="grey",lty=2,cex=2)

plot(x=diff.mat[t(pairs)],y=ec.logskew.3,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal coefficient",
    main="",col="black",pch=20)
lines(x=diff.mat[t(pairs)],y=tc.logskew.3,cex=1,lwd=2,col="red",lty=1)
abline(h=c(1,2),col="grey",lty=2,cex=2)

plot(x=diff.mat[t(pairs)],y=ec.logskew.4,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal coefficient",
    main="",col="black",pch=20)
lines(x=diff.mat[t(pairs)],y=tc.logskew.4,cex=1,lwd=2,col="red",lty=1)
abline(h=c(1,2),col="grey",lty=2,cex=2)
dev.off()

save.image(file="data/bi_simulation.RData")

## compute the likelihood ## 
## it seems all the values agrees with each other ##
par = Z.logskew.4$par[[1]]
data = Z.logskew.4$val[[1]]

val.1 = V(data,par[[1]])
val.1.1 = nVI(data,par[[1]],I = c(1,2))
val.1.2 = nVI(data,par[[1]],I = 1)
val.1.3 = nloglik(par,data,model="BR")
val.2 = V_logskew(data,par)
val.2.1 = intensity_logskew(data,par,log=FALSE)
val.2.2 = partialV_logskew(data,idx = 1,par)
val.2.3 = nloglik(par,data,model="logskew")
max((val.2.1 - val.1.1)^2)
max((val.2 - val.1)^2)
max((val.2.2 - val.1.2)^2)







