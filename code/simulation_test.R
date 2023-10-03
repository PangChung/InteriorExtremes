library(tmvnsim)
library(parallel)
library(evd)
library(doParallel)
source("code/simulation.R")
source("code/exponent_functions.R")
### testing the simulator ###
d <- 10
coord = as.matrix(expand.grid(1:d,1:d)/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),
                         as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
corr <- function(x,r=0.5,v=1) exp(- (sum(x^2)/r)^v)                          
cov.mat <- matrix(apply(diff.vector, 1, corr), ncol=nrow(coord)) + diag(1e-6,nrow(coord))       
chol(cov.mat)
nu = 2
par1 <- list(nu=nu,sigma=cov.mat)

# Simulate a truncated extremal-t max-stable process
system.time(Z.trunc <- simu_truncT(m=10000,par=par1,ncores=10))
hist(pgev(Z.trunc[,1],1,1,1),50,prob=TRUE)
image(1:10,1:10,z=matrix(log(Z.trunc[1,]),nrow=10),col=rev(heat.colors(10)))

# Simulate a log-skew normal based max-stable process
alpha = rep(0.5,nrow(coord))
par2 <- list(alpha=alpha,sigma=cov.mat)
system.time(Z.logskew <- simu_logskew(m=10000,par=par2,ncores=10))
hist(pgev(Z.logskew[,1],1,1,1),50,prob=TRUE)
image(1:10,1:10,z=matrix(log(Z.logskew[2,]),nrow=10),col=rev(heat.colors(10)) )

all.pairs <- combn(1:ncol(Z.trunc),2)
all.pairs.list = split(all.pairs,col(all.pairs))
ec.trunc <- apply(all.pairs,2,empirical_extcoef,data=Z.trunc)

tc.truncT1 <- true_extcoef(all.pairs,par=par1,model="truncT1")
tc.truncT2 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par1,model="truncT2"),mc.cores=10)

pdf("figures/extcoef_truncT.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.trunc,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main="Truncated extremal t processes",pch=20)
points(x=diff.mat[t(all.pairs)],y=tc.truncT1,type="p",cex=0.5,col="red",pch=20)
points(x=diff.mat[t(all.pairs)],y=tc.truncT2,type="p",cex=0.5,col="#0000ff72",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Emperical","Method 1","Method 2"),col=c("black","red","#0000ff72"),
    bty="n",pch=20,cex=1)
dev.off()

ec.logskew <- apply(all.pairs,2,empirical_extcoef,data=Z.logskew)
tc.logskew1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par2,model="logskew1"),mc.cores=10)
tc.logskew2 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par2,model="logskew2"),mc.cores=10)

pdf("figures/extcoef_logskew.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20)
points(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=0.5,col="#ff0000a8",pch=20)
points(x=diff.mat[t(all.pairs)],y=tc.logskew2,type="p",cex=0.5,col="#7eb3d8f2",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Emperical","Method 1","Method 2"),col=c("black","#ff0000a8","#7eb3d8f2"),
    bty="n", pch=20,cex=1)
dev.off()





