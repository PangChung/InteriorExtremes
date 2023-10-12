library(tmvnsim)
library(parallel)
library(evd)
library(gridExtra)
library(ggplot2)
source("code/simulation.R")
source("code/exponent_functions.R")
d <- 1000
coord = as.matrix(expand.grid(0,0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),
                         as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
corr <- function(x,r=1.5,v=0.3) exp(- (sqrt(sum(x^2))/r)^v)                          
cov.mat <- matrix(apply(diff.vector, 1, corr), ncol=nrow(coord)) 
chol(cov.mat)
nu = 2

par1 <- list(nu=nu,sigma=cov.mat)
system.time(Z.trunc <- bi.simu(m=10000,par=par1,ncores=10,model="truncT"))
ks.test(pgev(Z.trunc$val[[2]][2,],1,1,1),sort(runif(10000)))$p.value
ec.trunc <- unlist(lapply(Z.trunc$val,empirical_extcoef,idx=c(1,2)))
pairs <- rbind(1,2:d);pairs.list = split(pairs,col(pairs))
tc.truncT2 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=par1,model="truncT2"),mc.cores=10)

#pdf("figures/extcoef_truncT_bi.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(pairs)],y=ec.trunc,type="l",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main="Truncated extremal t processes",col="black")
lines(x=diff.mat[t(pairs)],y=tc.truncT2,cex=0.5,col="red")
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Theoretical"),col=c("black","red"),
    bty="n",lwd=1,cex=1)
#dev.off()

alpha = - 1 - coord[,2] + exp(sin(5*coord[,2]))
alpha = (alpha - min(alpha))/(max(alpha)-min(alpha))*2-1
par2.1 <- list(alpha=alpha,sigma=cov.mat)
alpha = 1 + 1.5*coord[,2] - exp(2*sin(10*coord[,2]))
alpha = (alpha - min(alpha))/(max(alpha)-min(alpha))*2-1
par2.2 <- list(alpha=alpha,sigma=cov.mat)
system.time(Z.logskew.1 <- bi.simu(m=10000,par=par2.1,ncores=10, model="logskew"))
system.time(Z.logskew.2 <- bi.simu(m=10000,par=par2.2,ncores=10, model="logskew"))
ec.logskew.1 <- unlist(lapply(Z.logskew.1$val,empirical_extcoef,idx=1:2))
ec.logskew.2 <- unlist(lapply(Z.logskew.2$val,empirical_extcoef,idx=1:2))
tc.logskew.1 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=par2.1,model="logskew2"),mc.cores=10)
tc.logskew.2 <- mcmapply(true_extcoef,pairs.list,MoreArgs=list(par=par2.2,model="logskew2"),mc.cores=10)
pdf("figures/extcoef_logskew_bi.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(pairs)],y=ec.logskew.1,type="p",cex=0.2,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main="Log-skew normal based max-stable processes",col="black",pch=20)
points(x=diff.mat[t(pairs)],y=ec.logskew.2,cex=0.2,ylim=c(1,2),col="grey",pch=20)
lines(x=diff.mat[t(pairs)],y=tc.logskew.1,cex=0.5,col="red",lty=1)
lines(x=diff.mat[t(pairs)],y=tc.logskew.2,cex=0.5,col="pink",lty=1)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical alpha 1","Empirical alpha 2"),col=c("black","grey"),
    bty="n",pch=20,cex=1)
legend("topright",legend=c("Theoretical alpha 1","Theoretical alpha 2"),col=c("red","pink"),
    bty="n",lwd=1,lty=1,cex=1)
dev.off()

save.image("data/bi_simulation.RData")
