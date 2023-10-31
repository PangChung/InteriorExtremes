library(tmvnsim)
library(parallel)
library(mvtnorm)
library(evd)
library(gridExtra)
library(partitions)
library(ggplot2)
library(Rfast)
library(partitions)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/MLE_BrownResnick.R")

### setting ###
d <- 10
coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),
                         as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
cov.mat <- cov.func(coord,c(0.5,1))
chol(cov.mat)
nu = 2
m = 1e+4
ncores=20
all.pairs <- combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))

### simulate the truncated extremal-t model ###
par1 <- list(nu=nu,sigma=cov.mat)
system.time(Z.trunc <- simu_truncT(m=m,par=par1,ncores=ncores))
which(apply(Z.trunc,1,anyDuplicated)>0)

png("figures/marginal_qqplot_truncT.png",width=d*1200,height=d*1200,res=200)
par(mfrow=c(d,d),mgp=c(2,1,0),mar=c(2,2,3,1))
y=(1:m)/(1+m)
for(idx in 1:nrow(Z.trunc)){
x = pgev(Z.trunc[idx,],1,1,1);
p.value = ks.test(x,y)$p.value
qqplot(x,y,main=paste0("Station ",idx, " P value ", round(p.value,2)),xlab="Data",ylab="Theoretical")
abline(0,1,col="red")
}
dev.off()
image(1:10,1:10,z=matrix(log(Z.trunc[,3]),nrow=10),col=rev(heat.colors(10)))

ec.trunc <- apply(all.pairs,2,empirical_extcoef,data=Z.trunc)
tc.truncT1 <- true_extcoef(all.pairs,par=par1,model="truncT1")

pdf("figures/extcoef_truncT.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.trunc,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main="Truncated extremal t processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=tc.truncT1,type="p",cex=0.5,col="red",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Theoretical"),col=c("black","red"),
    bty="n",pch=20,cex=1)
dev.off()

# Simulate a log-skew normal based max-stable process
alpha = alpha.func(coord,2)
#alpha = (alpha - min(alpha))/(max(alpha)-min(alpha))*10-5
par2 <- list(alpha=alpha,sigma=cov.mat)
system.time(Z.logskew <- simu_logskew(m=m,par=par2,ncores=ncores))
which(apply(Z.logskew,1,anyDuplicated)>0)

png("figures/marginal_qqplot_logskew.png",width=d*1200,height=d*1200,res=200)
par(mfrow=c(d,d),mgp=c(2,1,0),mar=c(2,2,3,1))
for(idx in 1:nrow(Z.logskew)){
z.order = order(Z.logskew[idx,],decreasing=FALSE)
x = pgev(Z.logskew[idx,z.order],1,1,1);y=sort(runif(length(x)))
p.value = ks.test(x,y)$p.value
qqplot(x,y,main=paste0("Station ",idx, " P value ", round(p.value,2)),xlab="Data",ylab="Theoretical")
abline(0,1,col="red")
}
dev.off()

image(1:10,1:10,z=matrix(log(Z.logskew[,1]),nrow=10),col=rev(heat.colors(10)) )


ec.logskew <- apply(all.pairs,2,empirical_extcoef,data=Z.logskew)
tc.logskew1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par2,model="logskew1"),mc.cores=ncores)

pdf("figures/extcoef_logskew.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=0.5,col="red",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Theoretical"),col=c("black","red"),
    bty="n", pch=20,cex=1)
dev.off()

## fit the model 
# fit the truncated extremal t model
system.time( fit.truncT <- fit.model(data=Z.trunc,loc=coord,init=fit.truncT$par,fixed=c(F,F,T),thres=0.9,model="truncT",ncores=ncores,maxit=500,lb=c(0.01,0.01),ub=c(10,2.0),bootstrap=FALSE,hessian=FALSE) )

# fit the log-skew based model
alpha = alpha.func(coord,-2)
#alpha = (alpha - min(alpha))/(max(alpha)-min(alpha))*10-5
par2 <- list(alpha=alpha,sigma=cov.mat)
system.time(Z.logskew <- simu_logskew(m=m,par=par2,ncores=ncores))
which(apply(Z.logskew,1,anyDuplicated)>0)
system.time( fit.logskew <- fit.model(data=Z.logskew,loc=coord,init=c(0.5,1,3),fixed=c(T,T,F),thres=0.95,model="logskew",method="Nelder-Mead",lb=c(0.1,0.1,-Inf),ub=c(10,1.9,Inf),bootstrap=FALSE,ncores=ncores,maxit=10000,hessian=TRUE) )

#system("say \'your program has finished\'")
#fit.result <- MCLE.BR(data=t(Z.logskew[1:10,1:100]),init=c(0.5,1),fixed=c(F,F),distmat=coord[1:10,],FUN = cov.func,index=combn(10,2),ncores=10,method="Nelder-Mead",maxit=1000,hessian=FALSE)
#idx.pairs <- which(all.pairs[1,]==45 | all.pairs[2,]==45)
cov.mat = cov.func(coord,fit.logskew$par[1:2])
alpha = alpha.func(coord,2)

fitted.extcoef.logskew1.1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=list(alpha=alpha,sigma=cov.mat),model="logskew1"),mc.cores=10)
plot(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=fitted.extcoef.logskew1.1,type="p",cex=0.5,col="red",pch=20)
boxplot(ec.logskew - fitted.extcoef.logskew1.1)
mean((fitted.extcoef.logskew1.1 - tc.logskew )^2)

cov.mat = cov.func(coord,fit.logskew$par[1:2])
alpha = alpha.func(coord,fit.logskew$par[-c(1:2)])
fitted.extcoef.logskew1.2 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=list(alpha=alpha,sigma=cov.mat),model="logskew1"),mc.cores=10)
plot(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=fitted.extcoef.logskew1.2,type="p",cex=0.5,col="red",pch=20)
boxplot(tc.logskew1 - fitted.extcoef.logskew1.2)
mean((fitted.extcoef.logskew1.2 - tc.logskew1 )^2)

plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=fitted.extcoef.logskew1,type="p",cex=0.5,col="red",pch=20)

save.image("data/simulation_test.RData")

system.time(for(i in 1:10) {x = combn(200,3)})
system.time(for(i in 1:10) {x = Rfast::comb_n(200,3)})