library(tmvnsim)
library(parallel)
library(mvtnorm)
library(evd)
library(gridExtra)
library(partitions)
library(ggplot2)
library(Rfast)
library(matrixStats)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/MLE_BrownResnick.R")

### setting ###
d <- 4
coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
cov.mat <- cov.func(coord,c(0.5,1))
chol(cov.mat)
nu = 2
m = 1e+4
ncores=10
all.pairs <- combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))

### simulate the truncated extremal-t model ###
par1 <- list(sigma=cov.mat,nu=nu)
system.time(Z.trunc <- simu_truncT(m=m,par=par1,ncores=ncores))
which(apply(Z.trunc,2,anyDuplicated)>0)

png("figures/marginal_qqplot_truncT.png",width=d*600,height=d*600,res=300)
par(mfrow=c(d,d),mgp=c(2,1,0),mar=c(2,2,3,1))
y=(1:m)/(1+m)
for(idx in 1:ncol(Z.trunc)){
print(idx)
x = pgev(Z.trunc[,idx],1,1,1)
p.value = ks.test(x,y)$p.value
qqplot(x,y,main=paste0("Station ",idx, " P value ", round(p.value,2)),xlab="Data",ylab="Theoretical")
abline(0,1,col="red")
}
dev.off()

image(1:10,1:10,z=matrix(log(Z.trunc[1,]),nrow=10),col=rev(heat.colors(10)))

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
alpha = alpha.func(coord,2)
#alpha = (alpha - min(alpha))/(max(alpha)-min(alpha))*10-5
par2 <- list(sigma=cov.mat,alpha=alpha)
system.time(Z.logskew <- simu_logskew(m=m,par=par2,ncores=ncores))

which(apply(Z.logskew,1,anyDuplicated)>0)

png("figures/marginal_qqplot_logskew.png",width=d*600,height=d*600,res=300)
par(mfrow=c(d,d),mgp=c(2,1,0),mar=c(2,2,3,1))
for(idx in 1:ncol(Z.logskew)){
z.order = order(Z.logskew[,idx],decreasing=FALSE)
x = pgev(Z.logskew[z.order,idx],1,1,1);y=sort(runif(length(x)))
p.value = ks.test(x,y)$p.value
qqplot(x,y,main=paste0("Station ",idx, " P value ", round(p.value,2)),xlab="Data",ylab="Theoretical")
abline(0,1,col="red")
}
dev.off()

image(1:10,1:10,z=matrix(log(Z.logskew[1,]),nrow=10),col=rev(heat.colors(10)) )


ec.logskew <- apply(all.pairs,2,empirical_extcoef,data=Z.logskew)
tc.logskew1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par2,model="logskew1"),mc.cores=ncores)

png("figures/extcoef_logskew.png",width=6*300,height=4*300,res=300)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=0.5,col="red",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Theoretical"),col=c("black","red"),
    bty="n", pch=20,cex=1)
dev.off()

######################
## fit the model #####
#######################
# fit the truncated extremal t model: the angular density approach
system.time( fit.truncT <- fit.model(data=Z.trunc,loc=coord,init=c(0.1,0.5,2),fixed=c(F,F,T),thres=0.9,model="truncT",ncores=ncores,maxit=500,lb=c(0.01,0.01),ub=c(10,2.0),bootstrap=FALSE,hessian=FALSE) )
# the composite likelihood approach : 
fit.truncT.comp <- MCLE(data=Z.trunc,init=c(0.1,0.5,2),fixed=c(F,F,T),loc=coord,FUN=cov.func,index=all.pairs,maxit=200,model="truncT",
                lb=c(0.1,0.1,-Inf),ub=c(10,2.5,Inf),ncores=ncores)


# fit the log-skew based model
system.time( fit.logskew <- fit.model(data=Z.logskew,loc=coord,init=c(0.5,1,3),fixed=c(T,T,F),thres=0.95,model="logskew",method="Nelder-Mead",lb=c(0.1,0.1,-Inf),ub=c(10,1.9,Inf),bootstrap=FALSE,ncores=ncores,maxit=10000,hessian=TRUE) )
fit.truncT.comp <- MCLE(data=Z.logskew,init=c(0.1,0.5,2),fixed=c(F,F,F),loc=coord,FUN=cov.func,index=all.pairs,maxit=200,model="logskew",
                lb=c(0.1,0.1,-Inf),ub=c(10,2.5,Inf),alpha.func=alpha.func,ncores=ncores,method="Nelder-Mead",hessian=FALSE

#system("say \'your program has finished\'")
#fit.result <- MCLE.BR(data=t(Z.logskew[1:10,1:100]),init=c(0.5,1),fixed=c(F,F),distmat=coord[1:10,],FUN = cov.func,index=combn(10,2),ncores=10,method="Nelder-Mead",maxit=1000,hessian=FALSE)
#idx.pairs <- which(all.pairs[1,]==45 | all.pairs[2,]==45)
cov.mat = cov.func(coord,fit.logskew$par[1:2])
alpha = alpha.func(coord,2)

fitted.extcoef.logskew1.1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=list(sigma=cov.mat,alpha=alpha),model="logskew1"),mc.cores=10)
plot(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=fitted.extcoef.logskew1.1,type="p",cex=0.5,col="red",pch=20)
boxplot(ec.logskew - fitted.extcoef.logskew1.1)
mean((fitted.extcoef.logskew1.1 - tc.logskew )^2)

cov.mat = cov.func(coord,fit.logskew$par[1:2])
alpha = alpha.func(coord,fit.logskew$par[-c(1:2)])
fitted.extcoef.logskew1.2 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=list(sigma=cov.mat,alpha=alpha),model="logskew1"),mc.cores=10)
plot(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=fitted.extcoef.logskew1.2,type="p",cex=0.5,col="red",pch=20)
boxplot(tc.logskew1 - fitted.extcoef.logskew1.2)
mean((fitted.extcoef.logskew1.2 - tc.logskew1 )^2)

plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=fitted.extcoef.logskew1,type="p",cex=0.5,col="red",pch=20)

save.image("data/simulation_test.RData")
