library(tmvnsim)
library(parallel)
library(mvtnorm)
library(evd)
library(gridExtra)
library(partitions)
library(ggplot2)
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

### simulate the truncated extremal-t model ###
par1 <- list(nu=nu,sigma=cov.mat)

Z.T.1 = TruncatedNormal::rtmvt(n=m,mu=cov.mat[,1],sigma=cov.mat,lb=rep(0,nrow(coord)),ub=rep(Inf,nrow(coord)),df=nu+1)

Z.T.2 <- tmvtsim(m,nrow(coord),lower=rep(0,nrow(coord)),upper=rep(Inf,nrow(coord)),means=cov.mat[,1],sigma=cov.mat,df=nu+1)[[1]]

idx = 35
hist((Z.T.1[,idx]-cov.mat[idx,1]),50,prob=TRUE)
abline(v=0,col="red")

hist((Z.T.2[idx,]-cov.mat[idx,1]),50,prob=TRUE)
abline(v=0,col="red")

system.time(Z.trunc <- simu_truncT(m=m,par=par1,ncores=10))
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

# Simulate a log-skew normal based max-stable process
alpha = alpha.func(coord,2)
#alpha = (alpha - min(alpha))/(max(alpha)-min(alpha))*10-5
par2 <- list(alpha=alpha,sigma=cov.mat)
system.time(Z.logskew <- simu_logskew(m=m,par=par2,ncores=10))
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

all.pairs <- combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))

ec.trunc <- apply(all.pairs,2,empirical_extcoef,data=Z.trunc)
tc.truncT1 <- true_extcoef(all.pairs,par=par1,model="truncT1")

ec.logskew <- apply(all.pairs,2,empirical_extcoef,data=Z.logskew)
tc.logskew1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par2,model="logskew1"),mc.cores=10)

idx = which.min(apply(diff.mat,2,sum))
idx.pairs = which(all.pairs[1,] == idx | all.pairs[2,] == idx)

pdf("figures/extcoef_truncT.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.trunc,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main="Truncated extremal t processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=tc.truncT1,type="p",cex=0.5,col="red",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Theoretical"),col=c("black","red"),
    bty="n",pch=20,cex=1)
dev.off()

pdf("figures/extcoef_logskew.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=0.5,col="red",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Theoretical"),col=c("black","red"),
    bty="n", pch=20,cex=1)
dev.off()

## plot the bivariate extremal coefficients map with respect to the central point ##
#for(idx.min in 1:100){
#idx.min = which.min(apply(diff.mat,2,sum))
idx.ind.1 = which(all.pairs[1,]==idx.min)
idx.ind.2 = which(all.pairs[2,]==idx.min)
idx.ind = c(idx.ind.1,idx.ind.2)
pairs.select = c(all.pairs[2,idx.ind.1],all.pairs[1,idx.ind.2])
dat <- data.frame(x=coord[pairs.select,1],y=coord[pairs.select,2],ec.truncT=ec.trunc[idx.ind],
    ec.logskew=ec.logskew[idx.ind],tc.logskew1=tc.logskew1[idx.ind],
    tc.truncT1=tc.truncT1[idx.ind])
dat[nrow(dat)+1,c("x","y")] <- coord[idx.min,]
limits.max <- max(c(ec.trunc,ec.logskew,tc.logskew1,tc.truncT1))

p1 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=ec.truncT)) + scale_fill_gradient(low="green",high="red",limits=c(1,limits.max),name="Empirical") +
        ggtitle("Truncated extremal t processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))
    
p2 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=tc.truncT1)) + scale_fill_gradient(low="green",high="red",limits=c(1,limits.max),name="Method 1") +
        ggtitle("Truncated extremal t processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))

p3 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=ec.logskew)) + scale_fill_gradient(low="green",high="red",limits=c(1,limits.max),name="Empirical") +
        ggtitle("Log-skew based max-stable processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))


p4 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=tc.logskew1)) + scale_fill_gradient(low="green",high="red",limits=c(1,limits.max),name="Method 1") +
        ggtitle("Log-skew based max-stable processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))

pdf(paste0("figures/extcoef_map_",idx.min,".pdf"),width=12,height=8)
grid.arrange(p1, p2, p3, p4, ncol=3, nrow=2, widths=c(1,1,1), heights=c(1,1))
dev.off()
#}

## fit the model 
# fit the truncated extremal t model
system.time( fit.truncT <- fit.model(data=Z.trunc,loc=coord,init=c(0.5,1,2),fixed=c(F,F,T),thres=0.9,model="truncT",ncores=10,maxit=500,method="L-BFGS-B",lb=c(0.01,0.01),ub=c(10,1.9),bootstrap=TRUE) )

# fit the log-skew based model
system.time( fit.logskew <- fit.model(data=Z.logskew,loc=coord,init=c(0.5,1,2),fixed=c(F,F,T),thres=0.9,model="logskew",method="L-BFGS-B",lb=c(0.1,0.1,-Inf),ub=c(10,1.9,Inf),bootstrap=TRUE,ncores=10,maxit=10000))

system("say \'your program has finished\'")
#fit.result <- MCLE.BR(data=t(Z.logskew[1:10,1:100]),init=c(0.5,1),fixed=c(F,F),distmat=coord[1:10,],FUN = cov.func,index=combn(10,2),ncores=10,method="Nelder-Mead",maxit=1000,hessian=FALSE)
#idx.pairs <- which(all.pairs[1,]==45 | all.pairs[2,]==45)
cov.mat = cov.func(coord,fit.logskew$par[1:2])
alpha = alpha.func(coord,2)

fitted.extcoef.logskew1.1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=list(alpha=alpha,sigma=cov.mat),model="logskew1"),mc.cores=10)
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=fitted.extcoef.logskew1.1,type="p",cex=0.5,col="red",pch=20)
jitter(ec.logskew - fitted.extcoef.logskew1.1)
mean((fitted.extcoef.logskew1.1[idx.pairs] - ec.logskew[idx.pairs] )^2)

cov.mat = cov.func(coord,fit.logskew$par[1:2])
alpha = alpha.func(coord,fit.logskew$par[-c(1:2)])
fitted.extcoef.logskew1.2 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=list(alpha=alpha,sigma=cov.mat),model="logskew1"),mc.cores=10)
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=fitted.extcoef.logskew1.2,type="p",cex=0.5,col="red",pch=20)
boxplot(ec.logskew - fitted.extcoef.logskew1.2)
mean((fitted.extcoef.logskew1.2[idx.pairs] - ec.logskew[idx.pairs] )^2)

plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="black")
points(x=diff.mat[t(all.pairs)],y=fitted.extcoef.logskew1,type="p",cex=0.5,col="red",pch=20)

save.image("data/simulation_test.RData")
