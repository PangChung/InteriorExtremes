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
### testing the simulator ###
d <- 4
coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),
                         as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
cov.mat <- cov.func(coord,c(0.5,1))
chol(cov.mat)
nu = 2
m = 10000

par1 <- list(nu=nu,sigma=cov.mat)
system.time(Z.trunc <- simu_truncT(m=m,par=par1,ncores=10))
which(apply(Z.trunc,1,anyDuplicated)>0)

png("figures/marginal_qqplot_truncT.png",width=d*1200,height=d*1200,res=300)
par(mfrow=c(d,d),mgp=c(2,1,0),mar=c(2,2,3,1))
y=(1:m)/(1+m)
for(idx in 1:nrow(Z.trunc)){
x = pgev(Z.trunc[idx,],1,1,1);
p.value = ks.test(x,y)$p.value
qqplot(x,y,main=paste0("Station ",idx, " P value ", round(p.value,2)),xlab="Data",ylab="Theoretical")
abline(0,1,col="red")
}
dev.off()

image(1:10,1:10,z=matrix(log(Z.trunc[,1]),nrow=10),col=rev(heat.colors(10)))

# Simulate a log-skew normal based max-stable process
alpha = alpha.func(coord,c(0,-10,10))
#alpha = 1 + 1.5*coord[,2] - exp(2*sin(10*coord[,2]))
#alpha = (alpha - min(alpha))/(max(alpha)-min(alpha))*10-5
par2 <- list(alpha=alpha,sigma=cov.mat)
system.time(Z.logskew <- simu_logskew(m=10000,par=par2,ncores=10))

png("figures/marginal_qqplot_logskew.png",width=d*1200,height=d*1200,res=300)
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
tc.truncT2 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par1,model="truncT2"),mc.cores=10)

ec.logskew <- apply(all.pairs,2,empirical_extcoef,data=Z.logskew)
tc.logskew1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par2,model="logskew1"),mc.cores=10)
tc.logskew2 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par2,model="logskew2"),mc.cores=10)

idx.pairs = which(all.pairs[1,]==1 | all.pairs[2,]==1)

pdf("figures/extcoef_truncT.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)][idx.pairs],y=ec.trunc[idx.pairs],type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main="Truncated extremal t processes",pch=20,col="#00000033")
points(x=diff.mat[t(all.pairs)][idx.pairs],y=tc.truncT1[idx.pairs],type="p",cex=0.5,col="#ff000033",pch=20)
points(x=diff.mat[t(all.pairs)][idx.pairs],y=tc.truncT2[idx.pairs],type="p",cex=0.5,col="#7eb3d833",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Method 1","Method 2"),col=c("#00000033","#ff000033","#7eb3d833"),
    bty="n",pch=20,cex=1)
dev.off()

pdf("figures/extcoef_logskew.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)][idx.pairs],y=ec.logskew[idx.pairs],type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="#00000033")
points(x=diff.mat[t(all.pairs)][idx.pairs],y=tc.logskew1[idx.pairs],type="p",cex=0.5,col="#ff000033",pch=20)
points(x=diff.mat[t(all.pairs)][idx.pairs],y=tc.logskew2[idx.pairs],type="p",cex=0.5,col="#7eb3d833",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Method 1","Method 2"),col=c("#00000033","#ff000033","#7eb3d833"),
    bty="n", pch=20,cex=1)
dev.off()

## plot the bivariate extremal coefficients map with respect to the central point ##
#for(idx.min in 1:100){
#idx.min = which.min(apply(diff.mat,2,sum))
idx.min = 55
idx.ind.1 = which(all.pairs[1,]==idx.min)
idx.ind.2 = which(all.pairs[2,]==idx.min)
idx.ind = c(idx.ind.1,idx.ind.2)
pairs.select = c(all.pairs[2,idx.ind.1],all.pairs[1,idx.ind.2])
dat <- data.frame(x=coord[pairs.select,1],y=coord[pairs.select,2],ec.truncT=ec.trunc[idx.ind],
    ec.logskew=ec.logskew[idx.ind],tc.logskew1=tc.logskew1[idx.ind],tc.logskew2=tc.logskew2[idx.ind],
    tc.truncT1=tc.truncT1[idx.ind],tc.truncT2=tc.truncT2[idx.ind])
dat[nrow(dat)+1,c("x","y")] <- coord[idx.min,]
limits.max <- max(c(ec.trunc,ec.logskew,tc.logskew1,tc.logskew2,tc.truncT1,tc.truncT2))
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

p3 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=tc.truncT2)) + scale_fill_gradient(low="green",high="red",limits=c(1,limits.max),name="Method 2") +
        ggtitle("Truncated extremal t processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))

p4 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=ec.logskew)) + scale_fill_gradient(low="green",high="red",limits=c(1,limits.max),name="Empirical") +
        ggtitle("Log-skew based max-stable processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))


p5 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=tc.logskew1)) + scale_fill_gradient(low="green",high="red",limits=c(1,limits.max),name="Method 1") +
        ggtitle("Log-skew based max-stable processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))

p6 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=tc.logskew2)) + scale_fill_gradient(low="green",high="red",limits=c(1,limits.max),name="Method 2") +
        ggtitle("Log-skew based max-stable processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))

pdf(paste0("figures/extcoef_map_",idx.min,".pdf"),width=12,height=8)
grid.arrange(p1, p2, p3, p4, p5, p6, ncol=3, nrow=2, widths=c(1,1,1), heights=c(1,1))
dev.off()
#}

## fit the model 
# fit the truncated extremal t model
system.time( fit.truncT <- fit.model(data=Z.trunc,loc=coord,init=c(0.5,1,2),fixed=c(F,F,T),thres=0.8,model="truncT",ncores=10,maxit=500) )
# fit the log-skew based model
system.time( fit.logskew <- fit.model(data=Z.logskew,loc=coord,init=c(0.5,1,0,-10,10),fixed=c(F,F,T,T,F),thres=0.80,model="logskew",ncores=10,maxit=10000) )
#fit.result <- MCLE.BR(data=t(Z.logskew[1:10,1:100]),init=c(0.5,1),fixed=c(F,F),distmat=coord[1:10,],FUN = cov.func,index=combn(10,2),ncores=10,method="Nelder-Mead",maxit=1000,hessian=FALSE)

cov.mat = cov.func(coord,fit.logskew$par[1:2])
alpha = alpha.func(coord,fit.logskew$par[3:5])
fitted.extcoef.logskew1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=list(alpha=alpha,sigma=cov.mat),model="logskew1"),mc.cores=10)
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="#00000033")
points(x=diff.mat[t(all.pairs)],y=fitted.extcoef.logskew1,type="p",cex=0.5,col="#ff000033",pch=20)

