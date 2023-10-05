library(tmvnsim)
library(parallel)
library(evd)
library(gridExtra)
library(ggplot2)
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
alpha = alpha.func(coord,c(0.1,-0.5,0.5))
par2 <- list(alpha=alpha,sigma=cov.mat)
system.time(Z.logskew <- simu_logskew(m=10000,par=par2,ncores=10))
hist(pgev(Z.logskew[,1],1,1,1),50,prob=TRUE)
image(1:10,1:10,z=matrix(log(Z.logskew[2,]),nrow=10),col=rev(heat.colors(10)) )

all.pairs <- combn(1:ncol(Z.trunc),2)
all.pairs.list = split(all.pairs,col(all.pairs))
ec.trunc <- apply(all.pairs,2,empirical_extcoef,data=Z.trunc)

tc.truncT1 <- true_extcoef(all.pairs,par=par1,model="truncT1")
tc.truncT2 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par1,model="truncT2"),mc.cores=10)

ec.logskew <- apply(all.pairs,2,empirical_extcoef,data=Z.logskew)
tc.logskew1 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par2,model="logskew1"),mc.cores=10)
tc.logskew2 <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=par2,model="logskew2"),mc.cores=10)

pdf("figures/extcoef_truncT.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.trunc,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main="Truncated extremal t processes",pch=20,col="#0000001A")
points(x=diff.mat[t(all.pairs)],y=tc.truncT1,type="p",cex=0.5,col="#ff00001A",pch=20)
points(x=diff.mat[t(all.pairs)],y=tc.truncT2,type="p",cex=0.5,col="#7eb3d81A",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Method 1","Method 2"),col=c("#0000001A","#ff00001A","#7eb3d81A"),
    bty="n",pch=20,cex=1)
dev.off()


pdf("figures/extcoef_logskew.pdf",width=6,height=4)
par(mfrow=c(1,1),mar=c(4,4,2,1),cex.main=1,cex.lab=1,mgp=c(2,1,0))
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",
    main = "Log-skew normal based max-stable processes",pch=20,col="#0000001A")
points(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=0.5,col="#ff00001A",pch=20)
points(x=diff.mat[t(all.pairs)],y=tc.logskew2,type="p",cex=0.5,col="#0f83d61A",pch=20)
abline(h=c(1,2),col="grey",lty=2,cex=2)
legend("topleft",legend=c("Empirical","Method 1","Method 2"),col=c("#0000001A","#ff00001A","#7eb3d81A"),
    bty="n", pch=20,cex=1)
dev.off()

## plot the bivariate extremal coefficients map with respect to the central point ##
#idx.min = which.min(apply(diff.mat,2,sum))
for(idx.min in 1:100){
#idx.min = 35
idx.ind.1 = which(all.pairs[1,]==idx.min)
idx.ind.2 = which(all.pairs[2,]==idx.min)
idx.ind = c(idx.ind.1,idx.ind.2)
pairs.select = c(all.pairs[2,idx.ind.1],all.pairs[1,idx.ind.2])
dat <- data.frame(x=coord[pairs.select,1],y=coord[pairs.select,2],ec.truncT=ec.trunc[idx.ind],
    ec.logskew=ec.logskew[idx.ind],tc.logskew1=tc.logskew1[idx.ind],tc.logskew2=tc.logskew2[idx.ind],
    tc.truncT1=tc.truncT1[idx.ind],tc.truncT2=tc.truncT2[idx.ind])
dat[nrow(dat)+1,c("x","y")] <- coord[idx.min,]
p1 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=ec.truncT)) + scale_fill_gradient(low="blue",high="red",limits=c(1,2),name="Empirical") +
        ggtitle("Truncated extremal t processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))
    
p2 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=tc.truncT1)) + scale_fill_gradient(low="blue",high="red",limits=c(1,2),name="Method 1") +
        ggtitle("Truncated extremal t processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))

p3 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=tc.truncT2)) + scale_fill_gradient(low="blue",high="red",limits=c(1,2),name="Method 2") +
        ggtitle("Truncated extremal t processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))

p4 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=ec.logskew)) + scale_fill_gradient(low="blue",high="red",limits=c(1,2),name="Empirical") +
        ggtitle("Log-skew based max-stable processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))


p5 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=tc.logskew1)) + scale_fill_gradient(low="blue",high="red",limits=c(1,2),name="Method 1") +
        ggtitle("Log-skew based max-stable processes") + theme(axis.text = element_text(size=10), 
        axis.title.x = element_text(size=14), 
        axis.title.y = element_text(size=14),
        plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(fill = "transparent", 
                                    color = "transparent", 
                                    linewidth = 0.5))

p6 <- ggplot(dat) + geom_tile(aes(x=x,y=y,fill=tc.logskew2)) + scale_fill_gradient(low="blue",high="red",limits=c(1,2),name="Method 2") +
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
}
