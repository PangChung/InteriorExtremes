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
nu = 10
par1 <- list(nu=nu,sigma=cov.mat)

# Simulate a truncated extremal-t max-stable process
system.time(Z.trunc <- simu_truncT(m=10000,par=par1,ncores=10))
hist(pgev(Z.trunc[,1],1,1,1),50,prob=TRUE)
image(1:10,1:10,z=matrix(log(Z.trunc[2,]),nrow=10),col=rev(heat.colors(10)) )

# Simulate a log-skew normal based max-stable process
alpha = rep(0.5,nrow(coord))
par2 <- list(alpha=alpha,sigma=cov.mat)
system.time(Z.logskew <- simu_logskew(m=10000,par=par2,ncores=10))
hist(pgev(Z.logskew[,1],1,1,1),50,prob=TRUE)
image(1:10,1:10,z=matrix(log(Z.logskew[2,]),nrow=10),col=rev(heat.colors(10)) )

# calculate empirical extremal coefficients
empirical_extcoef <- function(p,data){
    return(min(2,max(1,1/mean(1/pmax(data[,p[1]],data[,p[2]])))))
}
all.pairs <- combn(1:ncol(Z.trunc),2)
ec.trunc <- apply(all.pairs,2,empirical_extcoef,data=Z.trunc)
plot(x=diff.mat[t(all.pairs)],y=ec.trunc,type="p",cex=0.5)
ec.logskew <- apply(all.pairs,2,empirical_extcoef,data=Z.logskew)
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5)



