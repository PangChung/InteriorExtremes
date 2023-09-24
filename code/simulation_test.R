library(tmvnsim)
library(parallel)
library(evd)
library(doParallel)
source("code/simulation.R")
### testing the simulator ###
d <- 10
coord = as.matrix(expand.grid(1:d,1:d)/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),
                         as.vector(outer(coord[,2],coord[,2],'-'))) 
corr <- function(x,r=0.5,v=1) exp(- (sum(x^2)/r)^v)                          
cov.mat <- matrix(apply(diff.vector, 1, corr), ncol=nrow(coord)) + diag(1e-6,nrow(coord))       
chol(cov.mat)
nu = 10
par1 <- list(nu=nu,sigma=cov.mat)

# Simulate a truncated extremal-t max-stable process
system.time(Z <- simu_truncT(m=1000,par=par1,ncores=10))
#hist(pgev(Z[,1],1,1,1),50,prob=TRUE∏X

# Simulate a log-skew normal based max-stable process
alpha = rep(0,nrow(coord))
par2 <- list(alpha=alpha,sigma=cov.mat)
system.time(Z <- simu_logskew(m=1000,par=par2,ncores=10))
#hist(pgev(Z[,1],1,1,1),50,prob=TRUE)

