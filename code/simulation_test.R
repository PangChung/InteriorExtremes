library(tmvnsim)
library(parallel)
library(evd)
source("code/simulation.R")
### testing the simulator ###
d <- 5
rho <- 1
Sigma <- matrix(0, nrow=d, ncol=d)
Sigma <- 0.8^abs(row(Sigma) - col(Sigma))
nu = 10
par <- list(nu=nu,sigma=Sigma)

# Simulate a truncated extremal-t max-stable process
Z = simu_truncT(m=10000,par=par,parallel=TRUE,ncores=10)
hist(pgev(Z[,1],loc=0,scale=1,shape=1),50,prob=TRUE)


# Simulate a log-skew normal based max-stable process
alpha = rep(0.5,d)
par <- list(alpha=alpha,sigma=Sigma)
Z = simu_logskew(m=100000,par=par,parallel=TRUE,ncores=10)
hist(pgev(Z[,4],loc=0,scale=1,shape=1),50,prob=TRUE)



