###################################################################################
######## This file contains functions to simulate the max-stable processes ########  
###################################################################################

# Simulate a truncated extremal-t max-stable process
# @n : number of observations
# @par : parameters of the model
simu_truncT <- function(m,par,parallel=TRUE,ncores){    
    nu = par[[1]]
    n = nrow(sigma)
    phi = log(mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[1])
    gamma_1 = log(gamma((nu+1)/2))
    a_fun <- function(j,upper){
        sigma.11_j = (sigma[-j,-j] - sigma[-j,j] %*% sigma[j,-j]) 
        val = mvtnorm::pmvt(lower=rep(0,n-1),upper=upper,mean=sigma[-j,j]*sqrt(nu+1),sigma=sigma.11_j,df=nu+1)[1]
        return(list(log(val)),sigma11_j)
    }
    T_j = lapply(1:n,a_fun,upper=rep(Inf,n-1))
    T_j_val = sapply(T_j,function(x) x[[1]])
    sigma.list = lapply(T_j,function(x) x[[2]])
    a_j = T_j_val-log(phi)+log(2)*((nu-2)/2)+gamma_1-1/2*log(pi)
    
    Z = matrix(0,nrow=n,ncol=m)
    r = rexp(m)
    
    
    

}

# Example for full rank with strong dependence
d <- 500
rho <- 0.9
Sigma <- matrix(0, nrow=d, ncol=d)
Sigma <- rho^abs(row(Sigma) - col(Sigma))
D1 <- diag(1,d) # Full rank
set.seed(1203)
t0 <- system.time(ans.1 <- tmvmixnorm::rtmvn(n=1, Mean=rep(0,d), Sigma, lower=rep(0,d),int=rep(0.1,d),upper=rep(Inf,d), burn=50))[3]

t1 <- proc.time()
ans.2 <- mvtnorm::rmvnorm(n=1, mean=rep(0,d), Sigma, method="chol")
t2 <- (proc.time() - t1)[3]


system.time(ans.3 <- mvnfast::rmvt(n=10000, mu=rep(0,d),sigma=Sigma, df=10, ncores=1))
any(!apply(ans.3,1,function(x){any(x<0)}))

t0 <- system.time(ans.1 <- tmvmixnorm::rtmvt(n=1, nu=10, Mean=rep(0,d), Sigma, lower=rep(0,d),D=diag(1,d),int=rep(0.1,d),upper=rep(Inf,d), burn=0))[3]

# Simulate a log-skew normal based max-stable process
# @n : number of observations
# @par : parameters of the model
simu_logskew <- function(n,par,parallel=TRUE,ncores){  
    alpha = par[[1]]; sigma = par[[2]]
    n = nrow(sigma)
    
}
