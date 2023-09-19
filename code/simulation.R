###################################################################################
######## This file contains functions to simulate the max-stable processes ########  
###################################################################################

# Simulate a truncated extremal-t max-stable process
# @m : number of replicates
# @par : parameters of the model
    # @par[[1]] : nu: degrees of freedom
    # @par[[2]] : sigma: correlation matrix
simu_truncT <- function(m,par,parallel=FALSE,ncores=NULL){    
    nu = par[[1]];sigma = par[[2]]
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
    
    # Simulate the max-stable process
    func <- function(idx){
        r = rexp(1)
        z = rep(0,n)
        for(j in 1:n){
            while(1/r > z[j]){
                z_temp = rep(1,n)
                z_temp[-j] = tmvtsim(n=1,k=d-1, means=sigma[-j,j], sigma=sigma.list[[j]], lower=rep(0,d),upper=rep(Inf,d),df=nu+1)[[1]]
                z_temp[-j] = z_temp[-j]^nu * a[j] / a[-j]
                z_temp = z_temp / r
                if(any(! z_temp[0:j] < z[0:j])){
                    z = pmax(z,z_temp)
                }else{
                    r = r + rexp(1)
                }
            }
        }
        return(z)
    }
    if(parallel){
        Z = mclapply(1:m,func,mc.cores=ncores)
        
    }else{
       Z = lapply(1:m,func)
    }
    Z = matrix(unlist(Z),byrow=TRUE, nrow=m)
    return(Z)
}

# Example for full rank with strong dependence
d <- 500
rho <- 0.9
Sigma <- matrix(0, nrow=d, ncol=d)
Sigma <- rho^abs(row(Sigma) - col(Sigma))
D1 <- diag(1,d) # Full rank
set.seed(1203)
t0 <- system.time(ans.1 <- tmvmixnorm::rtmvn(n=1, Mean=rep(0,d), Sigma, lower=rep(0,d),int=rep(0.1,d),upper=rep(Inf,d), burn=0))[3]

t1 <- proc.time()
ans.2 <- mvtnorm::rmvnorm(n=1, mean=rep(0,d), Sigma, method="chol")
t2 <- (proc.time() - t1)[3]

t3 <- system.time(ans.4 <- tmvmixnorm::rtmvt(n=1, nu=10, Mean=rep(0,d), Sigma, lower=rep(0,d),D=diag(1,d),int=rep(0.1,d),upper=rep(Inf,d), burn=0))[3] ## worst performance

t4 <- system.time(ans.5 <- tmvnsim(n=1,k=d, means=rep(0,d), sigma=Sigma, lower=rep(0,d),upper=rep(Inf,d))[[1]])[3] ## certainly this is the best
t4 <- system.time(ans.5 <- tmvtsim(n=1,k=d,df=10, means=rep(0,d), sigma=Sigma, lower=rep(0,d),upper=rep(Inf,d))[[1]])[3] ## certainly this is the best

# Simulate a log-skew normal based max-stable process
# @m : number of replicates
# @par : parameters of the model
    # @par[[1]] : nu: degrees of freedom
    # @par[[2]] : sigma: correlation matrix
simu_logskew <- function(m,par,parallel=FALSE,ncores=NULL){    
    alpha = par[[1]];sigma = par[[2]]
    n = nrow(sigma)
    omega = diag(sqrt(diag(sigma)))
    omega_inv = diag(diag(omega)^(-1))
    sigma_bar = omega_inv %*% sigma %*% omega_inv
    chol.sigma = chol(sigma)
    inv.sigma = chol2inv(chol.sigma)
    logdet.sigma = sum(log(diag(chol.sigma)))*2
    delta = sigma_bar %*% alpha/(1+t(alpha) %*% sigma_bar %*% alpha)
    a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)
    sum.inv.sigma = sum(omega_inv)
    
    

    # Simulate the max-stable process
    func <- function(idx){
        r = rexp(1)
        z = rep(0,n)
        for(j in 1:n){
            while(1/r > z[j]){
                z_temp = rep(1,n)
                z_temp[-j] = tmvtsim(n=1,k=d-1, means=sigma[-j,j], sigma=sigma.list[[j]], lower=rep(0,d),upper=rep(Inf,d),df=nu+1)[[1]]
                z_temp[-j] = z_temp[-j]^nu * a[j] / a[-j]
                z_temp = z_temp / r
                if(any(! z_temp[0:j] < z[0:j])){
                    z = pmax(z,z_temp)
                }else{
                    r = r + rexp(1)
                }
            }
        }
        return(z)
    }
    if(parallel){
        Z = mclapply(1:m,func,mc.cores=ncores)
        
    }else{
       Z = lapply(1:m,func)
    }
    Z = matrix(unlist(Z),byrow=TRUE, nrow=m)
    return(Z)
}