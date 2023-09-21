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
    logphi = log(mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[1])
    gamma_1 = log(gamma((nu+1)/2))
    a_fun <- function(j,upper){
        sigma.11_j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F]) 
        val = mvtnorm::pmvt(lower=rep(0,n-1),upper=upper,delta=sigma[-j,j],sigma=sigma.11_j/(nu+1),df=nu+1)[1]
        return(list(log(val),sigma.11_j))
    }
    if(parallel) T_j <- mclapply(1:n,FUN=a_fun,upper=rep(Inf,n-1),mc.cores=ncores) else T_j <- lapply(1:n,FUN=a_fun,upper=rep(Inf,n-1))
    T_j_val = sapply(T_j,function(x) x[[1]])
    sigma.list = lapply(T_j,function(x) x[[2]])
    a_j = T_j_val-logphi+log(2)*((nu-2)/2)+gamma_1-1/2*log(pi)
    a = exp(a_j)
    # Simulate the max-stable process
    func <- function(idx){
        r = rexp(1)
        z = rep(1,n)
        z[-1] = tmvtsim(n=1,k=n-1, means=sigma[-1,1], sigma=sigma.list[[1]]/(nu+1), lower=rep(0,n-1),upper=rep(Inf,n-1),df=nu+1)[[1]]
        z[-1] = z[-1]^nu * a[1] / a[-1]
        z = z/r
        for(j in 2:n){
            r = rexp(1)
            while(1/r > z[j]){
                z_temp = rep(1,n)
                z_temp[-j] = tmvtsim(n=1,k=n-1, means=sigma[-j,j], sigma=sigma.list[[j]]/(nu+1), lower=rep(0,n-1),upper=rep(Inf,n-1),df=nu+1)[[1]]
                z_temp[-j] = z_temp[-j]^nu * a[j] / a[-j]
                z_temp = z_temp / r
                if(!any( z_temp[0:(j-1)] > z[0:(j-1)])){
                    z = pmax(z,z_temp)
                }
                r = r + rexp(1)
            }
        }
        return(z)
    }
    if(parallel) Z = mclapply(1:m,func,mc.cores=ncores) else Z = lapply(1:m,func)
    Z = matrix(unlist(Z),byrow=TRUE, nrow=m)
    return(Z)
}

# Simulate a log-skew normal based max-stable process
# @m : number of replicates
# @par : parameters of the model
    # @par[[1]] : nu: degrees of freedom
    # @par[[2]] : sigma: correlation matrix
simu_logskew <- function(m,par,parallel=FALSE,ncores=NULL){    
    alpha = par[[1]];sigma = par[[2]]
    n = nrow(sigma)
    omega = diag(sqrt(diag(sigma)))
    omega.inv = diag(diag(omega)^(-1))
    sigma.bar = omega.inv %*% sigma %*% omega.inv
    chol.sigma = chol(sigma)
    inv.sigma = chol2inv(chol.sigma)
    logdet.sigma = sum(log(diag(chol.sigma)))*2
    delta = c(sigma.bar %*% alpha)/c(1+ t(alpha) %*% sigma.bar %*% alpha)
    a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)
    sum.inv.sigma = sum(inv.sigma)
    ones = rep(1,n)
    
    b = c(t(alpha) %*% omega.inv %*% ones)
    alpha.hat = c(t(alpha) %*% omega.inv %*% omega * (1 + b^2)^(-1/2))
    alpha.0.hat = (1+b^2)^(-1/2) * b * (sum.inv.sigma)^(-1/2)
    tau.hat = alpha.0.hat * (1 + t(alpha.hat) %*% sigma.bar %*% alpha.hat)^(-1/2)

    func_paras <- function(idx){
        sigma.11 = sigma[idx,idx]
        mu.2.1 = sigma[-idx,idx]/sigma.11 * a[idx]
        sigma.2.1 = sigma[-idx,-idx] - sigma[-idx,idx,drop=F] %*% sigma[idx,-idx,drop=F]/sigma.11
        omega.2.1 = diag(sqrt(diag(sigma.2.1)))
        omega.2.1.inv = diag(diag(omega.2.1)^(-1))
        sigma.2.1.bar.true = omega.2.1.inv %*% sigma.2.1 %*% omega.2.1.inv
        alpha.2.1 = omega.2.1 %*% omega.inv[-idx,-idx] %*% alpha.hat[-idx]
        sigma.2.1.bar = sigma.bar[-idx,-idx,drop=F] - sigma.bar[-idx,idx,drop=F] %*% sigma.bar[idx,-idx,drop=F]
        alpha.1.2 = c((1 + alpha.hat[-idx] %*% sigma.2.1.bar %*% alpha.hat[-idx])^(-1/2) * (alpha.hat[idx] + sigma.bar[idx,-idx] %*% alpha.hat[-idx]))
        tau.2.1 = c(tau.hat * (1 + alpha.1.2^2)^(1/2) + alpha.1.2 * omega.inv[idx,idx]*a[idx])
        
        delta.2.1 = c(c(1 + t(alpha.2.1) %*% sigma.2.1.bar.true %*% alpha.2.1)^(-1/2) * sigma.2.1.bar.true %*% alpha.2.1)
        Psi = sigma.2.1.bar.true - delta.2.1 %*% t(delta.2.1)
        Psi.chol = chol(Psi)
        return(list(psi.chol = Psi.chol,delta = delta.2.1,omega = omega.2.1,tau = tau.2.1,mu = mu.2.1))
    }
    
    if(parallel){
        paras.list = mclapply(1:n,func_paras,mc.cores=ncores)
        
    }else{
       paras.list = lapply(1:n,func_paras)
    }

    # Simulate the max-stable process
    func <- function(idx){
        r = rexp(1)
        z = rep(1,n)
        u0 = rnorm(1)
        while(u0 <= paras.list[[1]]$tau){
            u0 = rnorm(1)
        }
        u1 = paras.list[[1]]$omega %*% (paras.list[[1]]$psi.chol %*% rnorm(n-1) + paras.list[[1]]$delta * u0) + paras.list[[1]]$mu
        z[-1] = exp(u1-a[-1])
        z = z / r
        for(j in 2:n){
            r = rexp(1)
            while(1/r > z[j]){
                z_temp = rep(1,n)
                u0 = rnorm(1)
                while(u0 <= paras.list[[j]]$tau){
                    u0 = rnorm(1)
                }
                u1 = paras.list[[j]]$omega %*% (paras.list[[j]]$psi.chol %*% rnorm(n-1) + paras.list[[j]]$delta * u0) + paras.list[[j]]$mu
                z_temp[-j] = exp(u1-a[-j])
                z_temp = z_temp / r
                if(!any(z_temp[0:(j-1)] > z[0:(j-1)])){
                    z = pmax(z,z_temp)
                }
                r = r + rexp(1)
            }
        }
        return(z)
    }
    if(parallel)  Z = mclapply(1:m,func,mc.cores=ncores) else Z = lapply(1:m,func)
    Z = matrix(unlist(Z),byrow=TRUE, nrow=m)
    return(Z)
}





