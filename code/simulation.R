###################################################################################
######## This file contains functions to simulate the max-stable processes ########  
###################################################################################

# Simulate a truncated extremal-t max-stable process
# @m : number of replicates
# @par : parameters of the model
    # @par[[1]] : nu: degrees of freedom
    # @par[[2]] : sigma: correlation matrix
simu_truncT <- function(m,par,ncores=NULL){    
    nu = par[[1]];sigma = par[[2]]
    n = nrow(sigma)
    logphi = log(mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[1])
    gamma_1 = log(gamma((nu+1)/2))
    a_fun <- function(j,upper){
        sigma.1.j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F]) 
        val = mvtnorm::pmvt(lower=rep(0,n-1),upper=upper,delta=sigma[-j,j],sigma=sigma.1.j/(nu+1),df=nu+1)[1]
        return(list(log(val),sigma.1.j))
    }
    if(!is.null(ncores)) T_j <- mclapply(1:n,FUN=a_fun,upper=rep(Inf,n-1),mc.cores=ncores) else T_j <- lapply(1:n,FUN=a_fun,upper=rep(Inf,n-1))
    T_j_val = sapply(T_j,function(x) x[[1]])
    sigma.list = lapply(T_j,function(x) x[[2]])
    a_j = T_j_val-logphi+log(2)*((nu-2)/2)+gamma_1-1/2*log(pi)
    a = exp(a_j)
    message("normalizing constants computed")
    # Simulate the max-stable process
    simu <- function(idx){
        r = rexp(1)
        r.hat = 1/r
        z = rep(NA,n)
        z[1] = 1
        k = nrow(sigma.list[[1]])
        z[-1] = tmvtsim(n=1,k=k, means=sigma[,1], sigma=sigma.list[[1]]/(nu+1), lower=rep(0,k),upper=rep(Inf,k),df=nu+1)[[1]]
        z[-1] = z[-1]^nu * a[1]/a[-1]
        z = z * r.hat
        for(j in 2:n){
            r = rexp(1)
            r.hat = 1/r
            while(r.hat > z[j]){
                z_temp = rep(NA,n);z_temp[j] = 1
                z_temp[-j] = tmvtsim(n=1,k=k, means=sigma[,j], sigma=sigma.list[[j]]/(nu+1), lower=rep(0,k),upper=rep(Inf,k),df=nu+1)[[1]]
                z_temp[-j] = z_temp[-j]^nu * a[j]/a[-j]
                z_temp = z_temp * r.hat
                if(!any( z_temp[1:(j-1)] > z[1:(j-1)])){
                    z = pmax(z,z_temp)
                }
                r = r + rexp(1)
                r.hat = 1/r
            }
        }
        return(z)
    }
    if(!is.null(ncores)) Z = mclapply(1:m,simu,mc.cores=ncores) else Z = lapply(1:m,simu)
    Z = matrix(unlist(Z),byrow=TRUE, nrow=m)
    return(Z)
}




# Simulate a log-skew normal based max-stable process
# @m : number of replicates
# @par : parameters of the model
    # @par[[1]] : alpha: slant parameter of the skew normal distribution
    # @par[[2]] : sigma: correlation matrix
simu_logskew <- function(m,par,ncores=NULL){  
    alpha = par[[1]];sigma = par[[2]]
    n = nrow(sigma)
    omega = diag(sqrt(diag(sigma)))
    omega.inv = diag(diag(omega)^(-1))
    sigma.bar = omega.inv %*% sigma %*% omega.inv
    sigma.bar.chol = chol(sigma.bar)
    delta = c(sigma.bar %*% alpha)/c(1+ t(alpha) %*% sigma.bar %*% alpha)
    a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)

    tau.new = delta * diag(omega)

    sigma.star = rbind(cbind(sigma.bar, delta), c(delta, 1))
    sigma.star.chol = chol(sigma.star)
    message("normalizing constants computed")
    # Simulate the max-stable process
    simu <- function(idx){
        r <- rexp(1)
        r.hat <- 1/r
        x <- c(t(sigma.star.chol) %*% rnorm(n+1))
        while(x[n+1] + tau.new[1] < 0){ x <- c(t(sigma.star.chol) %*% rnorm(n+1))}
        x0 = c(omega %*% x[1:n] + sigma[,1])
        z <- exp(x0-a);z <- z/z[1]
        z <- z * r.hat
        for(j in 2:n){
            r <- rexp(1)
            r.hat <- 1/r
            while(r.hat > z[j]){
                x <- c(t(sigma.star.chol) %*% rnorm(n+1))
                while(x[n+1] + tau.new[j] < 0){ x <- c(t(sigma.star.chol) %*% rnorm(n+1))}
                x0 <- c(omega %*% x[1:n] + sigma[,j])
                z_temp <- exp(x0-a);z_temp <- z_temp/z_temp[j]
                z_temp <- z_temp * r.hat
                if(!any(z_temp[1:(j-1)] > z[1:(j-1)])){
                    z <- pmax(z,z_temp)
                }
                r <- r + rexp(1)
                r.hat <- 1/r
            }
        }
        return(z)
    }

    if(!is.null(ncores)){ 
        Z = mclapply(1:m,simu,mc.cores=ncores)
    }else Z = lapply(1:m,simu)
    
    Z = matrix(unlist(Z),byrow=TRUE, nrow=m)
    return(Z)
}





