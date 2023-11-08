###################################################################################
######## This file contains functions to simulate the max-stable processes ########  
###################################################################################

# Simulate a truncated extremal-t max-stable process
# @m : number of replicates
# @par : parameters of the model
    # @par[[1]] : nu: degrees of freedom
    # @par[[2]] : sigma: correlation matrix
simu_truncT <- function(m,par,ncores=NULL){    
    nu = par[[2]];sigma = par[[1]]
    n = nrow(sigma)
    phi = mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[[1]]
    gamma_1 = gamma((nu+1)/2)
    a_fun <- function(j,upper){
        sigma.1.j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F]) 
        sigma.j = (sigma - sigma[,j,drop=F] %*% sigma[j,,drop=F] + diag(1e-5,n))/(nu+1) #causing numerical problems
        val = mvtnorm::pmvt(lower=rep(0,n-1)-sigma[-j,j],upper=upper-sigma[-j,j],sigma=sigma.1.j/(nu+1),df=nu+1)[[1]]
        return(list(val,sigma.j))
    }
    if(!is.null(ncores)) T_j <- mclapply(1:n,FUN=a_fun,upper=rep(Inf,n-1),mc.cores=ncores) else T_j <- lapply(1:n,FUN=a_fun,upper=rep(Inf,n-1))
    T_j_val = sapply(T_j,function(x) x[[1]])
    sigma.list = lapply(T_j,function(x) x[[2]])
    sigma.chol = lapply(sigma.list,chol)
    a = T_j_val/phi*2^((nu-2)/2)*gamma_1*(pi^(-1/2))
    # Simulate the max-stable process
    simu <- function(idx){
        r = rexp(1)
        r.hat = 1/r
        z = rep(NA,n)
        denom <- sqrt(rchisq(1, df=nu+1)/(nu+1))
        z = tmvnsim(n=1,k=n, means=rep(0,n),sigma=sigma.list[[1]],
                  lower=(rep(0,n)-sigma[,1])*denom,upper=rep(Inf,n))[[1]]
        z = z/denom + sigma[,1]          
        z = z^nu * a[1]/a
        z = z * r.hat
        for(j in 2:n){
            r = rexp(1)
            r.hat = 1/r
            while(r.hat > z[j]){
                z_temp = rep(NA,n)
                denom <- sqrt(rchisq(1, df=nu+1)/(nu+1))
                z_temp = tmvnsim(n=1,k=n, means=rep(0,n),sigma=sigma.list[[j]],
                  lower=(rep(0,n)-sigma[,j])*denom,upper=rep(Inf,n))[[1]]
                z_temp = z_temp/denom + sigma[,j]          
                z_temp = z_temp^nu*a[j]/a
                z_temp = z_temp * r.hat
                if(!any(z_temp[1:(j-1)] > z[1:(j-1)])){
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
#     @par[[1]] : sigma: correlation matrix
#     @par[[2]] : delta: slant parameter of the skew normal distribution
simu_logskew <- function(m,par,ncores=NULL){  
    delta = par[[2]];sigma = par[[1]]
    n = nrow(sigma)
    omega = diag(sqrt(diag(sigma)))
    omega.inv = diag(diag(omega)^(-1))
    sigma.bar = omega.inv %*% sigma %*% omega.inv
    a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)
    tau.new = delta * diag(omega)
    sigma.star = rbind(cbind(sigma.bar, delta), c(delta, 1))
    sigma.star.chol = chol(sigma.star)
    #message("normalizing constants computed")
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
                while(x[n+1] + tau.new[j] <= 0){ x <- c(t(sigma.star.chol) %*% rnorm(n+1))}
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
    Z = matrix(unlist(Z),byrow=TRUE,nrow=m)
    return(Z)
}

# Simulate a truncated extremal-t max-stable process
bi.simu <- function(m,par,ncores=10,model="truncT"){
    if(model == "truncT"){
        par.func <- function(idx){
            new.par <- list(par[[1]][c(1,idx),c(1,idx)],par[[2]])
        }
        new.par.list <- lapply(2:nrow(par[[1]]),par.func)
        val.list<- mclapply(new.par.list,simu_truncT,m=m,ncores=NULL,mc.cores=10,mc.preschedule = TRUE)
    }
    if(model == "logskew"){
        par.func <- function(idx){
            new.par <- list(par[[1]][c(1,idx),c(1,idx)],par[[2]][c(1,idx)])
        }
        new.par.list <- lapply(2:nrow(par[[1]]),par.func)
        val.list<- mclapply(new.par.list,simu_logskew,m=m,ncores=NULL,mc.cores=10,mc.preschedule = TRUE) 
    }
    return(list(val=val.list,par=new.par.list))
}


ks.test.new <- function(x,n=1000){
    y = runif(length(x))
    ks.test.boot <- function(x,y){
        func <- function(idx){
            ind.sample.x <- sample(1:length(x),n,replace=FALSE)
            ind.sample.y <- sample(1:length(y),n,replace=FALSE)
            val = ks.test(x[ind.sample.x],y[ind.sample.y])$statistic
            return(val)
        }
        val = unlist(lapply(1:1000,func))
        return(val)
    }
    stat.cut.thres <- min(quantile(ks.test.boot(x,x),0.95),quantile(ks.test.boot(y,y),0.95))
    val = median(ks.test.boot(x,y))
    return(val < stat.cut.thres)
}

