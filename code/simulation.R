###################################################################################
######## This file contains functions to simulate the max-stable processes ########  
###################################################################################

# Simulate a truncated extremal-t max-stable process
# @m : number of replicates
# @par : parameters of the model
    # @par[[1]] : nu: degrees of freedom
    # @par[[2]] : sigma: correlation matrix
simu_truncT <- function(m,par,ncores=NULL){    
    message("start computing matrices")
    nu = par[[2]];sigma = par[[1]]
    n = nrow(sigma)
    phi = mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[[1]]
    gamma_1 = gamma((nu+1)/2)
    a_fun <- function(j,upper){
        sigma.1.j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F]/sigma[j,j]) 
        sigma.j = (sigma - sigma[,j,drop=F] %*% sigma[j,,drop=F]/sigma[j,j])/(nu+1)/sigma[j,j] 
        val = mvtnorm::pmvt(lower=rep(0,n-1)-sigma[-j,j]/sigma[j,j],upper=upper-sigma[-j,j]/sigma[j,j],sigma=sigma.1.j/(nu+1),df=nu+1)[[1]]
        return(list(val,sigma.j))
    }
    if(!is.null(ncores)) T_j <- mclapply(1:n,FUN=a_fun,upper=rep(Inf,n-1),mc.cores=ncores) else T_j <- lapply(1:n,FUN=a_fun,upper=rep(Inf,n-1))
    T_j_val = sapply(T_j,function(x) x[[1]])
    sigma.list = lapply(T_j,function(x) x[[2]])
    sigma.chol = lapply(sigma.list,chol)
    a = T_j_val/phi*2^((nu-2)/2)*gamma_1*(pi^(-1/2))
    
    # Simulate the max-stable process
    idx.loc = cbind(1:m,1)
    idx.loc.sum = lapply(1:n,function(id){which(idx.loc[,2]==id)})
    r = rexp(m)
    r.hat = 1/r

    func <- function(idx.j,j){
        m.idx.j = length(idx.j)
        if(m.idx.j > 0 & j<=n){
            val = TruncatedNormal::rtmvt(n=m.idx.j,mu=sigma[,j]/sigma[j,j],sigma=sigma.list[[j]],df=nu+1,lb=rep(0,n),ub=rep(Inf,n))
            if(!is.matrix(val)){ val = matrix(val,nrow=m.idx.j)}
            val = t(t(val)^nu * a[j]/a)
            return(val)
        }
        return(NULL)
    }

    func.compare <- function(idx,j){
        val = !any(z[idx,1:(j-1)] < z.temp[idx,1:(j-1)]) 
        if(val){
            return(idx)
        }else{
            return(NULL)
        }
    }
    
    func.pmax <- function(idx){
        return(pmax(z[idx,],z.temp[idx,]))
    }

    if(!is.null(ncores)){
        z = mcmapply(func,idx.j=idx.loc.sum,j=1:n,mc.cores=ncores,SIMPLIFY=FALSE)
    }else{
        z = mapply(func,idx.j=idx.loc.sum,j=1:n,SIMPLIFY=FALSE)
    }
    
    message("start simulation")
    z = do.call(rbind,z) 
    z = z *  r.hat
    idx.loc[,2] = 2
    idx.loc.sum = lapply(1:n,function(id){which(idx.loc[,2]==id)})
    idx.new.r <- idx.pass.r <-  rep(TRUE,m)
    idx.finish <- idx.r.old <- rep(FALSE,m)
    count = 1
    while(any(!idx.finish)){
        r[idx.r.old] = r[idx.r.old ] + rexp(sum(idx.r.old))
        r.hat[idx.r.old] = 1/r[idx.r.old]
        if(sum(idx.new.r) != 0){
            r[idx.new.r] = rexp(sum(idx.new.r))
            r.hat[idx.new.r] = 1/r[idx.new.r]
        }
        idx.pass.r[!idx.finish] <- r.hat[!idx.finish] > z[idx.loc[!idx.finish,,drop=FALSE]]
        idx.pass.r[idx.finish] = FALSE
        idx.go.ahead = !idx.pass.r & !idx.finish
        if(sum(idx.go.ahead) > 0){
            idx.loc[idx.go.ahead,2] = idx.loc[idx.go.ahead,2] +  1
            idx.new.r[idx.go.ahead] = TRUE
        }
        if(sum(idx.pass.r)!=0){
            idx.new.r[idx.pass.r] = FALSE
            idx.loc.sum.simu = lapply(1:n,function(id){which(idx.loc[idx.pass.r,2]==id)})
            if(!is.null(ncores)){
                z.temp.0 = mcmapply(func,idx.j=idx.loc.sum.simu,j=1:n,mc.cores=ncores,SIMPLIFY=FALSE)
            }else{
                z.temp.0 = mapply(func,idx.j=idx.loc.sum.simu,j=1:n,SIMPLIFY=FALSE)
            }
            z.temp = z
            idx.temp = which(idx.pass.r)[unlist(idx.loc.sum.simu)]
            z.temp[idx.temp,] = do.call(rbind,z.temp.0) * r.hat[idx.temp]
            idx.temp = unlist(mapply(func.compare,idx=which(idx.pass.r),j=idx.loc[idx.pass.r,2],SIMPLIFY=FALSE))
            if(length(idx.temp)!=0){
                z[idx.temp,] = do.call(rbind,lapply(idx.temp,func.pmax))
                idx.loc[idx.temp,2] = idx.loc[idx.temp,2] + 1
                idx.new.r[idx.temp] = TRUE
            }
        }
        idx.finish = idx.loc[,2] > n
        idx.new.r = idx.new.r & !idx.finish
        idx.r.old = !idx.new.r & !idx.finish
        count = count + 1
        if(sum(idx.finish) %% 1000 == 0){
            print(paste0(c(count,sum(idx.finish),m),collapse = "/"))
        }
    }
    return(z)
}


# Simulate a log-skew normal based max-stable process
# @m : number of replicates
# @par : parameters of the model
#     @par[[1]] : sigma: correlation matrix
#     @par[[2]] : delta: slant parameter of the skew normal distribution
simu_logskew <- function(m,par,ncores=NULL){  
    delta = par[[2]];sigma = par[[1]]
    n = nrow(sigma)
    a = log(2) + diag(sigma)/2 + sapply(delta,pnorm,log.p=TRUE)
    sigma.star = rbind(cbind(sigma, delta), c(delta, 1))
    sigma.star.chol = chol(sigma.star)
    #message("normalizing constants computed")
    # Simulate the max-stable process
    simu <- function(idx){
        r <- rexp(1)
        r.hat <- 1/r
        x <- c(t(sigma.star.chol) %*% rnorm(n+1))
        while(x[n+1] + delta[1] < 0){ x <- c(t(sigma.star.chol) %*% rnorm(n+1))}
        x0 = x[1:n] + sigma[,1]
        z <- exp(x0-a);z <- z/z[1]
        z <- z * r.hat
        z.list = list()
        z.list[[1]] = z
        count = 1
        for(j in 2:n){
            r <- rexp(1)
            r.hat <- 1/r
            while(r.hat > z[j]){
                x <- c(t(sigma.star.chol) %*% rnorm(n+1))
                while(x[n+1] + delta[j] <= 0){ x <- c(t(sigma.star.chol) %*% rnorm(n+1))}
                x0 <- x[1:n] + sigma[,j]
                z_temp <- exp(x0-a);z_temp <- z_temp/z_temp[j]
                z_temp <- z_temp * r.hat
                if(!any(z_temp[1:(j-1)] > z[1:(j-1)])){
                    z <- pmax(z,z_temp)
                }
                r <- r + rexp(1)
                r.hat <- 1/r
                z.list[[count+1]] <- z_temp
                count = count + 1
            }
        }
        return(z.list)
    }
    if(!is.null(ncores)){ 
        Z = mclapply(1:m,simu,mc.cores=ncores,mc.set.seed = TRUE)
    }else {Z = lapply(1:m,simu)}
    Z = matrix(unlist(Z),byrow=TRUE,ncol=n)
    return(Z)
}

simu_logskew2 <- function(m,par,ncores=NULL){  
    delta = par[[2]];sigma = par[[1]]
    n = nrow(sigma)
    a = log(2) + diag(sigma)/2 + sapply(delta,pnorm,log.p=TRUE)
    sigma.star = rbind(cbind(sigma, delta), c(delta, 1))
    sigma.star.chol = chol(sigma.star)
    # Simulate the max-stable process
    idx.loc = cbind(1:m,1)
    idx.loc.sum = lapply(1:n,function(id){which(idx.loc[,2]==id)})
    r = rexp(m)
    r.hat = 1/r
    func <- function(idx.j,j){
        m.idx.j = length(idx.j)
        if(m.idx.j==0){return(NULL)}
        func.i <- function(i){
            x <- c(t(sigma.star.chol) %*% rnorm(n+1))
            while(x[n+1] + delta[j] <= 0){ x <- c(t(sigma.star.chol) %*% rnorm(n+1))}
            x0 <- c(x[1:n] + sigma[,j])
            val <- exp(x0-a);val <- val/val[j]
            return(val)
        }
        if(!is.null(ncores)){
            z = mclapply(1:m.idx.j,func.i,mc.cores=ncores)    
        }else{
            z = lapply(1:m.idx.j,func.i)
        }
        return(do.call(rbind,z))
    }

    func.compare <- function(idx,j){
        val = !any(z[idx,1:(j-1)] < z.temp[idx,1:(j-1)]) 
        if(val){
            return(idx)
        }else{
            return(NULL)
        }
    }
    
    func.pmax <- function(idx){
        return(pmax(z[idx,],z.temp[idx,]))
    }

    z = mapply(func,idx.j=idx.loc.sum,j=1:n,SIMPLIFY=FALSE)
    z = do.call(rbind,z) 
    z = z *  r.hat
    idx.loc[,2] = 2
    idx.loc.sum = lapply(1:n,function(id){which(idx.loc[,2]==id)})
    idx.new.r <- idx.pass.r <-  rep(TRUE,m)
    idx.finish <- idx.r.old <- rep(FALSE,m)
    count = 1
    while(any(!idx.finish)){
        r[idx.r.old] = r[idx.r.old ] + rexp(sum(idx.r.old))
        r.hat[idx.r.old] = 1/r[idx.r.old]
        if(sum(idx.new.r) != 0){
            r[idx.new.r] = rexp(sum(idx.new.r))
            r.hat[idx.new.r] = 1/r[idx.new.r]
        }
        idx.pass.r[!idx.finish] <- r.hat[!idx.finish] > z[idx.loc[!idx.finish,,drop=FALSE]]
        idx.pass.r[idx.finish] = FALSE
        idx.go.ahead = !idx.pass.r & !idx.finish
        if(sum(idx.go.ahead) > 0){
            idx.loc[idx.go.ahead,2] = idx.loc[idx.go.ahead,2] +  1
            idx.new.r[idx.go.ahead] = TRUE
        }
        if(sum(idx.pass.r)!=0){
            idx.new.r[idx.pass.r] = FALSE
            idx.loc.sum.simu = lapply(1:n,function(id){which(idx.loc[idx.pass.r,2]==id)})
            z.temp.0 = mapply(func,idx.j=idx.loc.sum.simu,j=1:n,SIMPLIFY=FALSE)
            z.temp = z
            idx.temp = which(idx.pass.r)[unlist(idx.loc.sum.simu)]
            z.temp[idx.temp,] = do.call(rbind,z.temp.0) * r.hat[idx.temp]
            idx.temp = unlist(mapply(func.compare,idx=which(idx.pass.r),j=idx.loc[idx.pass.r,2],SIMPLIFY=FALSE))
            if(length(idx.temp)!=0){
                z[idx.temp,] = do.call(rbind,lapply(idx.temp,func.pmax))
                idx.loc[idx.temp,2] = idx.loc[idx.temp,2] + 1
                idx.new.r[idx.temp] = TRUE
            }
        }
        idx.finish = idx.loc[,2] > n
        idx.new.r = idx.new.r & !idx.finish
        idx.r.old = !idx.new.r & !idx.finish
        count = count + 1
        if(sum(idx.finish) %% 1000 == 0){
            print(paste0(c(count,sum(idx.finish),m),collapse = "/"))
        }
    }
    return(z)
}

## simulate the r-Pareto process associated with the log-skew normal process ##
simu_Pareto_logskew <- function(m,par,riskr,ncores=NULL){
    delta = par[[2]];sigma = par[[1]]
    n = nrow(sigma)
    a = log(2) + diag(sigma)/2 + sapply(delta,pnorm,log.p=TRUE)
    sigma.star = rbind(cbind(sigma, delta), c(delta, 1))
    sigma.star.chol = chol(sigma.star)
    Z = matrix(NA,ncol=n,nrow=m)
    # Simulate the r-Pareto process 
    func <- function(m.i){
        func.i <- function(i){
            j = sample(1:n,1)
            x <- c(t(sigma.star.chol) %*% rnorm(n+1))
            while(x[n+1] + delta[j] <= 0){ x <- c(t(sigma.star.chol) %*% rnorm(n+1))}
            x0 <- c(x[1:n] + sigma[,j])
            val <- exp(x0-a);val <- val/val[j]
            return(val/sum(val))
        }
        if(!is.null(ncores)){
            z = mclapply(1:m.i,func.i,mc.cores=ncores)    
        }else{
            z = lapply(1:m.i,func.i)
        }
        return(do.call(rbind,z))
    }
    r = mev::rgp(m,1,1,1)
    z = func(m)*r
    idx.finish <- apply(z,1,riskr) > 1
    Z[idx.finish,] = z[idx.finish,]
    while(any(!idx.finish)){
        m.temp = sum(!idx.finish)
        z.temp = func(m.temp)*mev::rgp(m.temp,1,1,1)
        idx.finish.temp = apply(z.temp,1,riskr) > 1
        if(any(idx.finish.temp)){
            idx.finish[!idx.finish] = idx.finish.temp
            Z[!idx.finish[idx.finish.temp],] = z.temp[idx.finish.temp,]
        }
    }
    return(Z)
}

## simulate the r-Pareto process associated with the truncated process ##
simu_Pareto_truncT <- function(m,par,riskr,ncores=NULL){
    nu = par[[2]];sigma = par[[1]]
    n = nrow(sigma)
    phi = mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[[1]]
    gamma_1 = gamma((nu+1)/2)
    a_fun <- function(j,upper){
        sigma.1.j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F]/sigma[j,j]) 
        sigma.j = (sigma - sigma[,j,drop=F] %*% sigma[j,,drop=F]/sigma[j,j])/(nu+1)/sigma[j,j] 
        val = mvtnorm::pmvt(lower=rep(0,n-1)-sigma[-j,j]/sigma[j,j],upper=upper-sigma[-j,j]/sigma[j,j],sigma=sigma.1.j/(nu+1),df=nu+1)[[1]]
        return(list(val,sigma.j))
    }
    if(!is.null(ncores)) T_j <- mclapply(1:n,FUN=a_fun,upper=rep(Inf,n-1),mc.cores=ncores) else T_j <- lapply(1:n,FUN=a_fun,upper=rep(Inf,n-1))
    T_j_val = sapply(T_j,function(x) x[[1]])
    sigma.list = lapply(T_j,function(x) x[[2]])
    sigma.chol = lapply(sigma.list,chol)
    a = T_j_val/phi*2^((nu-2)/2)*gamma_1*(pi^(-1/2))
    Z = matrix(NA,ncol=n,nrow=m)
    # Simulate the r-Pareto process
    
    func <- function(m.i){
        j = sample(1:n,1)
        if(m.i > 0 & j<=n){
            val = TruncatedNormal::rtmvt(n=m.i,mu=sigma[,j]/sigma[j,j],sigma=sigma.list[[j]],df=nu+1,lb=rep(0,n),ub=rep(Inf,n))
            if(!is.matrix(val)){ val = matrix(val,nrow=m.i)}
            val = t(t(val)^nu * a[j]/a)
            return(val/sum(val))
        }
        return(NULL)
    }
    r = mev::rgp(m,1,1,1)
    z = func(m)*r
    idx.finish <- apply(z,1,riskr) > 1
    Z[idx.finish,] = z[idx.finish,]
    while(any(!idx.finish)){
        m.temp = sum(!idx.finish)
        z.temp = func(m.temp)*mev::rgp(m.temp,1,1,1)
        idx.finish.temp = apply(z.temp,1,riskr) > 1
        if(any(idx.finish.temp)){
            idx.finish[!idx.finish] <- idx.finish.temp
            Z[!idx.finish[idx.finish.temp],] <- z.temp[idx.finish.temp,]
        }
    }
    return(Z)
}

