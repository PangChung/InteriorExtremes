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

simu_logskew2 <- function(m,par,ncores=NULL){  
    delta = par[[2]];sigma = par[[1]]
    n = nrow(sigma)
    omega = diag(sqrt(diag(sigma)))
    omega.inv = diag(diag(omega)^(-1))
    sigma.bar = omega.inv %*% sigma %*% omega.inv
    a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)
    tau.new = delta * diag(omega)
    sigma.star = rbind(cbind(sigma.bar, delta), c(delta, 1))
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
            while(x[n+1] + tau.new[j] <= 0){ x <- c(t(sigma.star.chol) %*% rnorm(n+1))}
            x0 <- c(omega %*% x[1:n] + sigma[,j])
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
    return(Z)
}

# Simulate a truncated extremal-t max-stable process
bi.simu <- function(m,par,ncores=10,model="truncT",alpha.para=TRUE){
    if(model == "truncT"){
        par.func <- function(idx){
            new.par <- list(par[[1]][c(1,idx),c(1,idx)],par[[2]])
        }
        new.par.list <- lapply(2:nrow(par[[1]]),par.func)
        val.list<- mclapply(new.par.list,simu_truncT,m=m,ncores=NULL,mc.cores=10,mc.preschedule = TRUE)
    }
    if(model == "logskew"){
        par.func <- function(idx){
            idx = c(1,idx)
            if(alpha.para)
                new.par <- alpha2delta(list(par[[1]][idx,idx],par[[2]][idx]))
            else 
                new.par <- list(par[[1]][idx,idx],par[[2]][idx])
        }
        new.par.list <- lapply(2:nrow(par[[1]]),par.func)
        val.list<- mclapply(new.par.list,simu_logskew,m=m,ncores=NULL,mc.cores=10,mc.preschedule = TRUE) 
    }
    if(model == "logskew2"){
        par.func <- function(idx){
            if(alpha.para)
                new.par <- alpha2delta(list(par[[1]],par[[2]][c(idx,idx)]))
            else 
                new.par <- list(par[[1]],par[[2]][c(idx,idx)])
        }
        new.par.list <- lapply(1:length(par[[2]]),par.func)
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

