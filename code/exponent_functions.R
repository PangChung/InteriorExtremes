#########################################################
###### Intensity function for truncated extremal-t ######
#########################################################

## this function returns the intensity function of the
## truncated extremal-t max-stable processes
intensty_truncT <- function(x,par,parallel=TRUE,ncores=2,log=TRUE){
    sigma = par[[1]];nu = par[[2]]
    n = nrow(sigma)
    chol.sigma = chol(sigma)
    inv.sigma = chol2inv(chol.sigma)
    logdet.sigma = sum(log(diag(chol.sigma)))*2
    phi = log(mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[1])
    gamma_1 = log(gamma((nu+1)/2))
    gamma_n = log(gamma((nu+d)/2))
    a_fun <- function(j,upper){
        sigma.11_j = (sigma[-j,-j] - sigma[-j,j] %*% sigma[j,-j]) * (nu + 1)
        val = mvtnorm::pmvt(lower=rep(0,n-1),upper=upper,mean=sigma[-j,j],sigma=sigma.11_j,df=nu+1)[1]
        return(log(val))
    }
    T_j = unlist(lapply(1:n,a_fun,upper=rep(Inf,n-1)))
    a_j = T_j-log(phi)+log(2)*((nu-2)/2)+gamma_1-1/2*log(pi)
    if(!is.matrix(x)){ x = matrix(x,nrow=n,byrow=FALSE)} 
    func <- function(idx){
        log_x = log(x[,idx])
        x_j = (log_x + a_j) * 1/nu
        x_circ = exp(x_j)
        val = -log(phi) + 1/nu*sum(x_j) + log(nu) +  log(t(x_circ) %*% inv.sigma %*% x_circ) * 
            (-nu-d)/2 + gamma_n
        return(val)
    }
    if(parallel){
        val = unlist(parallel::mclapply(1:ncol(x),func,mc.cores = ncores))
    }else{
        val = unlist(lapply(1:ncol(x),func))
    }
    if(log) return(val)
    else return(exp(val))
} 

## this function returns the exponent function of the
## truncated extremal-t max-stable processes
V_truncT <- function(x,par,parallel=TRUE,ncores=2){
    sigma = par[[1]];nu = par[[2]]
    n = nrow(sigma)
    phi = mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)
    gamma_1 = gamma((nu+1)/2)

    sigma_fun <- function(j){
        sigma_j = (sigma[-j,-j] - sigma[-j,j] %*% sigma[j,-j]) * (nu + 1)
    }

    a_fun <- function(j,upper,sigma_j){
        val = mvtnorm::pmvt(lower=rep(0,n-1),upper=upper,delta=sigma[-j,j],sigma=sigma_j)
        return(val)
    }
    sigma_j = lapply(1:n,sigma_fun)
    T_j = mapply(FUN=a_fun,sigma_j=sigma_j,j=1:n,MoreArgs = list(upper=rep(Inf,n-1)),SIMPLIFY = TRUE)
    a_j = T_j/phi*2^((nu-2)/2)*gamma_1*pi^{-1/2}
    if(!is.matrix(x)){x <- matrix(x,nrow=n,byrow=FALSE)}
    func <- function(idx){
        x_j = x[,idx] * a_j
        x_upper = lapply(1:n,function(i){(x_j[i]/x_j)^{1/nu}})
        V_j = mapply(a_fun,j=1:n,upper=x_upper,sigma_j = sigma_j,SIMPLIFY = TRUE)
        val = sum(V_j/T_j/x[,idx])
    }
    if(parallel){
        val = unlist(parallel::mclapply(1:n,func,mc.cores = ncores))
    }
    else{
        val = unlist(lapply(1:n,func))
    }
    return(val)
}


## this function returns the partial derivatives of the exponent function
## for the truncated extremal-t max-stable processes
partialV_truncT <- function(x,idx,par,parallel=TRUE,ncores=2,log=TRUE){
    sigma = par[[1]];nu = par[[2]]
    n = nrow(sigma)
    phi = mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)
    chol.sigma.11 = chol(sigma[idx,idx])
    inv.sigma.11 = chol2inv(chol.sigma)
    logdet.sigma.11 = sum(log(diag(chol.sigma)))*2
    sigma_T = sigma[-idx,-idx] - sigma[-idx,idx] %*% inv.sigma %*% sigma[-idx,idx] 
    k = length(idx)
    gamma_1 = gamma((nu+1)/2)
    gamma_k = gamma((k+nu)/2)
    a_fun <- function(j,upper){
        sigma_j = (sigma[-j,-j] - sigma[-j,j] %*% sigma[j,-j]) * (nu + 1)
        val = mvtnorm::pmvt(lower=rep(0,n-1),upper=upper,delta=sigma[-j,j],sigma=sigma_j)
        return(val)
    }
    T_j = unlist(lapply(1:n,FUN=a_fun,upper=rep(Inf,n-1))
    a_j = T_j/phi*2^((nu-2)/2)*gamma_1*pi^{-1/2}

    if(!is.matrix(x)){x <- matrix(x,nrow=n,byrow=FALSE)}
    
    func <- function(idx_j){
        x_j = x[,idx_j] * a_j
        x_log = log(x_j)
        Q_sigma = (t(x_j) %*% inv.sigma %*% x_j)^(1/2)
        val = 1/nu*sum(x_log[idx]-x[idx,idx_j]) - phi + (nu-2)/2 * log(2) + (1-k)*log(nu) -
            k/2*log(pi) - 1/2*logdet.sigma -(k+nu)*log(Q_sigma)+log(gamma_k)
        upper = (x_j[-idx])^(1/nu)*sqrt(k+nu)/Q_sigma
        loc = sigma[-idx,idx] %*% inv.sigma.11 %*% x_j[idx]^(1/nu)*sqrt(k+nu)/Q_sigma
        val = val + log(mvtnorm::pmvt(lower=rep(0,n-k),upper=upper,delta=loc,sigma=sigma_T,df=k+nu))
        return(val)
    }
    if(parallel){
        val = unlist(parallel::mclapply(1:n,func,mc.cores = ncores))
    }
    else{
        val = unlist(lapply(1:n,func))
    }
    if(log){
        return(val)
    }
    else{
        return(exp(val))
    }
}

###############################################################################
###### Intensity function for log skew-normal based max-stable processes ######
###############################################################################

## this function computes the intensity function 
## for the log skew-normal based max-stable processes
intensity_logskew <- function(x,par,parallel=TRUE,ncores=2,log=TRUE){
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
    if(!is.matrix(x)){x <- matrix(x,nrow=n,byrow=FALSE)}
    b = (t(alpha) %*% omega_inv %*% rep(1,n))^2/sum_inv_sigma)
    beta =  t(alpha) %*% omega %*% (diag(rep(1,n)) + rep(1,n) %*% t(colSums(inv.sigma))) * (1+b)^(-1/2)
    A = inv.sigma %*% rep(1,n) %*% t(rep(1,n)) %*% inv.sigma/sum.inv.sigma - inv.sigma
    delta.hat = (1+b)^(-1/2)*sqrt(b)
    func <- function(idx){
        x_log = log(x[,idx])
        x_circ = x_log + a
        val = (-n-3)/2 * log(2) - (d-1)/2*log(pi) + pnorm(t(beta) %*% x_circ + delta.hat/sqrt(sum.inv.sigma),log.p=TRUE) -
            1/2*logdet.sigma - 1/2*log(sum.inv.sigma) - sum(x_log)
        val = val - 1/2 * t(x_circ) %*% A %*% x_circ - (t(rep(1,n)) %*% inv.sigma %*% x_circ - 1)/2/sum.inv.sigma
        return(val)
    }
    if(parallel){
        val = unlist(parallel::mclapply(1:ncol(x),func,mc.cores = ncores))
    }
    else{
        val = unlist(lapply(1:ncol(x),func))
    }
    if(log)
        return(val)
    else
        return(exp(val))    
}

## this function computes the exponent function 
## for the log skew-normal based max-stable processes

V_logskew <- function(x,par,parallel=TRUE,ncores=2){
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
    if(!is.matrix(x)){x <- matrix(x,nrow=n,byrow=FALSE)}
    b = t(alpha) %*% omega_inv %*% sigma
    I.mat1 = diag(rep(1,n))
    I.mat2 = diag(rep(1,n-1))
    func <- function(j){        
        A.j = cbind(I.mat2[,0:(j-1)],rep(-1,n-1),I.mat2[,j:(n-1)])
        sigma.j = t(A.j) %*% sigma A.j
        sigma.j.chol = chol(sigma.j)
        sigma.j.inv = chol2inv(sigma.j.chol)
        omega.j = sqrt(diag(diag(sigma.j)))
        omega.j.inv = diag(diag(omega.j)^(-1))
        sigma.j.bar = omega.j.inv %*% sigma.j %*% omega.j.inv
        alpha.hat = (1 - t(delta) %*% omega %*% A.j %*% sigma.j.inv %*% t(A.j) %*% omega %*% delta)^(-1/2) %*% omega.j %*% 
            sigma.j.inv %*% t(A.j) %*% omega %*% delta   
        u.j = t(A.j) %*% sigma %*% I.mat1[,j]
        b1 = t(alpha.hat) %*% sigma.j.bar %*% alpha.hat
        b2 = t(alpha) %*% omega_inv %*% sigma %*% I.mat1[,j] (1+b.j)^(-1/2)
        b3 = -(1+b1)^(-1/2)*sigma.j.bar %*% alpha.hat)
        sigma_circ = cbind(rbind(sigma.j.bar,b3),rbind(b3,1))
        func_temp <- function(i){
            xi = x[,i]
            mu = rbind(omega.j.inv %*% (a[-j] - a[j] + log(xi[-j]/xi[j])-u.j),b2)
            val = pnorm(mu)^(-1)/xi[j] * mvtnorm::pmvnorm(lower=rep(-Inf,n+1),upper=mu,sigma=sigma_circ)[1]
            return(val)
        }
        val = unlist(lapply(1:ncol(x),func_temp))
        return(val)
    }
    if(parallel){
        val = parallel::mclapply(1:n,func,mc.cores = ncores)
        val = matrix(unlist(val),nrow=n,byrow=TRUE)
        val = apply(val,2,sum)
    }
    else{
        val = lapply(1:n,func)
        val = matrix(unlist(val),nrow=n,byrow=TRUE)
        val = apply(val,2,sum)
    }
    return(val)
}

## this function returns the partial derivatives of the exponent function
## for the truncated extremal-t max-stable processes
partialV_logskew <- function(x,idx,par,parallel=TRUE,ncores=2,log=TRUE){
    alpha = par[[1]];sigma = par[[2]]
    n = nrow(sigma)
    sigma.chol = chol(sigma)
    sigma.inv = chol2inv(sigma.chol)
    sum.sigma.inv = sum(sigma.inv)
    
    sigma.tilde.inv = sigma.inv %*% rep(1,n) %*% t(rep(1,n)) %*% sigma.inv/sum.sigma.inv - sigma.inv
    sigma.tilde.chol = chol(sigma.tilde.inv)
    sigma.tilde = chol2inv(sigma.tilde.chol)
    mu.tilde = sigma.tilde.inv %*% sigma.inv %*% rep(1,n)/sum.sigma.inv/2
    
    alpha.tilde = 
}