#########################################################
###### Intensity function for truncated extremal-t ######
#########################################################

## Returns the normalizing constant for the truncated extremal-t max-stable processes ##
a_fun <- function(par,ncores=NULL){
    sigma = par[[1]];nu=par[[2]]
    n = nrow(sigma)
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    logphi = log(mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),sigma=sigma)[[1]])
    if(is.null(ncores)){
    val = sapply(1:n,function(j){
                sigma_j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F]/sigma[j,j])/(nu + 1)/sigma[j,j]
                mvtnorm::pmvt(lower=-sigma[-j,j]/sigma[j,j],upper=rep(Inf,n-1),sigma=sigma_j,df=nu+1)[[1]]})
    } else {
    val = mcmapply(function(j){
                sigma_j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F]/sigma[j,j])/(nu + 1)/sigma[j,j]
                mvtnorm::pmvt(lower=-sigma[-j,j]/sigma[j,j],upper=rep(Inf,n-1),sigma=sigma_j,df=nu+1)[[1]]},j=1:n,mc.cores=ncores,mc.set.seed = FALSE)
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
    return(list(log(val),logphi))
}

## this function returns the intensity function of the
## truncated extremal-t max-stable processes
intensity_truncT <- function(x,par,T_j=NULL,ncores=NULL,log=TRUE){
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    sigma = par[[1]];nu = par[[2]]
    n = ncol(x)
    if(n==1) return(1/(x^2))
    tryCatch(chol.sigma <- chol(sigma),error=function(e) browser())
    inv.sigma = chol2inv(chol.sigma)
    logdet.sigma = sum(log(diag(chol.sigma)))*2

    gamma_1 = log(gamma((nu+1)/2))
    gamma_n = log(gamma((nu+n)/2))

    a_fun <- function(j,upper=rep(Inf,n-1)){
        sigma_j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F]/sigma[j,j])/(nu + 1)/sigma[j,j]
        val = mvtnorm::pmvt(lower=-sigma[-j,j]/sigma[j,j],upper=rep(Inf,n-1),sigma=sigma_j,df=nu+1)[[1]]
        return(log(val))
    }

    if(is.null(T_j)){
        logphi = log(mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),sigma=sigma)[[1]])
        if(!is.null(ncores)) T_j = unlist(mclapply(1:n,a_fun,mc.cores=ncores,mc.set.seed = FALSE,mc.preschedule = TRUE)) else T_j = unlist(lapply(1:n,a_fun))
    }else{
        T_j = T_j[[1]];logphi = T_j[[2]]
    }
    
    a_j = T_j - logphi+ (nu-2)/2 * log(2) + gamma_1 - 1/2*log(pi)

    const = -logphi - (n-1)*log(nu) - logdet.sigma/2 - n/2 * log(pi) + (nu-2)/2 * log(2) + gamma_n
    func <- function(idx){
        log_x = log(x[idx,])
        x_j = (log_x + a_j) * 1/nu
        x_circ = exp(x_j)
        val = sum(x_j) - sum(log_x) + log(c(t(x_circ) %*% inv.sigma %*% x_circ)) * (-nu-n)/2 + const
        return(val)
    }
    if(!is.null(ncores)){
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores,mc.preschedule = TRUE))

    }else{
        val = unlist(lapply(1:nrow(x),func))
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
    if(log) return(val)
    else return(exp(val))
} 

## this function returns the exponent function of the
## truncated extremal-t max-stable processes
V_truncT <- function(x,par,T_j=NULL,ncores=NULL){
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    sigma = par[[1]];nu = par[[2]]
    n = ncol(x)
    if(n==1){return(1/x)}
    gamma_1 = gamma((nu+1)/2)
    sigma_fun <- function(j){
        sigma_j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F]/sigma[j,j])/ (nu + 1)/sigma[j,j]
    }
    a_fun <- function(j,upper,sigma_j){
        val = mvtnorm::pmvt(lower=-sigma[-j,j]/sigma[j,j],upper=upper-sigma[-j,j]/sigma[j,j],sigma=sigma_j,df=nu+1)[[1]]
        return(val)
    }
    idx.finite = which(apply(x,2,function(x.i){any(is.finite(x.i))}))
    sigma_j = lapply(1:n,sigma_fun)
    if(is.null(T_j)){
        phi = mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),sigma=sigma)[[1]]
        T_j = rep(1,n)
        if(!is.null(ncores)){ 
            T_j[idx.finite] = mcmapply(a_fun,sigma_j=sigma_j[idx.finite],j=idx.finite,MoreArgs = list(upper=rep(Inf,n-1)),SIMPLIFY = TRUE,mc.cores=ncores,mc.set.seed = FALSE)
        }else{
            T_j[idx.finite] = mapply(a_fun,sigma_j=sigma_j[idx.finite],j=idx.finite,MoreArgs = list(upper=rep(Inf,n-1)),SIMPLIFY = TRUE)
        } 
    }else{
        phi = exp(T_j[[2]])
        T_j = exp(T_j[[1]])
        
    }
    a_j = rep(1,n)
    a_j[idx.finite] = T_j[idx.finite]/phi*2^(nu/2-1)*gamma_1/sqrt(pi)*diag(sigma)[idx.finite]^(nu/2)    
    func <- function(idx){
        x_j = x[idx,] * a_j
        idx.finite.j = which(is.finite(x_j))
        x_upper = lapply(idx.finite.j,function(i){(x_j[-i]/x_j[i])^(1/nu)})
        V_j = mapply(a_fun,j=idx.finite.j,upper=x_upper,sigma_j = sigma_j[idx.finite.j],SIMPLIFY = TRUE)
        tryCatch({val <- sum(V_j/T_j[idx.finite.j]/x[idx,idx.finite.j])},error=function(e) browser())
    }
    if(!is.null(ncores)){
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores,mc.preschedule = TRUE))
    }else{
        val = unlist(lapply(1:nrow(x),func))
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
    return(val)
}

## this function returns the partial derivatives of the exponent function
## for the truncated extremal-t max-stable processes
partialV_truncT <- function(x,idx,par,T_j=NULL,ncores=NULL,log=TRUE){
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    sigma = par[[1]];nu = par[[2]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    k = length(idx)
    if(k==0){
        val = V_truncT(x,par,T_j=T_j,ncores=ncores)
        if(log) return(log(val))
        else return(val)
    }
    if(k==n){
       val = intensity_truncT(x,par,T_j=T_j,ncores,log)
       return(val)
    }
    chol.sigma.11 = chol(sigma[idx,idx])
    inv.sigma.11 = chol2inv(chol.sigma.11)
    logdet.sigma.11 = sum(log(diag(chol.sigma.11)))*2
    sigma_T = (sigma[-idx,-idx] - sigma[-idx,idx,drop=F] %*% inv.sigma.11 %*% sigma[idx,-idx,drop=F])/(nu+k)
    gamma_1 = log(gamma((nu+1)/2))
    gamma_k = log(gamma((k+nu)/2))
    a_fun <- function(j){
        sigma_j = (sigma[-j,-j] - sigma[-j,j,drop=FALSE] %*% sigma[j,-j,drop=FALSE])/(nu + 1)
        val = mvtnorm::pmvt(lower=-sigma[-j,j],upper=rep(Inf,n-1),sigma=sigma_j,df=nu+1)[[1]]
        return(log(val))
    }
    if(is.null(T_j)){
        logphi = log(mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),sigma=sigma)[[1]])
        if(!is.null(ncores)) T_j = unlist(mclapply(1:n,a_fun,mc.cores=ncores,mc.set.seed = FALSE,mc.preschedule = TRUE)) else T_j = unlist(lapply(1:n,a_fun))
        
    }else{
        logphi = T_j[[2]]
        T_j = T_j[[1]]
    }
    a_j = T_j - logphi + (nu-2)/2 * log(2) + gamma_1 - 1/2 * log(pi)
    const = - logphi + (nu-2)/2 * log(2) + (1-k)*log(nu) - k/2*log(pi) - 1/2*logdet.sigma.11 + gamma_k + 1/nu* sum(a_j[idx])
    func <- function(idx_j){
        x_j = x[idx_j,]
        logx = log(x_j)
        x_j.circ = (x_j * exp(a_j))^(1/nu)
        Q_sigma = c((x_j.circ[idx] %*% inv.sigma.11 %*% x_j.circ[idx]))
        val = (1/nu-1) * sum(logx[idx]) - (k+nu)*log(Q_sigma)/2 + const
        loc = c(sigma[-idx,idx] %*% inv.sigma.11 %*% (x_j.circ[idx]))
        upper = x_j.circ[-idx]
        val = val + log(mvtnorm::pmvt(lower=-loc,upper=upper-loc,sigma=sigma_T*Q_sigma,df=k+nu)[[1]])
        return(val)
    }
    if(!is.null(ncores)){
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores,mc.preschedule = TRUE))
    }
    else{
        val = unlist(lapply(1:nrow(x),func))
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
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

# this function computes the intensity function 
# for the log skew-normal based max-stable processes
intensity_logskew <- function(x,par,alpha.para=TRUE,log=TRUE){
    sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(n==1) return(1/(x^2))
    omega2 = diag(sigma)
    chol.sigma = chol(sigma)
    inv.sigma = chol2inv(chol.sigma)
    logdet.sigma = sum(log(diag(chol.sigma)))*2
    if(alpha.para){
        alpha = par[[2]]
        delta = c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha))
    }else{
        delta = par[[2]]
        alpha = c(1 - delta %*% inv.sigma %*% delta)^(-1/2) * c(inv.sigma %*% delta)
    }
    a = log(2) + pnorm(delta,log.p=TRUE)
    q = rowSums(inv.sigma)
    sum.q = sum(q);sum.alpha = sum(alpha)
    q.mat = matrix(q,n,n,byrow=TRUE)
    x.log = log(x)
    x.circ = x.log + matrix(a,nrow=nrow(x),ncol=n,byrow=TRUE)
    beta = (1+sum.alpha^2/sum.q)^(-0.5)
    tau.tilde = apply(x.circ,1,function(x.i)  beta * sum((alpha - sum.alpha*q/sum.q) * (x.i + omega2/2))+ beta*sum.alpha/sum.q)
    A = inv.sigma - q %*% t(q)/sum.q
    val = -(n-3)/2 * log(2) - (n-1)/2*log(pi)-1/2*logdet.sigma - 1/2*log(sum.q) - 1/2 * (sum(q*omega2)-1)/sum.q - 1/8*c(omega2 %*% A %*% omega2) 
    val = val - rowSums(x.log) - 1/2 * apply(x.circ,1,function(x.i) c(x.i %*% A %*% x.i) + sum(x.i * (2*q/sum.q + c(A %*% omega2)))) + pnorm(tau.tilde,log.p=TRUE)
    if(log)
        return(val)
    else
        return(exp(val))    
}

## slant paramter nomalized by the variance, where the variance is constant and sum(alpha)=0
# this function computes the intensity function 
# for the log skew-normal based max-stable processes
intensity_logskew_constraint <- function(x,par,alpha.para=TRUE,log=TRUE){
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(n==1) return(1/(x^2))
    omega2 = diag(sigma)
    chol.sigma = chol(sigma)
    inv.sigma = chol2inv(chol.sigma)
    logdet.sigma = sum(log(diag(chol.sigma)))*2
    if(alpha.para){
        alpha = par[[2]]
        delta = c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha))
    }else{
        delta = par[[2]]
        alpha = c(1 - delta %*% inv.sigma %*% delta)^(-1/2) * c(inv.sigma %*% delta)
    }
    a = log(2) + pnorm(delta,log.p=TRUE)
    q = rowSums(inv.sigma)
    sum.q = sum(q)
    q.mat = matrix(q,n,n,byrow=TRUE)
    x.log = log(x)
    x.circ = x.log + matrix(a,nrow=nrow(x),ncol=n,byrow=TRUE)
    tau.tilde = apply(x.circ,1,function(x.i)  sum(alpha * (x.i + omega2/2)))
    A = inv.sigma - q %*% t(q)/sum.q
    val = -(n-3)/2 * log(2) - (n-1)/2*log(pi)-1/2*logdet.sigma - 1/2*log(sum.q) - 1/2 * (sum(q*omega2)-1)/sum.q - 1/8*c(omega2 %*% A %*% omega2) 
    val = val - rowSums(x.log) - 1/2 * apply(x.circ,1,function(x.i) c(x.i %*% A %*% x.i) + sum(x.i * (2*q/sum.q + c(A %*% omega2)))) + pnorm(tau.tilde,log.p=TRUE)
    assign(".Random.seed", oldSeed, envir=globalenv())
    if(log)
        return(val)
    else
        return(exp(val))    
}


## this function computes the exponent function 
## for the log skew-normal based max-stable processes
V_logskew <- function(x,par,alpha.para=TRUE){
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(n==1) return(1/x)
    omega2 = diag(sigma)
    inv.sigma = chol2inv(chol(sigma))
    if(alpha.para){
        alpha = par[[2]]
        delta = c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha))
    }else{
        delta = par[[2]]
        alpha = c(1 - delta %*% inv.sigma %*% delta)^(-1/2) * c(inv.sigma %*% delta)
    }
    phi.delta = pnorm(delta,log.p=TRUE)
    a = log(2) + phi.delta + omega2/2
    I.mat = diag(n-1)
    A.j = lapply(1:n,function(j){
        if(j<n) return(cbind(I.mat[,0:(j-1)],rep(-1,n-1),I.mat[,j:(n-1)]))
        else  return(cbind(I.mat[,0:(j-1)],rep(-1,n-1))) })
    mu.j = lapply(1:n, function(j) c(A.j[[j]] %*% sigma[,j]))    
    sigma.j = lapply(1:n,function(j) A.j[[j]] %*% sigma %*% t(A.j[[j]]))
    x.circ = log(x) + matrix(a,nrow=nrow(x),ncol=n,byrow=TRUE)
    mu.val.j = lapply(1:n,function(j) cbind(x.circ[,-j]- x.circ[,j] - matrix(mu.j[[j]],ncol=n-1,nrow=nrow(x),byrow=FALSE),delta[j]))
    sigma.val.j = lapply(1:n,function(j) {b = c(A.j[[j]] %*% delta); unname(cbind(rbind(sigma.j[[j]],-b),c(-b,1)))})
    val = unlist(lapply(1:nrow(x),function(i) sum(unlist(lapply(1:n,function(j) mvtnorm::pmvnorm(lower=rep(-Inf,n),upper=mu.val.j[[j]][i,],sigma=sigma.val.j[[j]])[[1]]/exp(phi.delta[j])/x[i,j])))))
    assign(".Random.seed", oldSeed, envir=globalenv())
    return(val)
}


V_bi_logskew <- function(x,par,alpha.para=TRUE){
    sigma = par[[1]]
    if(alpha.para){
        alpha = par[[2]]
        delta = c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha))
    }else{
        delta = par[[2]]
    }
    phi.delta.log = pnorm(delta,log.p=TRUE)
    phi.delta = exp(phi.delta.log)
    vario.sqrt = sigma[1,1] + sigma[2,2] - 2*sigma[1,2]
    if(!is.matrix(x)) x <- matrix(x,ncol=2)
    sigma.1 = matrix(c(vario.sqrt,delta[1]-delta[2],delta[1]-delta[2],1),nrow=2)
    sigma.2 = matrix(c(vario.sqrt,delta[2]-delta[1],delta[2]-delta[1],1),nrow=2)
    mu.1 = cbind(phi.delta.log[2]-phi.delta.log[1]+log(x[,2]/x[,1]) + vario.sqrt/2,delta[1])
    mu.2 = cbind(phi.delta.log[1]-phi.delta.log[2]+log(x[,1]/x[,2]) + vario.sqrt/2,delta[2])
    val = unlist(lapply(1:nrow(x),function(i) 1/x[i,1]/phi.delta[1]*mvtnorm::pmvnorm(lower=rep(-Inf,2),upper=mu.1[i,],sigma=sigma.1)[[1]] + 1/x[i,2]/phi.delta[2]*mvtnorm::pmvnorm(lower=rep(-Inf,2),upper=mu.2[i,],sigma=sigma.2)[[1]]))
    return(val)
}


## this function returns the partial derivatives of the exponent function
## for the truncated extremal-t max-stable processes
partialV_logskew <- function(x,idx,par,alpha.para=TRUE,log=FALSE){
    # set a random seed
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(length(idx)==0){
        val = V_logskew(x,par,alpha.para)
        if(log) return(log(val))
        else return(val)
    }
    if(length(idx)==n){
        val = intensity_logskew(x,par,alpha.para,log)
        return(val)
    }
    sigma.inv = chol2inv(chol(sigma))
    omega2 = diag(sigma)
    if(alpha.para){
        alpha = par[[2]]
        delta = c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha))
    }else{
        delta = par[[2]]
        alpha = c(1 - delta %*% sigma.inv %*% delta)^(-1/2) * c(sigma.inv %*% delta)
    }

    a = log(2)  + pnorm(delta,log.p=TRUE)
    q = rowSums(sigma.inv)
    sum.q = sum(q);sum.alpha = sum(alpha)
    A =  sigma.inv - q %*% t(q)/sum.q

    sigma.tilde.inv = A[-idx,-idx,drop=FALSE]
    sigma.tilde = chol2inv(chol(sigma.tilde.inv))

    beta = (1+sum.alpha^2/sum.q)^(-0.5) * (alpha-sum.alpha*q/sum.q)
    alpha.tilde = beta[-idx]
    x.circ  = log(x) + matrix(a,nrow=nrow(x),ncol=n,byrow=TRUE)
    mu.tilde0 = q[-idx]/sum.q + (A %*% omega2/2)[-idx]
    mu.tilde = - t(sigma.tilde %*%  (A[-idx,idx] %*% t(x.circ[,idx,drop=FALSE]) + matrix(mu.tilde0,ncol=nrow(x),nrow=length(mu.tilde0),byrow=FALSE))) # dim: nrow(x)*(-idx) 
    omega2.mat = matrix(omega2/2,nrow=nrow(x),ncol=n,byrow=TRUE)
    alpha0.tilde =  c((x.circ[,idx,drop=FALSE] + omega2.mat[,idx,drop=FALSE]) %*% beta[idx] + (mu.tilde + omega2.mat[,-idx,drop=FALSE]) %*% beta[-idx]) + (1+sum.alpha^2/sum.q)^(-0.5) * sum.alpha/sum.q # nrow(x) 
    b0  = c((1+ alpha.tilde %*% sigma.tilde %*% alpha.tilde)^(-0.5))
    tau.tilde = alpha0.tilde * b0
    phi.tau.tilde = pnorm(tau.tilde)
    intensity.marginal = intensity_logskew(x[,idx,drop=FALSE],par=list(sigma=sigma[idx,idx],delta[idx]),alpha.para=FALSE,log=FALSE) # nrow(x)
    b = c(alpha.tilde %*% sigma.tilde * b0)
    sigma.val = unname(cbind(rbind(sigma.tilde, -b),c(-b,1)))
    val = unlist(lapply(1:nrow(x),function(i) mvtnorm::pmvnorm(lower=rep(-Inf,ncol(mu.tilde)+1),upper=c(x.circ[i,-idx]-mu.tilde[i,],tau.tilde[i]),sigma=sigma.val)[[1]]/phi.tau.tilde[i]*intensity.marginal[i]))
    assign(".Random.seed", oldSeed, envir=globalenv())
    if(log){
        return(log(val))
    }
    else{
        return(val)
    }
}

# return a covariance matrix for the HR model conditional on the location i
vario.func.i <- function(loc,par,i){ ##return a covariance matrix
    lambda = par[1];alpha = par[2]
    if(!is.matrix(loc)){stop("loc must be a matrix")}
    if(nrow(loc)==2){
        loc.i = loc[-i,] - loc[i,]
    }else{
        loc.i =  t(apply(loc[-i,],1, function(x) x - loc[i,]))
    }
    if(!is.matrix(loc.i)){loc.i = matrix(loc.i,ncol=2,nrow=1)}
    n = nrow(loc.i)
    if(n==1){
        val=(sqrt(sum(loc.i[1,]^2))/lambda)^alpha
        return(list(2*val,val))
    }
    vario <- function(coord){
        if(!is.matrix(coord)) {val <- (sqrt(sum(coord^2))/lambda)^alpha}
        else {val <- (sqrt(sum((coord[1,]-coord[2,])^2))/lambda)^alpha}
        return(val)
    }
    all.pairs = Rfast::comb_n(1:n,2)
    gamma.vec=unlist(mclapply(1:ncol(all.pairs),function(i) {idx = all.pairs[,i];vario(loc[idx,])} ,mc.cores=ncores,mc.preschedule = TRUE))
    gamma.origin = unlist(mclapply(1:n,function(i) vario(loc[i,]),mc.cores=ncores,mc.preschedule = TRUE))
    cov.mat = matrix(0,n,n);
    cov.mat[t(all.pairs)] <- gamma.vec
    cov.mat[t(all.pairs[2:1,])] <- gamma.vec
    cov.mat = matrix(gamma.origin,n,n,byrow=TRUE) + matrix(gamma.origin,n,n,byrow=FALSE) - cov.mat      
    return(cov.mat)
}

# intensity function for HR model with Frechet margins
intensity_HR <- function(data,par,i){
     cov.mat <- cov.transform.i(par[[1]],i)
     gamma.origin <- diag(cov.mat)/2
     if(length(gamma.origin)==1){val1 = dnorm(log(data[-i])-log(data[i]),mean=-gamma.origin,sd=sqrt(cov.mat),log=TRUE)
     }else{val1 = mvtnorm::dmvnorm(x=log(data[-i])-log(data[i]),mean=-gamma.origin,sigma=cov.mat,log=TRUE)}
     val = val1-sum(log(data))-log(data[i])
     return(exp(val))
}

# exponent function for HR model with Frechet margins
V_HR <- function(data,par,i){
    set.seed(342424)
    cov.mat <- cov.transform.i(par[[1]],i)
    gamma.origin <- diag(cov.mat)/2
    if(length(gamma.origin)==1){
        func <- function(r){
            1 - pnorm(log(data[-i])+log(r),mean=-gamma.origin,sd=sqrt(cov.mat))
        }
    }else{
        func <- function(r){
            func.i <- function(r.i){
                oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
                set.seed(19873436)
                val = 1 - mvtnorm::pmvnorm(lower=rep(-Inf,nrow(cov.mat)),upper=log(data[-i])+log(r.i),mean=-gamma.origin,sigma=cov.mat)[[1]]
                assign(".Random.seed", oldSeed, envir=globalenv())
                return(val)
            }
            vapply(r,func.i,numeric(1))
        }
    }
    val = integrate(func,lower=1/data[i],upper=Inf,rel.tol=1e-5,subdivisions=1e4)$value + 1/data[i]
    return(val)
}

cov.transform.i <- function(cov.mat,i){
    d = nrow(cov.mat)
    A = matrix(0,ncol=d,nrow=d-1)
    A[,-i] = diag(d-1)
    A[,i] = -1
    cov.mat.i = A %*% cov.mat %*% t(A)
}

# intensity function for skewed HR model with Frechet margins
intensity_skewedHR <- function(data,par,i){
    gamma.origin <- diag(par[[1]])/2
    delta = par[[2]]
    a = log(2) + gamma.origin + pnorm(par[[2]],log.p=TRUE) - par[[1]][i,]
    a = a[-i] - a[i]
    cov.mat.i <- cov.transform.i(par[[1]],i)
    cov.mat.inv = chol2inv(chol(cov.mat.i))
    
    tau = delta[i]
    delta = delta[-i]-delta[i]
    alpha = c(1 - t(delta) %*% cov.mat.inv %*% delta)^(-1/2) * c(cov.mat.inv %*% delta)
    alpha0 = tau * c(1  + t(alpha) %*% cov.mat.i %*% alpha)^0.5

    if(length(alpha)==1){
        x = log(data[-i])-log(data[i])
        val1 = dnorm(x,mean=-a,sd=sqrt(cov.mat.i),log=TRUE) + pnorm(sum(alpha * (x + a))+alpha0,log=TRUE) - pnorm(tau,log=TRUE) 
    }else{
        x = log(data[-i])-log(data[i])
        val1 = mvtnorm::dmvnorm(x,mean=-a,sigma=cov.mat.i,log=TRUE) + pnorm(sum(alpha * (x + a)) + alpha0,log=TRUE) - pnorm(tau,log=TRUE)
    }
    val = val1-sum(log(data))-log(data[i])
    return(exp(val))
}

V_skewedHR <-  function(data,par,i){
    gamma.origin <- diag(par[[1]])/2
    delta = par[[2]]
    a = log(2) + gamma.origin + pnorm(par[[2]],log.p=TRUE) - par[[1]][i,]
    a = a[-i] - a[i]
    cov.mat.i <- cov.transform.i(par[[1]],i)
    tau = delta[i]
    delta = delta[-i]-delta[i]
    sigma  = unname(cbind(rbind(cov.mat.i,-delta),c(-delta,1)))
    p.tau = pnorm(tau)
    x = log(data[-i]) + a
    func <- function(r){
        func.i <- function(r.i){
            oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
            set.seed(19873436)
            val = 1 - mvtnorm::pmvnorm(lower=rep(-Inf,nrow(sigma)),upper=c(x+log(r.i),tau),sigma=sigma)[[1]]/p.tau
            assign(".Random.seed", oldSeed, envir=globalenv())
            return(val)
        }
        sapply(r,func.i,simplify = TRUE)
    }
    val = integrate(func,lower=1/data[i],upper=Inf,rel.tol=1e-5,subdivisions=1e5)$value + 1/data[i]
    return(val)
}

# calculate empirical extremal coefficients: returns the MLE estimator (see page 374 of the lecture notes).
empirical_extcoef <- function(idx,data){
    return(min(2,max(1,1/mean(1/pmax(data[,idx[1]],data[,idx[2]])))))
}

alpha2delta <- function(par){
    alpha = par[[2]];sigma = par[[1]]
    delta = c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha))
    return(list(sigma,delta))
}

delta2alpha <- function(par){
    delta = par[[2]];sigma = par[[1]]
    chol.sigma = chol(sigma)
    inv.sigma = chol2inv(chol.sigma)
    alpha = c(1 - delta %*% inv.sigma %*% delta)^(-1/2) * c(inv.sigma %*% delta)
    return(list(sigma,alpha))
}

# calculate true extremal coefficients
true_extcoef <- function(idx,par,model="logskew1",T_j=NULL){
    if(model=="logskew1"){
        delta = par[[2]];sigma = par[[1]]
        val = V_logskew(rep(1,length(idx)),list(sigma=sigma[idx,idx],delta[idx]),alpha.para=FALSE)
        
    }
    if(model == "truncT1"){
        n = nrow(par[[1]])
        x = matrix(unlist(lapply(idx,function(id){x.id = rep(Inf,n); x.id[id] = 1;x.id})),ncol=n,byrow=TRUE)
        val = V_truncT(x,par,T_j=T_j,ncores=min(detectCores(),nrow(x)))
    }
    
    if(model == "truncT2"){
        x = rep(1,length(idx))
        val = V_truncT(x,list(sigma=par[[1]][idx,idx],nu = par[[2]]),T_j=NULL,ncores=NULL)
    }
    
    if(model == "BR"){
        x = matrix(1,ncol=length(idx))
        if(length(idx)==2){
            val = V.biv(x,sigma=par[[1]][idx,idx])
        }else{
            val = V(x,sigma=par[[1]][idx,idx])
        }
    }
    return(val)
}

# cov.func <- function(loc,par){
#     r = par[1];v = par[2]
#     n = nrow(loc)
#     diff.vector <- cbind(as.vector(outer(loc[,1],loc[,1],'-')),
#         as.vector(outer(loc[,2],loc[,2],'-')))
#     cov.mat <- matrix(exp(-(sqrt(diff.vector[,1]^2 + diff.vector[,2]^2)/r)^v), ncol=n) + diag(1e-6,n) 
#     #cov.mat <- diag(seq(1,2,length.out=n)) %*% cov.mat %*% diag(seq(1,2,length.out=n))
#     return(cov.mat)
# }

cov.func <- function(distmat,par,v=1){
    r = par[1];shape=par[2]
    if(length(par)==3){v=par[3]}
    n=nrow(distmat)
    cov.mat <- v*exp(-(distmat/r)^shape)
    return(cov.mat + .Machine$double.eps * diag(n))
}

alpha.func <- function(par,b.mat=basis){
    alpha <- c(par %*% t(b.mat))
    return(alpha)
}

## inference for simulated data ##  
fit.model <- function(data,loc,init,fixed=NULL,model="truncT",maxit=100,FUN=NULL,basis=NULL,alpha.func=NULL,
                    ncores=NULL,method="L-BFGS-B",lb=NULL,ub=NULL,hessian=FALSE,opt=FALSE,trace=FALSE,step2=TRUE,idx.para=1:2,pareto=FALSE,partial=FALSE){
    t0 <- proc.time()
    if(is.null(fixed)){fixed = rep(FALSE,length(init))}
    if(is.null(lb)){lb=rep(-Inf,length(init))}
    if(is.null(ub)){ub=rep(Inf,length(init))}
    fixed2 = fixed
    if(model == "logskew"){
    ## 5 parameters: 2 for the covariance function; 3 for the slant parameter
        # fixed2[-idx.para] = TRUE
        object.func <- function(par,opt=TRUE,ncore=NULL){
            if(any(par < lb[!fixed2]) | any(par > ub[!fixed2])){return(Inf)}
            par2 = init; par2[!fixed2] = par
            par.1 = par2[idx.para];par.2 = par2[-idx.para]
            sigma = FUN(loc,par.1,ncore)
            alpha = alpha.func(par=par.2,b.mat= basis)
            if(!partial){
                val = intensity_logskew_constraint(data,par=list(sigma,alpha),log=TRUE) 
            }else{
                delta = c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha))
                omega2 = diag(sigma)
                a = log(2) + pnorm(delta,log.p=TRUE)
                computeScores <- function(i){
                    ind = data[[i]][[1]]
                    x = data[[i]][[2]]
                    n = length(x)
                    if(n==1) return(1/(x^2))
                    sigma.i = sigma[ind,ind]
                    chol.sigma = chol(sigma.i)
                    inv.sigma = chol2inv(chol.sigma)
                    logdet.sigma = sum(log(diag(chol.sigma)))*2
                    delta.i = delta[ind]
                    omega2.i = omega2[ind]
                    alpha.i = c(1 - delta.i %*% inv.sigma %*% delta.i)^(-1/2) * c(inv.sigma %*% delta.i)
                    a.i = a[ind]
                    q = rowSums(inv.sigma)
                    sum.q = sum(q);sum.alpha = sum(alpha)
                    q.mat = matrix(q,n,n,byrow=TRUE)
                    x.log = log(x)
                    x.circ = x.log + a.i
                    beta = (1+sum.alpha^2/sum.q)^(-0.5)
                    tau.tilde =  beta * sum((alpha.i - sum.alpha*q/sum.q) * (x.circ + omega2.i/2))+ beta*sum.alpha/sum.q
                    A = inv.sigma - q %*% t(q)/sum.q
                    val = -(n-3)/2 * log(2) - (n-1)/2*log(pi)-1/2*logdet.sigma - 1/2*log(sum.q) - 1/2 * (sum(q*omega2.i)-1)/sum.q - 1/8*c(omega2.i %*% A %*% omega2.i) 
                    val = val - x.log - 1/2 * (c(x.circ %*% A %*% x.circ) + sum(x.circ * (2*q/sum.q + c(A %*% omega2.i)))) + pnorm(tau.tilde,log.p=TRUE) 
                    return(val)
                }
                if(!is.null(ncore)){
                    val = unlist(parallel::mclapply(1:length(data),computeScores,mc.cores = ncore,mc.preschedule = TRUE))
                }
                else{
                    val = unlist(lapply(1:length(data),computeScores))
                }
            }
            return(-mean(val)) 
        }
    }
    if(model == "truncT"){
    ## 3 parameters: 2 for the covariance function; 1 for the df parameter
        object.func <- function(par,opt=TRUE,ncore=NULL){
            #if(trace) print(par)
            par2 = init; par2[!fixed] = par
            par.1 = par2[idx.para];nu = par2[-idx.para]
            if(any(par < lb[!fixed]) | any(par > ub[!fixed])){return(Inf)}
            cov.mat = cov.func(loc,par.1)
            para.temp = list(sigma=cov.mat,nu=nu)
            val = intensity_truncT(data,par=para.temp,T_j=a_fun(para.temp,ncores=ncore),log=TRUE,ncores=ncore) 
            return(-mean(val)) 
        }
    }
    if(opt){
        if(method=="L-BFGS-B"){
            opt.result = optim(init[!fixed2],lower=lb[!fixed2],upper=ub[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace),hessian=hessian,ncore=ncores)
        }else{
            opt.result = optim(init[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace,reltol=1e-8),hessian=hessian,ncore=ncores)
        }
        if(model=="logskew" & any(!fixed[-idx.para]) & step2 & !is.null(ncores)){
            n.alpha = sum(!fixed[-idx.para])
            ncores.2 = ceiling(ncores/4)
            if(n.alpha==2){
                a = seq(0,2*pi,length.out=4)
                a = cbind(cos(a),sin(a))
            } else {
                a = matrix(rnorm(ncores*n.alpha),ncol=n.alpha)
                a = sweep(a,1,sqrt(rowSums(a^2)),FUN="/")
                a[,1] = pmax(a[,1],1)
            }
            init[!fixed2] = opt.result$par
            fixed2[-idx.para] = fixed[-idx.para]
            fixed2[idx.para] = TRUE;
            init.list = split(a,row(a)) 
            if(method=="L-BFGS-B"){
                opt.result2 = mcmapply(optim,par=init.list,MoreArgs = list(fn=object.func,lower=lb[!fixed2],upper=ub[!fixed2],method=method,control=list(maxit=maxit,trace=FALSE),hessian=FALSE,ncore=ncores.2),mc.cores=4,mc.set.seed=FALSE,SIMPLIFY=FALSE)
            }else{
                opt.result2 = mcmapply(optim,par=init.list,MoreArgs = list(fn=object.func,method=method,control=list(maxit=maxit,trace=FALSE,reltol=1e-8),hessian=FALSE,ncore=ncores.2),mc.cores=4,mc.set.seed=FALSE,SIMPLIFY=FALSE)
            }
            opt.values <- unlist(lapply(opt.result2,function(x){tryCatch(x$value,error=function(e){return(Inf)})}))
            opt.result = opt.result2[[which.min(opt.values)]]
            init[!fixed2] = opt.result$par
            fixed2 = fixed;fixed2[-idx.para]=TRUE
            if(method=="L-BFGS-B"){
                opt.result = optim(init[!fixed2],lower=lb[!fixed2],upper=ub[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace),hessian=hessian,ncore=ncores)
            }else{
                opt.result = optim(init[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace,reltol=1e-8),hessian=hessian,ncore=ncores)
            }
            opt.result$others = opt.result2
        }
    }else{
        opt.result = list(par=init,value=object.func(init,opt=FALSE))
    }
    opt.result$time <- proc.time() - t0
    if(hessian){
        h = 1e-01
        opt.result$grad = numDeriv::jacobian(object.func,opt.result$par,opt=FALSE,method="simple")
        opt.result$hessian = numDeriv::hessian(object.func,opt.result$par,opt=TRUE)
        opt.result$K = var(opt.result$grad)
    }
    par2 = init; par2[!fixed2] = opt.result$par
    opt.result$par = par2
    return(opt.result)
}

# function to create list of lists: 
# x: number of lists in each layer 
# n: number of layers
create_lists <- function(x){
    if(length(x)==0){
        return(list())
    }else{
        return(lapply(1:x[1],function(i){create_lists(x[-1])}))
    }
}

vario.func <- function(loc,par,ncores=NULL){ ##return a covariance matrix
    lambda = par[1];alpha = par[2]
    if(!is.matrix(loc)){loc = matrix(loc,nrow=1)}
    n = nrow(loc)
    if(n==1){
        val=(sqrt(sum(loc[1,]^2))/lambda)^alpha
        return(val)
    }
    vario <- function(coord){
        if(!is.matrix(coord)) {val <- (sqrt(sum(coord^2))/lambda)^alpha}
        else {val <- (sqrt(sum((coord[1,]-coord[2,])^2))/lambda)^alpha}
        return(val)
    }
    all.pairs = Rfast::comb_n(1:n,2)
    if(is.null(ncores)){
        gamma.vec=unlist(lapply(1:ncol(all.pairs),function(i) {idx = all.pairs[,i];vario(loc[idx,])}))
        gamma.origin = unlist(lapply(1:n,function(i) vario(loc[i,])))
    }else{
        gamma.vec=unlist(mclapply(1:ncol(all.pairs),function(i) {idx = all.pairs[,i];vario(loc[idx,])} ,mc.cores=ncores,mc.preschedule = TRUE))
        gamma.origin = unlist(mclapply(1:n,function(i) vario(loc[i,]),mc.cores=ncores,mc.preschedule = TRUE))
    }
    cov.mat = matrix(0,n,n);
    cov.mat[t(all.pairs)] <- gamma.vec
    cov.mat[t(all.pairs[2:1,])] <- gamma.vec
    cov.mat = matrix(gamma.origin,n,n,byrow=TRUE) + matrix(gamma.origin,n,n,byrow=FALSE) - cov.mat
    return(cov.mat + .Machine$double.eps * diag(n))
}

vario.func2 <- function(loc,par,ncores=NULL){
    lambda = par[1];alpha = par[2];theta = par[3];a = par[4]
    if(!is.matrix(loc)){loc = matrix(loc,nrow=1)}
    Omega = matrix(c(cos(theta),a*sin(theta),-sin(theta),a*cos(theta)),nrow=2,ncol=2)
    n = nrow(loc)
    Omega = t(Omega) %*% Omega
    if(n==1){
        val=(sqrt( c(loc[1,] %*% Omega %*% loc[1,]) )/lambda)^alpha
        return(val)
    }
    vario <- function(coord){
        if(!is.matrix(coord)){
            val <- (sqrt( c(coord %*% Omega %*% coord) )/lambda)^alpha
        }else{
            coord = coord[1,] - coord[2,]
            val <- (sqrt( c(coord %*% Omega %*% coord) )/lambda)^alpha
        }
        return(val)
    }
    all.pairs = Rfast::comb_n(1:n,2)
    if(is.null(ncores)){
        gamma.vec = unlist(lapply(1:ncol(all.pairs),function(i) {idx = all.pairs[,i];vario(loc[idx,])}))
        gamma.origin = unlist(lapply(1:n,function(i) vario(loc[i,])))
    }else{
        gamma.vec=unlist(mclapply(1:ncol(all.pairs),function(i) {idx = all.pairs[,i];vario(loc[idx,])} ,mc.cores=ncores,mc.preschedule = TRUE))
        gamma.origin = unlist(mclapply(1:n,function(i) vario(loc[i,]),mc.cores=ncores,mc.preschedule = TRUE))
    }
    cov.mat = matrix(0,n,n);
    cov.mat[t(all.pairs)] <- gamma.vec
    cov.mat[t(all.pairs[2:1,])] <- gamma.vec
    cov.mat = matrix(gamma.origin,n,n,byrow=TRUE) + matrix(gamma.origin,n,n,byrow=FALSE) - cov.mat
    return(cov.mat + .Machine$double.eps * diag(n))
}
