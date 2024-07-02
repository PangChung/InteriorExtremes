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
                mvtnorm::pmvt(lower=-sigma[-j,j]/sigma[j,j],upper=rep(Inf,n-1),sigma=sigma_j,df=nu+1)[[1]]},j=1:n,mc.cores=ncores)
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
        if(!is.null(ncores)) T_j = unlist(mclapply(1:n,a_fun,mc.cores=ncores)) else T_j = unlist(lapply(1:n,a_fun))
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
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))

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
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))
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
        if(!is.null(ncores)) T_j = unlist(mclapply(1:n,a_fun,mc.cores=ncores)) else T_j = unlist(lapply(1:n,a_fun))
        
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
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))
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
# intensity_logskew <- function(x,par,alpha.para=TRUE,ncores=NULL,log=TRUE){
#     oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
#     set.seed(747380)
#     sigma = par[[1]]
#     if(!is.matrix(x)){x <- matrix(x,nrow=1)}
#     n = ncol(x)
#     if(n==1) return(1/(x^2))
#     omega = diag(sqrt(diag(sigma)))
#     omega_inv = diag(diag(omega)^(-1))
#     sigma_bar = omega_inv %*% sigma %*% omega_inv
#     chol.sigma = chol(sigma)
#     inv.sigma = chol2inv(chol.sigma)
#     logdet.sigma = sum(log(diag(chol.sigma)))*2
#     if(alpha.para){
#         alpha = par[[2]]
#         delta = c(sigma_bar %*% alpha)/sqrt(c(1+alpha %*% sigma_bar %*% alpha))
#     }else{
#         inv.sigma.bar = omega %*% inv.sigma %*% omega
#         delta = par[[2]]
#         alpha = c(1 - delta %*% inv.sigma.bar %*% delta)^(-1/2) * c(inv.sigma.bar %*% delta)
#     }
#     a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)
#     sum.inv.sigma = sum(inv.sigma)
#     one.mat = matrix(1,n,n)
#     one.vec = rep(1,n)
#     b = c(alpha %*% omega_inv %*% one.vec/sqrt(sum.inv.sigma))
#     beta.hat =  c(alpha %*% omega_inv %*% (diag(n) - one.mat %*% inv.sigma/sum.inv.sigma) * (1+b^2)^(-1/2))
#     A = inv.sigma - inv.sigma %*% one.mat %*% inv.sigma/sum.inv.sigma 
#     delta.hat = (1+b^2)^(-1/2)*b
#     func <- function(idx){
#         x_log = log(x[idx,])
#         x_circ = x_log + a
#         val = -(n-3)/2 * log(2) - (n-1)/2*log(pi) + pnorm(beta.hat %*% x_circ + delta.hat/sqrt(sum.inv.sigma),log.p=TRUE) -
#             1/2*logdet.sigma - 1/2*log(sum.inv.sigma) - sum(x_log)
#         val = val - 1/2 * c(x_circ %*% A %*% x_circ) - c(one.vec %*% inv.sigma %*% x_circ)/sum.inv.sigma + 1/2/sum.inv.sigma
#         return(val)
#     }
    
#     if(!is.null(ncores)){
#         val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))
#     }
#     else{
#         val = unlist(lapply(1:nrow(x),func))
#     }
#     assign(".Random.seed", oldSeed, envir=globalenv())
#     if(log)
#         return(val)
#     else
#         return(exp(val))    
# }

## slant parameter without normlized by the variance
# intensity_logskew <- function(x,par,alpha.para=TRUE,ncores=NULL,log=TRUE){
#     oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
#     set.seed(747380)
#     sigma = par[[1]]
#     if(!is.matrix(x)){x <- matrix(x,nrow=1)}
#     n = ncol(x)
#     if(n==1) return(1/(x^2))
#     chol.sigma = chol(sigma)
#     inv.sigma = chol2inv(chol.sigma)
#     logdet.sigma = sum(log(diag(chol.sigma)))*2
#     if(alpha.para){
#         alpha = par[[2]]
#         delta = c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha))
#     }else{
#         delta = par[[2]]
#         alpha = c(1 - delta %*% inv.sigma %*% delta)^(-1/2) * c(inv.sigma %*% delta)
#     }
#     a = log(2) + diag(sigma)/2 + sapply(delta,pnorm,log.p=TRUE)
#     sum.inv.sigma = sum(inv.sigma)
#     one.mat = matrix(1,n,n)
#     one.vec = rep(1,n)
#     beta.hat =  c(alpha %*%  (diag(n) - one.mat %*% inv.sigma/sum.inv.sigma))
#     A = inv.sigma - inv.sigma %*% one.mat %*% inv.sigma/sum.inv.sigma 
#     func <- function(idx){
#         x_log = log(x[idx,])
#         x_circ = x_log + a
#         val = -(n-3)/2 * log(2) - (n-1)/2*log(pi) + pnorm(beta.hat %*% x_circ,log.p=TRUE) -
#             1/2*logdet.sigma - 1/2*log(sum.inv.sigma) - sum(x_log)
#         val = val - 1/2 * c(x_circ %*% A %*% x_circ) - c(one.vec %*% inv.sigma %*% x_circ)/sum.inv.sigma + 1/2/sum.inv.sigma
#         return(val)
#     }
    
#     if(!is.null(ncores)){
#         val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))
#     }
#     else{
#         val = unlist(lapply(1:nrow(x),func))
#     }
#     assign(".Random.seed", oldSeed, envir=globalenv())
#     if(log)
#         return(val)
#     else
#         return(exp(val))    
# }

## slant paramter nomalized by the variance, where the variance is constant and sum(alpha)=0
intensity_logskew <- function(x,par,alpha.para=TRUE,ncores=NULL,log=TRUE){
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(n==1) return(1/(x^2))
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
    a = log(2) + diag(sigma)/2 + sapply(delta,pnorm,log.p=TRUE)
    sum.inv.sigma = sum(inv.sigma)
    one.mat = matrix(1,n,n)
    one.vec = rep(1,n)
    A = inv.sigma - inv.sigma %*% one.mat %*% inv.sigma/sum.inv.sigma 
    func <- function(idx){
        x_log = log(x[idx,])
        x_circ = x_log + a
        val = -(n-3)/2 * log(2) - (n-1)/2*log(pi) + pnorm(c(alpha %*% x_circ),log.p=TRUE) -
            1/2*logdet.sigma - 1/2*log(sum.inv.sigma) - sum(x_log)
        val = val - 1/2 * c(x_circ %*% A %*% x_circ) - c(one.vec %*% inv.sigma %*% x_circ)/sum.inv.sigma + 1/2/sum.inv.sigma
        return(val)
    }
    
    if(!is.null(ncores)){
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))
    }
    else{
        val = unlist(lapply(1:nrow(x),func))
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
    if(log)
        return(val)
    else
        return(exp(val))    
}

## this function computes the exponent function 
## for the log skew-normal based max-stable processes
V_logskew <- function(x,par,alpha.para=TRUE,ncores=NULL){
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(n==1) return(1/x)
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
    a = log(2) + diag(sigma)/2 + sapply(delta,pnorm,log.p=TRUE)
    I.mat1 = diag(rep(1,n))
    I.mat2 = diag(rep(1,n-1))
    func <- function(j){        
        if(j<n){
            A.j = cbind(I.mat2[,0:(j-1)],rep(-1,n-1),I.mat2[,j:(n-1)])
        }else{
            A.j = cbind(I.mat2[,0:(j-1)],rep(-1,n-1))
        }
        sigma.j = A.j %*% sigma %*% t(A.j)
        u.j = c(A.j %*% sigma[,j])
        b3 = c(A.j %*% delta)
        sigma_circ = unname(cbind(rbind(sigma.j,-b3),c(-b3,1)))
        func_temp <- function(i){
            xi = x[i,]
            mu = c(a[-j] - a[j] + log(xi[-j]/xi[j])-u.j,delta[j])
            val = pnorm(delta[j])^(-1)/xi[j] * mvtnorm::pmvnorm(lower=rep(-Inf,n),upper=mu,sigma=sigma_circ)[[1]]
            return(val)
        }
        val = unlist(lapply(1:nrow(x),func_temp))
        return(val)
    }
    if(!is.null(ncores)){
        val = parallel::mclapply(1:n,func,mc.cores = ncores)
        val = matrix(unlist(val),nrow=n,byrow=TRUE)
        val = apply(val,2,sum)
    }
    else{
        val = lapply(1:n,func)
        val = matrix(unlist(val),nrow=n,byrow=TRUE)
        val = apply(val,2,sum)
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
    return(val)
}

V_bi_logskew <- function(x,par){
    delta = par[[2]];sigma = par[[1]]
    phi.delta = pnorm(delta)
    phi.delta.log = log(phi.delta)
    vario.sqrt = sigma[1,1] + sigma[2,2] - 2*sigma[1,2]
    if(!is.matrix(x)) x <- matrix(x,ncol=2)
    sigma.1 = matrix(c(vario.sqrt,delta[1]-delta[2],delta[1]-delta[2],1),nrow=2)
    sigma.2 = matrix(c(vario.sqrt,delta[2]-delta[1],delta[2]-delta[1],1),nrow=2)

    mu.1 = cbind(phi.delta.log[2]-phi.delta.log[1]+log(x[,2]/x[,1]) + vario.sqrt/2,delta[1])
    mu.2 = cbind(phi.delta.log[1]-phi.delta.log[2]+log(x[,1]/x[,2]) + vario.sqrt/2,delta[2])
    
    val = unlist(lapply(1:nrow(x),function(i) 1/x[i,1]/phi.delta[1]*mvtnorm::pmvnorm(lower=rep(-Inf,2),upper=mu.1[i,],sigma=sigma.1)[[1]] + 1/x[i,2]/phi.delta[2]*mvtnorm::pmvnorm(lower=rep(-Inf,2),upper=mu.2[i,],sigma=sigma.2)[[1]]))
    return(val)
}

## this function returns negative partial derivatives of the exponent function
## for the skewed-Brown-Resnick max-stable processes
partialV_logskew <- function(x,idx,par,alpha.para=TRUE,ncores=NULL,log=FALSE){
    # set a random seed
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(length(idx)==0){
        val = V_logskew(x,par,alpha.para,ncores=ncores)
        if(log) return(log(val))
        else return(val)
    }
    if(length(idx)==n){
        val = intensity_logskew(x,par,alpha.para,ncores,log)
        return(val)
    }
    ones <- rep(1,n)
    one.mat <- matrix(1,n,n)
    I <- diag(1,n)
    sigma.chol = chol(sigma)
    sigma.inv = chol2inv(sigma.chol)
    sum.sigma.inv = sum(sigma.inv)
    
    if(alpha.para){
        alpha = par[[2]]
        delta = c(sigma %*% alpha)/sqrt(c( 1 + alpha %*% sigma %*% alpha))
    }else{
        delta = par[[2]]
        alpha = c(1 - delta %*% sigma.inv %*% delta)^(-1/2) * c(sigma.inv %*% delta)
    }

    a = log(2) + sapply(delta,pnorm,log.p=TRUE)
    omega.tilde = diag(sigma)
    H =  sigma.inv - (sigma.inv %*% one.mat %*% sigma.inv/sum.sigma.inv)
    sigma.tilde = chol2inv(chol(H[-idx,-idx,drop=FALSE]))

    alpha.tilde = alpha[-idx]
    b1 =c((1 + alpha.tilde %*% sigma.tilde %*% alpha.tilde)^(-1/2))

    delta.tilde = b1 * c(sigma.tilde %*% alpha.tilde)
    scale.val = unname(cbind(rbind(sigma.tilde, -delta.tilde),c(-delta.tilde,1)))

    func <- function(i){
        xi.log = log(x[i,])
        xi.tilde = xi.log + a
        mu.tilde = c(-sigma.tilde %*% (H[-idx,idx,drop=FALSE] %*% xi.tilde[idx] + ((sigma.inv %*% ones)/sum.sigma.inv + H%*%omega.tilde/2)[-idx]))
        tau.tilde = c(b1 * (alpha[idx] %*% (xi.tilde[idx] + omega.tilde[idx]/2)  + alpha[-idx] %*% (mu.tilde + omega.tilde[-idx]/2)))
        mu.val = c(xi.tilde[-idx] - mu.tilde, tau.tilde)
        phi = pnorm(tau.tilde)
        intensity.marginal = c(intensity_logskew(x[i,idx],par=list(sigma[idx,idx],delta[idx]),alpha.para=FALSE,ncores=NULL,log=FALSE))
        val = intensity.marginal * max(mvtnorm::pmvnorm(upper=mu.val,sigma=scale.val)[[1]],0)/phi
        return(val)
    }
    if(!is.null(ncores)){
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))
    }
    else{
        val = unlist(lapply(1:nrow(x),func))
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
    if(log){
        return(log(val))
    }
    else{
        return(val)
    }
}

biv_condnorm <- function(sigma,mu){
    condsigma = sigma[1,1] - sigma[1,2]^2
    fun <- function(x){
        return(dnorm(x,sd=1)*pnorm(mu[1],mean=sigma[1,2]*x,sd=sqrt(condsigma)))
    }
    return(integrate(fun,-Inf,mu[2])$value)
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
        if(length(idx)==2){
            val = V_logskew(rep(1,length(idx)),list(sigma=sigma[idx,idx],delta[idx]),alpha.para=FALSE,ncores=NULL)
        }else{
            val = V_bi_logskew(rep(1,length(idx)),delta[idx],sigma[idx,idx])
        }
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

cov.func <- function(distmat,par){
    r = par[1];v = par[2];shape=par[3]
    n=nrow(distmat)
    cov.mat <- v*exp(-(distmat/r)^shape)
    return(cov.mat + .Machine$double.eps * diag(n))
}

alpha.func <- function(par,b.mat=basis){
    alpha <- c(par %*% t(b.mat))
    return(alpha)
}

## inference for simulated data ##  
fit.model <- function(data,loc,init,fixed=NULL,thres = 50,model="truncT",maxit=100,FUN=NULL,basis=NULL,alpha.func=NULL,
                    ncores=NULL,method="L-BFGS-B",lb=NULL,ub=NULL,hessian=FALSE,opt=FALSE,trace=FALSE,step2=TRUE,idx.para=1:2){
    t0 <- proc.time()
    n = ncol(data)
    data.sum = apply(data,1,sum)
    idx.thres = data.sum > thres*n #& data.sum < 1000*n 
    print(paste("#sampels: ",sum(idx.thres)))
    
    data = sweep(data[idx.thres,],1,data.sum[idx.thres],FUN="/")
    if(is.null(fixed)){fixed = rep(FALSE,length(init))}
    if(is.null(lb)){lb=rep(-Inf,length(init))}
    if(is.null(ub)){ub=rep(Inf,length(init))}
    fixed2 = fixed
    if(model == "logskew"){
    ## 5 parameters: 2 for the covariance function; 3 for the slant parameter
        # fixed2[-idx.para] = TRUE
        object.func <- function(par,opt=TRUE,ncore=NULL){
            #if(trace) print(par)
            par2 = init; par2[!fixed2] = par
            par.1 = par2[idx.para];par.2 = par2[-idx.para]
            cov.mat = FUN(loc,par.1)
            alpha = alpha.func(par=par.2,b.mat= basis)
            if(any(par < lb[!fixed2]) | any(par > ub[!fixed2])){return(Inf)}
            para.temp = list(sigma=cov.mat,alpha=alpha)
            val = intensity_logskew(data,par=para.temp,log=TRUE,ncores=ncore) 
            if(opt) return(-mean(val)) else return(-val)
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
            val = intensity_truncT(data,par=para.temp,T_j=a_fun(para.temp,ncores=ncores),log=TRUE,ncores=ncore) 
            if(opt) return(-mean(val)) else return(-val)
        }
    }
    if(model == "BR"){
    ## 3 parameters: 2 for the covariance function;
        object.func <- function(par,opt=TRUE,ncore=NULL){
            par2 = init; par2[!fixed] = par
            par.1 = par2[idx.para]
            cov.mat = FUN(loc,par.1)
            if(any(par < lb[!fixed]) | any(par > ub[!fixed])){return(Inf)}
            val = nVI(data,cov.mat,1:n,logval=TRUE)
            if(opt) return(-mean(val)) else return(-val)
        }
    }
    if(opt){
        if(method=="L-BFGS-B"){
            opt.result = optim(init[!fixed2],lower=lb[!fixed2],upper=ub[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace),hessian=hessian,ncore=ncores)
        }else{
            opt.result = optim(init[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace),hessian=hessian,ncore=ncores)
        }
        if(model=="logskew" & any(!fixed[-idx.para]) & step2){
            n.alpha = sum(!fixed[-idx.para])
            if(n.alpha==2){
                a = seq(0,2*pi,length.out=ncores)
                a = cbind(cos(a),sin(a))
            } else {
                a = matrix(rnorm(ncores*n.alpha),ncol=n.alpha)
                a = sweep(a,1,sqrt(rowSums(a^2)),FUN="/")
            }
            init[!fixed2] = opt.result$par
            fixed2[-idx.para] = fixed[-idx.para]
            fixed2[idx.para] = TRUE
            a = cbind(matrix(init[idx.para],ncol=length(idx.para),nrow=nrow(a),byrow=T),a)[,!fixed2]
            init.list = split(a,row(a)) 
            if(method=="L-BFGS-B"){
                opt.result2 = mcmapply(optim,par=init.list,MoreArgs = list(fn=object.func,lower=lb[!fixed2],upper=ub[!fixed2],method=method,control=list(maxit=maxit,trace=FALSE),hessian=FALSE),mc.cores=ncores,mc.set.seed=FALSE,SIMPLIFY=FALSE)
            }else{
                opt.result2 = mcmapply(optim,par=init.list,MoreArgs = list(fn=object.func,method=method,control=list(maxit=maxit,trace=FALSE),hessian=FALSE),mc.cores=ncores,mc.set.seed=FALSE,SIMPLIFY=FALSE)
            }
            opt.values <- unlist(lapply(opt.result2,function(x){x$value}))
            opt.result = opt.result2[[which.min(opt.values)]]
            init[!fixed2] = opt.result$par
            fixed2 = fixed;fixed2[-idx.para]=TRUE
            if(method=="L-BFGS-B"){
                opt.result = optim(init[!fixed2],lower=lb[!fixed2],upper=ub[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace),hessian=hessian)
            }else{
                opt.result = optim(init[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace),hessian=hessian)
            }
            opt.result$others = opt.result2
        }
    }else{
        opt.result = list(par=init,value=object.func(init,opt=FALSE))
    }
    opt.result$time <- proc.time() - t0
    if(hessian){
        h = 1e-01
        # par.mat.grad01 = matrix(opt.result$par,nrow=length(opt.result$par),ncol=length(opt.result$par),byrow=TRUE) + diag(h/2,length(opt.result$par))
        # par.mat.grad10 = matrix(opt.result$par,nrow=length(opt.result$par),ncol=length(opt.result$par),byrow=TRUE) - diag(h/2,length(opt.result$par))
        # val.object.grad = lapply(1:sum(!fixed),function(i){(object.func(par.mat.grad01[i,],opt=FALSE) - object.func(par.mat.grad10[i,],opt=FALSE))/h})
        # val.object.grad = matrix(unlist(val.object.grad),ncol = sum(!fixed),byrow=FALSE)
        opt.result$grad = numDeriv::jacobian(object.func,opt.result$par,opt=FALSE,method="simple")
        opt.result$hessian = numDeriv::hessian(object.func,opt.result$par,opt=TRUE)
        opt.result$K = var(opt.result$grad)
    }
    par2 = init; par2[!fixed2] = opt.result$par
    opt.result$par = par2
    return(opt.result)
}

partialV_logskew_num <- function(x,idx,par,alpha.para=TRUE,ncores=NULL,log=FALSE){
    # set a random seed
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(length(idx)==0){
        val = V_logskew(x,par,alpha.para,ncores=ncores)
        if(log) return(log(val))
        else return(val)
    }
    if(length(idx)==n){
        val = intensity_logskew(x,par,alpha.para,ncores,log)
        return(val)
    }
    ones <- rep(1,n)
    one.mat <- matrix(1,n,n)
    I <- diag(1,n)
    sigma.chol = chol(sigma)
    sigma.inv = chol2inv(sigma.chol)
    sum.sigma.inv = sum(sigma.inv)
    
    if(alpha.para){
        alpha = par[[2]]
        delta = c(sigma %*% alpha)/sqrt(c( 1 + alpha %*% sigma %*% alpha))
    }else{
        delta = par[[2]]
        alpha = c(1 - delta %*% sigma.inv %*% delta)^(-1/2) * c(sigma.inv %*% delta)
    }

    a = log(2) + diag(sigma)/2 + sapply(delta,pnorm,log.p=TRUE)

    func <- function(i){
        xi.log = log(x[i,])
        xi.tilde = xi.log + a
        mu.tilde = c(-sigma.tilde %*% (H[-idx,idx,drop=FALSE] %*% xi.log[idx] + ((sigma.inv %*% ones)/sum.sigma.inv + H%*%a)[-idx]))
        tau.tilde = c(b1 * (alpha[idx] %*% xi.log[idx] + alpha %*% a + alpha[-idx] %*% mu.tilde))
        mu.val = c(xi.log[-idx] - mu.tilde, tau.tilde)
        phi = pnorm(tau.tilde)
        intensity.marginal = c(intensity_logskew(x[i,idx],par=list(sigma[idx,idx],delta[idx]),alpha.para=FALSE,ncores=NULL,log=FALSE))
        val = intensity.marginal * max(mvtnorm::pmvnorm(upper=mu.val,sigma=scale.val)[[1]],0)/phi
        return(val)
    } 

    if(!is.null(ncores)){
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))
    }
    else{
        val = unlist(lapply(1:nrow(x),func))
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
    if(log){
        return(log(val))
    }
    else{
        return(val)
    }
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

vario.func <- function(loc,par){ ##return a covariance matrix
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
    all.pairs = combn(1:n,2)
    all.pairs.list = split(all.pairs,col(all.pairs))
    gamma.vec = unlist(lapply(all.pairs.list,function(idx) vario(loc[idx,])))
    gamma.origin = sapply(1:n,function(i) vario(loc[i,]))
    cov.mat = diag(2*gamma.origin)
    cov.mat[t(all.pairs)] <- sapply(1:length(gamma.vec),function(i){
        idx = all.pairs[,i]
        return(gamma.origin[idx[1]] + gamma.origin[idx[2]] - gamma.vec[i])})
    cov.mat[t(all.pairs[2:1,])] <- cov.mat[t(all.pairs)]         
    return(cov.mat + .Machine$double.eps * diag(n))
}

## fit the r-Pareto processes ##
fit.model.pareto <- function(data,loc,init,fixed=NULL,thres = 50,model="truncT",maxit=100,FUN=NULL,basis=NULL,alpha.func=NULL,
                    ncores=NULL,method="L-BFGS-B",lb=NULL,ub=NULL,hessian=FALSE,opt=FALSE,trace=FALSE,step2=TRUE,idx.para=1:2){
    t0 <- proc.time()
    n = ncol(data)
    data.sum = apply(data,1,sum)
    idx.thres = data.sum > thres 
    print(paste("#sampels: ",sum(idx.thres)))
    if(sum(idx.thres)<3){idx.thres = which(order(data.sum,decreasing=TRUE)<=3)}
    data = data[idx.thres,]
    if(is.null(fixed)){fixed = rep(FALSE,length(init))}
    if(is.null(lb)){lb=rep(-Inf,length(init))}
    if(is.null(ub)){ub=rep(Inf,length(init))}
    fixed2 = fixed
    if(model == "logskew"){
    ## 5 parameters: 2 for the covariance function; 3 for the slant parameter
        # fixed2[-idx.para] = TRUE
        object.func <- function(par,opt=TRUE,ncore=NULL){
            #if(trace) print(par)
            par2 = init; par2[!fixed2] = par
            par.1 = par2[idx.para];par.2 = par2[-idx.para]
            cov.mat = FUN(loc,par.1)
            alpha = alpha.func(par=par.2,b.mat= basis)
            if(any(par < lb[!fixed2]) | any(par > ub[!fixed2])){return(Inf)}
            para.temp = list(sigma=cov.mat,alpha=alpha)
            val = intensity_logskew(data,par=para.temp,log=TRUE,ncores=ncore) 
            if(opt) return(-mean(val)) else return(-mean(val))
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
            val = intensity_truncT(data,par=para.temp,T_j=a_fun(para.temp,ncores=ncores),log=TRUE,ncores=ncore) 
            if(opt) return(-mean(val)) else return(-mean(val))
        }
    }
    if(model == "BR"){
    ## 3 parameters: 2 for the covariance function;
        object.func <- function(par,opt=TRUE,ncore=NULL){
            par2 = init; par2[!fixed] = par
            par.1 = par2[idx.para]
            cov.mat = FUN(loc,par.1)
            if(any(par < lb[!fixed]) | any(par > ub[!fixed])){return(Inf)}
            val = nVI(data,cov.mat,1:n,logval=TRUE)
            if(opt) return(-mean(val)) else return(-mean(val))
        }
    }
    if(opt){
        if(method=="L-BFGS-B"){
            opt.result = optim(init[!fixed2],lower=lb[!fixed2],upper=ub[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace),hessian=hessian,ncore=ncores)
        }else{
            opt.result = optim(init[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace),hessian=hessian,ncore=ncores)
        }
        if(model=="logskew" & any(!fixed[-idx.para]) & step2){
            n.alpha = sum(!fixed[-idx.para])
            if(n.alpha==2){
                a = seq(0,2*pi,length.out=ncores)
                a = cbind(cos(a),sin(a))
            } else {
                a = matrix(rnorm(ncores*n.alpha),ncol=n.alpha)
                a = sweep(a,1,sqrt(rowSums(a^2)),FUN="/")
            }
            init[!fixed2] = opt.result$par
            fixed2[-idx.para] = fixed[-idx.para]
            fixed2[idx.para] = TRUE
            a = cbind(matrix(init[idx.para],ncol=length(idx.para),nrow=nrow(a),byrow=T),a)[,!fixed2]
            init.list = split(a,row(a)) 
            if(method=="L-BFGS-B"){
                opt.result2 = mcmapply(optim,par=init.list,MoreArgs = list(fn=object.func,lower=lb[!fixed2],upper=ub[!fixed2],method=method,control=list(maxit=maxit,trace=FALSE),hessian=FALSE),mc.cores=ncores,mc.set.seed=FALSE,SIMPLIFY=FALSE)
            }else{
                opt.result2 = mcmapply(optim,par=init.list,MoreArgs = list(fn=object.func,method=method,control=list(maxit=maxit,trace=FALSE),hessian=FALSE),mc.cores=ncores,mc.set.seed=FALSE,SIMPLIFY=FALSE)
            }
            opt.values <- unlist(lapply(opt.result2,function(x){x$value}))
            opt.result = opt.result2[[which.min(opt.values)]]
            init[!fixed2] = opt.result$par
            fixed2 = fixed;fixed2[-idx.para]=TRUE
            if(method=="L-BFGS-B"){
                opt.result = optim(init[!fixed2],lower=lb[!fixed2],upper=ub[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace),hessian=hessian)
            }else{
                opt.result = optim(init[!fixed2],object.func,method=method,control=list(maxit=maxit,trace=trace),hessian=hessian)
            }
            opt.result$others = opt.result2
        }
    }else{
        return(object.func(init[!fixed],opt,ncores))
    }
    if(hessian){
        h = 1e-3
        par.mat.grad = matrix(opt.result$par,nrow=length(opt.result$par),ncol=length(opt.result$par),byrow=TRUE) + diag(h,length(opt.result$par))
        val.object = object.func(opt.result$par,opt=FALSE)
        val.object.grad = apply(par.mat.grad,1,function(x){(object.func(x,opt=FALSE) - val.object)/h})
        opt.result$K = var(val.object.grad)
        opt.result$hessian.inv = solve(opt.result$hessian)
        opt.result$sigma = opt.result$hessian.inv %*% opt.result$K %*% opt.result$hessian.inv
    }
    par2 = init; par2[!fixed2] = opt.result$par
    opt.result$par = par2
    opt.result$time <- proc.time() - t0
    return(opt.result)
}
