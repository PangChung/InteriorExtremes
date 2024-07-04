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
intensity_logskew <- function(x,par,alpha.para=TRUE,log=TRUE){
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
    sum.q = sum(q);sum.alpha = sum(alpha)
    q.mat = matrix(q,n,n,byrow=TRUE)
    x.log = log(x)
    x.circ = x.log + matrix(a,nrow=nrow(x),byrow=TRUE)
    beta = (1+sum.alpha^2/sum.q)^(-0.5)
    tau.tilde = apply(x.circ,1,function(x.i)  beta * sum((alpha - sum.alpha*q/sum.q) * (x.i + omega2/2))+ beta*sum.alpha/sum.q)
    A = inv.sigma - q %*% t(q)/sum.q
    val = -(n-3)/2 * log(2) - (n-1)/2*log(pi)-1/2*logdet.sigma - 1/2*log(sum.q) - 1/2 * (sum(q*omega2)-1)/sum.q - 1/8*c(omega2 %*% A %*% omega2) 
    val = val - rowSums(x.log) - 1/2 * apply(x.circ,1,function(x.i) c(x.i %*% A %*% x.i) + sum(x.i * (2*q/sum.q + c(A %*% omega2))))

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
    A.j = lapply(1:n,function(j){
        if(j<n) return(cbind(I.mat2[,0:(j-1)],rep(-1,n-1),I.mat2[,j:(n-1)]))
        else  return(cbind(I.mat2[,0:(j-1)],rep(-1,n-1))) })
    mu.j = lapply(1:n, function(j) c(A.j[[j]] %*% sigma[,j]))    
    sigma.j = lapply(1:n,function(j) A.j[[j]] %*% sigma %*% t(A.j[[j]]))
    x.circ = log(x) + matrix(a,nrow=nrow(x),byrow=TRUE)
    mu.val.j = lapply(1:n,function(j) cbind(x.circ[,-j]- matrix(x.circ[,j]+mu.j[[j]],ncol=n-1,nrow=nrow(x),byrow=TRUE),delta[j]))
    sigma.val.j = lapply(1:n,function(j) b = A.j[[j]] %*% delta; unname(cbind(rbind(sigma.j[[j]],-b),c(-b,1))))
    val = unlist(lapply(1:nrow(x),function(i) sum(unlist(lapply(1:n,function(j) mvtnorm::pmvnorm(lower=rep(-Inf,n),upper=mu.val.j[[j]][i,],sigma=sigma.val.j[[j]])[[1]]/exp(phi.delta[j])/x[i,j])))))
    assign(".Random.seed", oldSeed, envir=globalenv())
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
        val = V_logskew(x,par,alpha.para,ncores=ncores)
        if(log) return(log(val))
        else return(val)
    }
    if(length(idx)==n){
        val = intensity_logskew(x,par,alpha.para,ncores,log)
        return(val)
    }
    sigma.inv = chol2inv(chol(sigma))
    omega2 = diag(sigma)
    if(alpha.para){
        alpha = par[[2]]
        delta = c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha))
    }else{
        delta = par[[2]]
        alpha = c(1 - delta %*% inv.sigma %*% delta)^(-1/2) * c(inv.sigma %*% delta)
    }

    a = log(2)  + pnorm(delta,log.p=TRUE)
    q = rowSums(sigma.inv)
    sum.q = sum(q);sum.alpha = sum(alpha)
    A =  sigma.inv - q %*% t(q)/sum.q

    sigma.tilde.inv = A[-idx,-idx,drop=FALSE]
    sigma.tilde = chol2inv(chol(sigma.tilde.inv))

    beta = (1+sum.alpha^2/sum.q)^(-0.5) * (alpha-sum.alpha*q/sum.q)
    alpha.tilde = beta[-idx]
    x.circ  = log(x) + matrix(a,nrow=nrow(x),byrow=TRUE)
    mu.tilde0 = q[-idx]/sum.q + (A %*% omega2/2)[-idx]
    mu.tilde = - t(sigma.tilde %*%  (A[-idx,idx] %*% t(x.circ[,idx,drop=FALSE]) + matrix(mu.tilde0,ncol=nrow(x),nrow=length(mu.tilde0),byrow=FALSE))) # dim: nrow(x)*(-idx) 
    omega2.mat = matrix(omega2/2,nrow=nrow(x),ncol=n,byrow=TRUE)
    alpha0.tilde =  c((x.circ[,idx,drop=FALSE]+omega2.mat[,idx,drop=FALSE]) %*% beta[idx] + (mu.tilde + omega2.mat[,-idx,drop=FALSE]) %*% beta[-idx]) + (1+sum.alpha^2/sum.q)^(-0.5) * sum.alpha/sum.q # nrow(x) 
    b0  = c((1+ alpha.tilde %*% sigma.tilde %*% alpha.tilde)^(-0.5))
    tau.tilde = alpha0.tilde * b0
    phi.tau.tilde = pnorm(tau.tilde)
    intensity.marginal = intensity_logskew(x[,idx],par=list(sigma=sigma[idx,idx],delta[idx]),alpha.para=FALSE,log=FALSE) # nrow(x)
    b = alpha.tilde %*% sigma.tilde %*% b0
    sigma.val = unname(cbind(rbind(sigma.tilde, -b2),c(-b2,1)))
    val = lapply(1:nrow(x),function(i) mvtnorm::pmvnorm(upper=c(mu.tilde[-i,],tau.tilde[i]),sigma=sigma.val)[[1]]/phi.tau.tilde[i]*intensity.marginal[i]) 
    assign(".Random.seed", oldSeed, envir=globalenv())
    if(log){
        return(log(val))
    }
    else{
        return(val)
    }
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
        val = V_logskew(rep(1,length(idx)),list(sigma=sigma[idx,idx],delta[idx]),alpha.para=FALSE,ncores=NULL)
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

cov.func <- function(distmat,par){
    r = par[1];v = par[2];n=nrow(distmat)
    cov.mat <- v*exp(-(distmat/r)) + diag(1e-6,n) 
    #cov.mat <- diag(seq(1,2,length.out=n)) %*% cov.mat %*% diag(seq(1,2,length.out=n))
    return(cov.mat)
}

alpha.func <- function(par,b.mat=basis){
    alpha <- c(par %*% t(b.mat))
    return(alpha)
}

