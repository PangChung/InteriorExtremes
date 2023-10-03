#########################################################
###### Intensity function for truncated extremal-t ######
#########################################################

## this function returns the intensity function of the
## truncated extremal-t max-stable processes
intensty_truncT <- function(x,par,parallel=TRUE,ncores=2,log=TRUE){
    sigma = par[[2]];nu = par[[1]]
    n = nrow(sigma)
    chol.sigma = chol(sigma)
    inv.sigma = chol2inv(chol.sigma)
    logdet.sigma = sum(log(diag(chol.sigma)))*2
    logphi = log(mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[1])
    gamma_1 = log(gamma((nu+1)/2))
    gamma_n = log(gamma((nu+d)/2))
    a_fun <- function(j,upper){
        sigma.11_j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F]) * (nu + 1)
        val = mvtnorm::pmvt(lower=rep(0,n-1),upper=upper,delta=sigma[-j,j],sigma=sigma.11_j,df=nu+1)[1]
        return(log(val))
    }
    T_j = unlist(lapply(1:n,a_fun,upper=rep(Inf,n-1)))
    a_j = T_j-logphi+log(2)*((nu-2)/2)+gamma_1-1/2*log(pi)
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
    sigma = par[[2]];nu = par[[1]]
    n = nrow(sigma)
    phi = mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)
    gamma_1 = gamma((nu+1)/2)
    sigma_fun <- function(j){
        sigma_j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F])/ (nu + 1)
    }
    a_fun <- function(j,upper,sigma_j){
        val = mvtnorm::pmvt(lower=rep(0,n-1),upper=upper,delta=sigma[-j,j],sigma=sigma_j)
        return(val)
    }
    sigma_j = lapply(1:n,sigma_fun)
    T_j = mapply(a_fun,sigma_j=sigma_j,j=1:n,MoreArgs = list(upper=rep(Inf,n-1)),SIMPLIFY = TRUE)
    a_j = T_j/phi*2^((nu-2)/2)*gamma_1*pi^(-1/2)
    if(!is.matrix(x)){x <- matrix(x,nrow=n,byrow=FALSE)}

    func <- function(idx){
        x_j = x[,idx] * a_j
        x_upper = lapply(1:n,function(i){ if(is.finite(x_j[i])) (x_j[-i]/x_j[i])^{1/nu} else rep(0,length(x_j[-i])) })
        V_j = mapply(a_fun,j=1:n,upper=x_upper,sigma_j = sigma_j,SIMPLIFY = TRUE)
        val = sum(V_j/T_j/x[,idx])
    }
    if(parallel){
        val = unlist(parallel::mclapply(1:ncol(x),func,mc.cores = ncores))
    }
    else{
        val = unlist(lapply(1:ncol(x),func))
    }
    return(val)
}


## this function returns the partial derivatives of the exponent function
## for the truncated extremal-t max-stable processes
partialV_truncT <- function(x,idx,par,parallel=TRUE,ncores=2,log=TRUE){
    sigma = par[[2]];nu = par[[1]]
    n = nrow(sigma)
    phi = mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)
    chol.sigma.11 = chol(sigma[idx,idx])
    inv.sigma.11 = chol2inv(chol.sigma)
    logdet.sigma.11 = sum(log(diag(chol.sigma)))*2
    sigma_T = sigma[-idx,-idx] - sigma[-idx,idx,drop=F] %*% inv.sigma %*% sigma[-idx,idx,drop=F] 
    k = length(idx)
    gamma_1 = gamma((nu+1)/2)
    gamma_k = gamma((k+nu)/2)
    a_fun <- function(j,upper){
        sigma_j = (sigma[-j,-j] - sigma[-j,j] %*% sigma[j,-j]) * (nu + 1)
        val = mvtnorm::pmvt(lower=rep(0,n-1),upper=upper,delta=sigma[-j,j],sigma=sigma_j)
        return(val)
    }
    T_j = unlist(lapply(1:n,FUN=a_fun,upper=rep(Inf,n-1)))
    a_j = T_j/phi*2^((nu-2)/2)*gamma_1*pi^(-1/2)

    if(!is.matrix(x)){x <- matrix(x,nrow=n,byrow=FALSE)}
    
    func <- function(idx_j){
        x_j = x[,idx_j] * a_j
        x_log = log(x_j)
        Q_sigma = (t(x_j) %*% inv.sigma %*% x_j)^(1/2)
        val = 1/nu*sum(x_log[idx]-x[idx,idx_j,drop=F]) - phi + (nu-2)/2 * log(2) + (1-k)*log(nu) -
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
    b = c((t(alpha) %*% omega_inv %*% rep(1,n))^2/sum_inv_sigma)
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
    delta = c(sigma_bar %*% alpha/c(1+t(alpha) %*% sigma_bar %*% alpha))
    a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)
    if(!is.matrix(x)){x <- matrix(x,nrow=n,byrow=FALSE)}
    b = t(alpha) %*% omega_inv %*% sigma
    I.mat1 = diag(rep(1,n))
    I.mat2 = diag(rep(1,n-1))
    func <- function(j){        
        if(j<n) A.j = t(cbind(I.mat2[,0:(j-1)],rep(-1,n-1),I.mat2[,j:(n-1)])) 
        else A.j = t(cbind(I.mat2[,0:(j-1)],rep(-1,n-1))) 
        sigma.j = t(A.j) %*% sigma %*% A.j
        sigma.j.chol = chol(sigma.j)
        sigma.j.inv = chol2inv(sigma.j.chol)
        omega.j = sqrt(diag(diag(sigma.j),nrow=ncol(A.j)))
        omega.j.inv = diag(diag(omega.j)^(-1),nrow=ncol(A.j))
        sigma.j.bar = omega.j.inv %*% sigma.j %*% omega.j.inv
        alpha.hat = (1 - t(delta) %*% omega %*% A.j %*% sigma.j.inv %*% t(A.j) %*% omega %*% delta)^(-1/2) %*% omega.j %*% 
            sigma.j.inv %*% t(A.j) %*% omega %*% delta   
        u.j = t(A.j) %*% sigma %*% I.mat1[,j]
        b1 = c(t(alpha.hat) %*% sigma.j.bar %*% alpha.hat)
        b2 = c(t(alpha) %*% omega_inv %*% sigma %*% I.mat1[,j]*(1+b1)^(-1/2))
        b3 = c(-(1+b1)^(-1/2)*sigma.j.bar %*% alpha.hat)
        sigma_circ = unname(cbind(rbind(sigma.j.bar,b3),rbind(b3,1)))
        func_temp <- function(i){
            xi = x[,i]
            mu = c(rbind(omega.j.inv %*% (a[-j] - a[j] + log(xi[-j]/xi[j])-u.j),b2))
            val = pnorm(b2)^(-1)/xi[j] * mvtnorm::pmvnorm(lower=rep(-Inf,n),upper=mu,sigma=sigma_circ)[1]
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
#TODO: equation 46, 50 and page 130 eqaution 5.28
partialV_logskew <- function(x,idx,par,parallel=TRUE,ncores=2,log=TRUE){
    alpha = par[[1]];sigma = par[[2]]
    n = nrow(sigma)
    sigma.chol = chol(sigma)
    sigma.inv = chol2inv(sigma.chol)
    sum.sigma.inv = sum(sigma.inv)
    omega = diag(sqrt(diag(sigma)))
    omega.inv = diag(diag(omega)^(-1))
    sigma.bar = omega.inv %*% sigma %*% omega.inv
    ones <- rep(1,n)
    
    sigma.bar.11.inv = chol2inv(chol(sigma.bar[idx,idx]))
    sigma.2.1.bar = sigma.bar[-idx,-idx] - sigma.bar[-idx,idx,drop=F] %*% sigma.bar.11.inv %*% sigma.bar[idx,-idx,drop=F]
    sigma.tilde.inv = sigma.inv %*% ones %*% t(ones) %*% sigma.inv/sum.sigma.inv - sigma.inv
    sigma.tilde.inv.chol = chol(sigma.tilde.inv)
    sigma.tilde = chol2inv(chol(sigma.tilde.inv.chol))
    omega.tilde = diag(sqrt(diag(sigma.tilde)))
    omega.tilde.inv = diag(diag(omega.tilde)^(-1))
    sigma.tilde.bar = omega.tilde.inv %*% sigma.tilde %*% omega.tilde.inv
    sigma.tilde.11 = sigma.tilde[idx,idx]
    sigma.tilde.11.chol = chol(sigma.tilde.11)
    sigma.tilde.11.inv = chol2inv(sigma.tilde.11.chol)
    sigma.tilde.11.bar = sigma.tilde.bar[idx,idx]
    sigma.tilde.11.bar.inv = chol2inv(chol(sigma.tilde.11.bar))
    sigma.tilde.2.1.bar = sigma.tilde.bar[-idx,-idx] - sigma.tilde.bar[-idx,idx,drop=F] %*% sigma.tilde.11.bar.inv %*% sigma.tilde.bar[idx,-idx,drop=F]
    sigma.tilde.2.1 = sigma.tilde[-idx,-idx] - sigma.tilde[-idx,idx,drop=F] %*% sigma.tilde.11.inv %*% sigma.tilde[idx,-idx,drop=F]
    omega.tilde.2.1 = diag(sqrt(diag(sigma.tilde.2.1)))
    omega.tilde.2.1.inv = diag(diag(omega.tilde.2.1)^(-1))
    sigma.tilde.2.1.bar.true = omega.tilde.2.1.inv %*% sigma.tilde.2.1.bar %*% omega.tilde.2.1.inv
    mu.tilde = sigma.tilde.inv %*% sigma.inv %*% ones /sum.sigma.inv/2
    
    delta = sigma_bar %*% alpha/(1+t(alpha) %*% sigma_bar %*% alpha)
    a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)

    b = c((t(alpha) %*% omega.inv %*% ones)^2/sum.sigma.inv)
    alpha.tilde =  t(alpha) %*% omega %*% (diag(ones) + ones %*% t(colSums(sigma.inv))) * (1+b)^(-1/2)
    delta.hat = (1+b)^(-1/2)*sqrt(b)

    alpha.tilde.0 = delta.hat/sqrt(sum.sigma.inv) + t(alpha.tilde) %*% mu.tilde

    tau.tilde = alpha.tilde.0 * (1 + alpha.tilde %*% sigma.tilde.bar %*% alpha.tilde)^(-1/2)
    
    

    alpha2.1 = omega.tilde.2.1 %*% omega.tilde.inv[-idx,-idx] %*% alpha.tilde[-idx]

    alpha1.2 = (1 + t(alpha.tilde[-idx]) %*% sigma.tilde.2.1.bar %*% alpha.tilde[-idx])^(-1/2) * (alpha.tilde[idx] - sigma.tilde.11.bar.inv %*% 
                sigma.tilde.bar[idx,-idx] %*% alpha.tilde[-idx] )
    alpha.k = (1 + t(alpha[-idx]) %*%  sigma.2.1.bar %*% alpha[-idx])^(-1/2) * (alpha[idx] - sigma.11.bar.inv %*% sigma.bar[idx,-idx,drop=F] %*% alpha[-idx] )

    b1 = t(alpha2.1) %*% sigma.tilde.2.1.bar.true %*% alpha2.1
    b2 = (1+ b1)^(-1/2) * sigma.tilde.2.1.bar.true %*% alpha2.1
    func <- function(i){
        xi = x[,i]
        tau2.1 = tau.tilde * (1 + t(alpha1.2) %*% sigma.tilde.11.bar %*% alpha1.2) + alpha1.2 %*% omega.tilde.inv[idx,idx] %*% (log(xi[idx]) + a[idx] - mu.tilde[idx])
        mu.2.1 = c(mu.tilde[-idx] - sigma.tilde[-idx,idx,drop=F] %*% sigma.tilde.11.inv %*% (log(x[idx]) + a[idx] - mu.tilde[idx]) ) 
        mu.val = rbind(omega2.1.inv %*% (log(xi[-idx]) + a[-idx] - mu.2.1),tau2.1 * (1+b1)^(-1/2))
        scale.val = cbind(rbind(sigma.2.1.bar, b2),rbind(b2,1))
        phi = pnorm(b2)
        intensity.marginal = intensity_logskew(xi[idx],par=list(alpha=alpha.k,sigma=sigma[idx,idx]),parallel=FALSE,log=FALSE)
        val = intensity.marginal * phi^(-1) * mvtnorm::pmvnorm(lower=rep(-Inf,length(mu.val)),upper=mu.val,sigma=scale.val)[1] 
        return(val)
    }
    if(parallel){
        val = unlist(parallel::mclapply(1:ncol(x),func,mc.cores = ncores))
    }
    else{
        val = unlist(lapply(1:ncol(x),func))
    }
    if(log){
        return(log(val))
    }
    else{
        return(val)
    } 
}

# calculate empirical extremal coefficients
empirical_extcoef <- function(p,data){
    return(min(2,max(1,1/mean(1/pmax(data[,p[1]],data[,p[2]])))))
}

true_extcoef <- function(idx,par,model="logskew1"){
    if(model=="logskew1"){
        alpha = par[[1]];sigma = par[[2]]
        n = nrow(sigma)
        omega = diag(sqrt(diag(sigma)))
        omega.inv = diag(diag(omega)^(-1))
        sigma.bar = omega.inv %*% sigma %*% omega.inv
        sigma.bar.11 = sigma.bar[idx,idx]
        sigma.bar.11.chol = chol(sigma.bar.11)
        sigma.bar.11.inv = chol2inv(sigma.bar.11.chol)
        sigma.22.1 = sigma.bar[-idx,-idx] - sigma.bar[-idx,idx] %*% sigma.bar.11.inv %*% sigma.bar[idx,-idx]
        alpha.new = c(c(1 + t(alpha[-idx]) %*% sigma.22.1 %*% alpha[-idx])^(-1/2) * (alpha[idx] - sigma.bar.11.inv %*% sigma.bar[idx,-idx] %*% alpha[-idx]))
        sigma.new = sigma[idx,idx]
        val = V_logskew(rep(1,length(idx)),list(alpha=alpha.new,sigma=sigma.new),parallel=FALSE)
    }

    if(model == "logskew2"){
        alpha = par[[1]];sigma = par[[2]]
        val = V_logskew(rep(1,length(idx)),list(alpha=alpha[idx],sigma=sigma[idx,idx]),parallel=FALSE)
    }
    
    if(model == "truncT1"){
        x = matrix(Inf,nrow=nrow(par[[2]]),ncol=ncol(all.pairs))
        all.pairs.new = cbind(c(all.pairs),rep(1:ncol(all.pairs),each=nrow(all.pairs)))
        x[all.pairs.new] = 1
        val = V_truncT(x,par,parallel=TRUE,ncores=parallel::detectCores())
    }

    if(model == "truncT2"){
        x = rep(1,length(idx))
        val = V_truncT(x,list(nu = par[[1]],sigma=par[[2]][idx,idx]),parallel=FALSE)
    }
    return(val)

}