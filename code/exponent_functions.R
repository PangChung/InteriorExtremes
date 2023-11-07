#########################################################
###### Intensity function for truncated extremal-t ######
#########################################################

## this function returns the intensity function of the
## truncated extremal-t max-stable processes
intensity_truncT <- function(x,par,ncores=NULL,log=TRUE){
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    sigma = par[[1]];nu = par[[2]]
    n = ncol(x)
    if(n==1) return(1/(x^2))
    chol.sigma = chol(sigma)
    inv.sigma = chol2inv(chol.sigma)
    logdet.sigma = sum(log(diag(chol.sigma)))*2
    logphi = log(mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[1])
    gamma_1 = log(gamma((nu+1)/2))
    gamma_n = log(gamma((nu+d)/2))
    a_fun <- function(j,upper){
        sigma.11_j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F])/(nu + 1)
        val = mvtnorm::pmvt(lower=rep(0,n-1)-sigma[-j,j],upper=upper-sigma[-j,j],sigma=sigma.11_j,df=nu+1)[1]
        return(log(val))
    }
    if(!is.null(ncores)) T_j = unlist(mclapply(1:n,a_fun,upper=rep(Inf,n-1),mc.cores=ncores)) else T_j = unlist(lapply(1:n,a_fun,upper=rep(Inf,n-1)))
    a_j = T_j - logphi+log(2)*((nu-2)/2)+gamma_1-1/2*log(pi)
    if(!is.matrix(x)){ x = matrix(x,ncol=n,byrow=TRUE)} 
    func <- function(idx){
        log_x = log(x[idx,])
        x_j = (log_x + a_j) * 1/nu
        x_circ = exp(x_j)
        val = -logphi + sum(x_j) - sum(log_x) - (n-1)*log(nu) - logdet.sigma/2 - n/2 * log(2*pi) +  log(t(x_circ) %*% inv.sigma %*% x_circ) * 
            (-nu-n)/2 + gamma_n
        return(val)
    }
    if(!is.null(ncores)){
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))

    }else{
        val = unlist(lapply(1:nrow(x),func))
    }
    if(log) return(val)
    else return(exp(val))
} 

## this function returns the exponent function of the
## truncated extremal-t max-stable processes
V_truncT <- function(x,par,ncores=NULL){
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    sigma = par[[1]];nu = par[[2]]
    n = ncol(x)
    phi = mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[[1]]
    gamma_1 = gamma((nu+1)/2)
    sigma_fun <- function(j){
        sigma_j = (sigma[-j,-j] - sigma[-j,j,drop=F] %*% sigma[j,-j,drop=F])/ (nu + 1)
    }
    a_fun <- function(j,upper,sigma_j){
        val = mvtnorm::pmvt(lower=rep(0,n-1)-sigma[-j,j],upper=upper-sigma[-j,j],sigma=sigma_j,df=nu+1)[[1]]
        return(val)
    }
    sigma_j = lapply(1:n,sigma_fun)
    idx.finite <- which(apply(x,2,function(xi){any(is.finite(xi))}))
    T_j = rep(1,n)
    if(!is.null(ncores)) T_j[idx.finite] = mcmapply(a_fun,sigma_j=sigma_j[idx.finite],j=idx.finite,MoreArgs = list(upper=rep(Inf,n-1)),SIMPLIFY = TRUE,mc.cores=ncores)
    else T_j[idx.finite] = mapply(a_fun,sigma_j=sigma_j[idx.finite],j=idx.finite,MoreArgs = list(upper=rep(Inf,n-1)),SIMPLIFY = TRUE)
    a_j = rep(1,n)
    a_j[idx.finite] = T_j[idx.finite]/phi*2^((nu-2)/2)*gamma_1*pi^(-1/2)
    
    func <- function(idx){
        x_j = x[idx,] * a_j
        x_upper = lapply(1:n,function(i){ if(is.finite(x_j[i])) (x_j[-i]/x_j[i])^{1/nu} else rep(0,length(x_j[-i])) })
        is.finite.ind = which(is.finite(x_j))
        V_j = rep(0,n)
        V_j[is.finite.ind] = mapply(a_fun,j=is.finite.ind,upper=x_upper[is.finite.ind],sigma_j = sigma_j[is.finite.ind],SIMPLIFY = TRUE)
        val = sum(V_j[is.finite.ind]/T_j[is.finite.ind]/x[idx,is.finite.ind])
    }
    if(!is.null(ncores)){
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))
    }else{
        val = unlist(lapply(1:nrow(x),func))
    }
    return(val)
}

## this function returns the partial derivatives of the exponent function
## for the truncated extremal-t max-stable processes
partialV_truncT <- function(x,idx,par,ncores=NULL,log=TRUE){
    sigma = par[[1]];nu = par[[2]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(length(idx)==0){
        val = V_truncT(x,par,ncores=ncores)
        if(log) return(log(val))
        else return(val)
    }
    if(length(idx)==n){
       val = intensity_truncT(x,par,ncores,log)
       return(val)
    }
    phi = mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),mean=rep(0,n),sigma=sigma)[[1]]
    chol.sigma.11 = chol(sigma[idx,idx])
    inv.sigma.11 = chol2inv(chol.sigma.11)
    logdet.sigma.11 = sum(log(diag(chol.sigma.11)))*2
    sigma_T = sigma[-idx,-idx] - sigma[-idx,idx,drop=F] %*% inv.sigma.11 %*% sigma[-idx,idx,drop=F] 
    k = length(idx)
    gamma_1 = gamma((nu+1)/2)
    gamma_k = gamma((k+nu)/2)
    a_fun <- function(j,upper){
        sigma_j = (sigma[-j,-j] - sigma[-j,j] %*% sigma[j,-j])/(nu + 1)
        val = mvtnorm::pmvt(lower=rep(0,n-1)-sigma[-j,j],upper=upper-sigma[-j,j],sigma=sigma_j,df=nu+1)[[1]]
        return(val)
    }
    if(!is.null(ncores)) T_j = unlist(mclapply(1:n,a_fun,upper=rep(Inf,n-1),mc.cores=ncores)) else T_j = unlist(lapply(1:n,a_fun,upper=rep(Inf,n-1)))
    a_j = T_j/phi*2^((nu-2)/2)*gamma_1*pi^(-1/2)
    func <- function(idx_j){
        x_j = x[idx_j,] * a_j
        x_log = log(x_j)
        Q_sigma = (x_j[idx] %*% inv.sigma.11 %*% x_j[idx])^(1/2)
        val = c(1/nu*sum(x_log[idx]- log(x[idx_j,idx,drop=F])) - phi + (nu-2)/2 * log(2) + (1-k)*log(nu) -
            k/2*log(pi) - 1/2*logdet.sigma.11 -(k+nu)*log(Q_sigma)+log(gamma_k))
        upper = c((x_j[-idx])^(1/nu)*sqrt(k+nu)/Q_sigma)
        loc = c(sigma[-idx,idx] %*% inv.sigma.11 %*% x_j[idx]^(1/nu)*sqrt(k+nu)/Q_sigma)
        val = val + log(mvtnorm::pmvt(lower=rep(0,n-k)-loc,upper=upper-loc,sigma=sigma_T,df=k+nu))[[1]]
        return(val)
    }
    if(!is.null(ncores)){
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))
    }
    else{
        val = unlist(lapply(1:nrow(x),func))
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

# this function computes the intensity function 
# for the log skew-normal based max-stable processes
intensity_logskew <- function(x,par,ncores=NULL,log=TRUE){
    alpha = par[[2]];sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(n==1) return(1/(x^2))
    omega = diag(sqrt(diag(sigma)))
    omega_inv = diag(diag(omega)^(-1))
    sigma_bar = omega_inv %*% sigma %*% omega_inv
    chol.sigma = chol(sigma)
    inv.sigma = chol2inv(chol.sigma)
    logdet.sigma = sum(log(diag(chol.sigma)))*2
    delta = c(sigma_bar %*% alpha)/sqrt(c(1+alpha %*% sigma_bar %*% alpha))
    a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)
    sum.inv.sigma = sum(inv.sigma)
    one.mat = matrix(1,n,n)
    one.vec = rep(1,n)
    b = c((alpha %*% omega_inv %*% one.vec)^2/sum.inv.sigma)
    beta =  c(alpha %*% omega_inv %*% (diag(n) - one.mat %*% inv.sigma/sum.inv.sigma) * (1+b)^(-1/2))
    A = inv.sigma - inv.sigma %*% one.mat %*% inv.sigma/sum.inv.sigma 
    delta.hat = (1+b)^(-1/2)*sqrt(b)
    func <- function(idx){
        x_log = log(x[idx,])
        x_circ = x_log + a
        val = -(n-3)/2 * log(2) - (n-1)/2*log(pi) + pnorm(beta %*% x_circ + delta.hat/sqrt(sum.inv.sigma),log.p=TRUE) -
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
    if(log)
        return(val)
    else
        return(exp(val))    
}

# intensity_logskew <- function(x,par,ncores=NULL,log=TRUE){
#     delta = par[[2]];sigma = par[[1]]
#     if(!is.matrix(x)){x <- matrix(x,nrow=1)}
#     n = ncol(x)
#     if(n==1) return(1/(x^2))
#     omega = diag(sqrt(diag(sigma)))
#     omega_inv = diag(diag(omega)^(-1))
#     sigma_bar = omega_inv %*% sigma %*% omega_inv
#     chol.sigma = chol(sigma)
#     inv.sigma = chol2inv(chol.sigma)
#     inv.sigma.bar = omega %*% inv.sigma %*% omega
#     alpha = c(c(1 - delta %*% inv.sigma.bar %*% delta)^(-1/2) * inv.sigma.bar %*% delta)
#     logdet.sigma = sum(log(diag(chol.sigma)))*2
#     #delta = c(sigma_bar %*% alpha)/sqrt(c(1+alpha %*% sigma_bar %*% alpha))
#     a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)
#     sum.inv.sigma = sum(inv.sigma)
#     one.mat = matrix(1,n,n)
#     one.vec = rep(1,n)
#     b = c((alpha %*% omega_inv %*% one.vec)^2/sum.inv.sigma)
#     beta =  c(alpha %*% omega_inv %*% (diag(n) - one.mat %*% inv.sigma/sum.inv.sigma) * (1+b)^(-1/2))
#     A = inv.sigma - inv.sigma %*% one.mat %*% inv.sigma/sum.inv.sigma 
#     delta.hat = (1+b)^(-1/2)*sqrt(b)
#     func <- function(idx){
#         x_log = log(x[idx,])
#         x_circ = x_log + a
#         val = -(n-3)/2 * log(2) - (n-1)/2*log(pi) + pnorm(beta %*% x_circ + delta.hat/sqrt(sum.inv.sigma),log.p=TRUE) -
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
#     if(log)
#         return(val)
#     else
#         return(exp(val))    
# }

## this function computes the exponent function 
## for the log skew-normal based max-stable processes
V_logskew <- function(x,par,ncores=NULL){
    alpha = par[[2]];sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,nrow=1)}
    n = ncol(x)
    if(n==1) return(1/x)
    omega = diag(sqrt(diag(sigma)))
    omega_inv = diag(diag(omega)^(-1))
    sigma_bar = omega_inv %*% sigma %*% omega_inv
    chol.sigma = chol(sigma)
    inv.sigma = chol2inv(chol.sigma)
    logdet.sigma = sum(log(diag(chol.sigma)))*2
    b = c(alpha %*% sigma_bar %*% alpha)
    delta = c(sigma_bar %*% alpha)/sqrt(c(1+b))
    a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)
    I.mat1 = diag(rep(1,n))
    I.mat2 = diag(rep(1,n-1))
    func <- function(j){        
        if(j<n){
        A.j = cbind(I.mat2[,0:(j-1)],rep(-1,n-1),I.mat2[,j:(n-1)])
        }else{
        A.j = cbind(I.mat2[,0:(j-1)],rep(-1,n-1))
        }
        sigma.j = A.j %*% sigma %*% t(A.j)
        sigma.j.chol = chol(sigma.j)
        sigma.j.inv = chol2inv(sigma.j.chol)
        omega.j = sqrt(diag(diag(sigma.j),nrow=n-1))
        omega.j.inv = diag(diag(omega.j)^(-1),nrow=n-1)
        sigma.j.bar = omega.j.inv %*% sigma.j %*% omega.j.inv
        alpha.hat = (1 - delta %*% omega %*% t(A.j) %*% sigma.j.inv %*% A.j %*% omega %*% delta)^(-1/2) %*% omega.j %*% 
            sigma.j.inv %*% A.j %*% omega %*% delta   
        u.j = A.j %*% sigma %*% I.mat1[,j]
        b1 = c(alpha.hat %*% sigma.j.bar %*% alpha.hat)
        b3 = c(-(1+b1)^(-1/2)*sigma.j.bar %*% alpha.hat)
        sigma_circ = unname(cbind(rbind(sigma.j.bar,b3),rbind(b3,1)))
        func_temp <- function(i){
            xi = x[i,]
            mu = c(omega.j.inv %*% (a[-j] - a[j] + log(xi[-j]/xi[j])-u.j),delta[j]*omega[j,j])
            val = pnorm(delta[j]*omega[j,j])^(-1)/xi[j] * mvtnorm::pmvnorm(lower=rep(-Inf,n),upper=mu,sigma=sigma_circ,algorithm = "TVPACK")[[1]]
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
    return(val)
}

## this function returns the partial derivatives of the exponent function
## for the truncated extremal-t max-stable processes
partialV_logskew <- function(x,idx,par,ncores=NULL,log=FALSE){
    alpha = par[[2]];sigma = par[[1]]
    if(!is.matrix(x)){x <- matrix(x,ncol=1)}
    n = ncol(x)
    if(length(idx)==0){
        val = V_logskew(x,par,ncores=ncores)
        if(log) return(log(val))
        else return(val)
    }
    if(length(idx)==n){
        val = intensity_logskew(x,par,ncores,log)
        return(val)
    }
    ones <- rep(1,n)
    one.mat <- matrix(1,n,n)
    
    sigma.chol = chol(sigma)
    sigma.inv = chol2inv(sigma.chol)
    sum.sigma.inv = sum(sigma.inv)
    omega = diag(sqrt(diag(sigma)))
    omega.inv = diag(diag(omega)^(-1))
    sigma.bar = omega.inv %*% sigma %*% omega.inv
    
    delta = c(sigma.bar %*% alpha)/sqrt(c(1+alpha %*% sigma.bar %*% alpha))
    a = log(2) + diag(sigma)/2 + sapply(diag(omega)*delta,pnorm,log.p=TRUE)
    H =  sigma.inv - (sigma.inv %*% one.mat %*% sigma.inv/sum.sigma.inv)

    sigma.tilde.inv = H[-idx,-idx,drop=FALSE]
    sigma.tilde.inv.chol = chol(sigma.tilde)
    sigma.tilde = chol2inv(sigma.tilde.inv.chol)
    omega.tilde = diag(sqrt(diag(sigma.tilde)))
    omega.tilde.inv = diag(diag(omega.tilde)^(-1))  
    sigma.tilde.bar = omega.tilde.inv %*% sigma.tilde %*% omega.tilde.inv

    b = c((alpha %*% omega.inv %*% ones)^2/sum.sigma.inv)
    beta =  c(alpha %*% omega.inv %*% (diag(n) - one.mat %*% sigma.inv/sum.sigma.inv) * (1+b)^(-1/2))    
    delta.hat = (1+b)^(-1/2)*sqrt(b)
    alpha.tilde = c(beta[-idx] %*% omega.tilde) 
    b1 =c((1 + alpha.tilde %*% sigma.tilde.bar %*% alpha.tilde)^(-1/2))
    func <- function(i){
        xi = x[i,]
        xi.tilde = log(xi) + a
        mu.tilde = -sigma.tilde %*% (H[-idx,idx] %*% xi.tilde[idx] + (sigma.inv %*% ones)[-idx]/sum.sigma.inv)
        tau.tilde = b1 * (beta[idx] %*% xi.tilde[idx] + delta.hat * sum.sigma.inv^(-1/2) + beta[-idx] %*% mu.tilde)
        b2 = -b1 %*% sigma.tilde.bar %*% alpha.tilde
        scale.val = cbind(rbind(sigma.2.1.bar, b2),c(b2,1))
        mu.val = c(omega.tilde.inv %*% xi.tilde[-idx] - mu.tilde, tau.tilde * b1)
        rownames(scale.val) <- colnames(scale.val) <-  NULL
        phi = pnorm(c(tau.tilde * b1))
        intensity.marginal = intensity_logskew(xi[idx],par=list(alpha=alpha.k,sigma=sigma[idx,idx]),ncores=NULL,log=FALSE)
        val = intensity.marginal/phi * mvtnorm::pmvnorm(lower=rep(-Inf,length(mu.tilde)),upper=mu.tilde,sigma=scale.val)[[1]]
        return(val)
    }
    if(!is.null(ncores)){
        val = unlist(parallel::mclapply(1:nrow(x),func,mc.cores = ncores))
    }
    else{
        val = unlist(lapply(1:nrow(x),func))
    }
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

# calculate true extremal coefficients
true_extcoef <- function(idx,par,model="logskew1"){
    if(model=="logskew1"){
        alpha = par[[2]];sigma = par[[1]]
        n = nrow(sigma)
        omega = diag(sqrt(diag(sigma)))
        omega.inv = diag(diag(omega)^(-1))
        sigma.bar = omega.inv %*% sigma %*% omega.inv
        sigma.bar.11 = sigma.bar[idx,idx]
        sigma.bar.11.chol = chol(sigma.bar.11)
        sigma.bar.11.inv = chol2inv(sigma.bar.11.chol)
        sigma.22.1 = sigma.bar[-idx,-idx] - sigma.bar[-idx,idx] %*% sigma.bar.11.inv %*% sigma.bar[idx,-idx]
        alpha.new = c(c(1 + alpha[-idx] %*% sigma.22.1 %*% alpha[-idx])^(-1/2) * (alpha[idx] - sigma.bar.11.inv %*% sigma.bar[idx,-idx] %*% alpha[-idx]))
        sigma.new = sigma[idx,idx]
        val = V_logskew(rep(1,length(idx)),list(sigma=sigma.new,alpha=alpha.new),ncores=NULL)
    }

    if(model == "logskew2"){
        alpha = par[[2]];sigma = par[[1]]
        val = V_logskew(rep(1,length(idx)),list(sigma=sigma[idx,idx],alpha=alpha[idx]),ncores=NULL)
    }
    if(model == "truncT1"){
        x = matrix(Inf,nrow=ncol(idx),ncol=nrow(par[[1]]))
        all.pairs.new = cbind(rep(1:ncol(idx),each=nrow(idx)),c(idx))
        x[all.pairs.new] = 1
        val = V_truncT(x,par,ncores=parallel::detectCores())
    }
    
    if(model == "truncT2"){
        x = rep(1,length(idx))
        val = V_truncT(x,list(sigma=par[[1]][idx,idx],nu = par[[2]]),ncores=NULL)
    }
    return(val)

}


cov.func <- function(loc,par){
    r = par[1];v = par[2]
    n = nrow(loc)
    diff.vector <- cbind(as.vector(outer(loc[,1],loc[,1],'-')),
        as.vector(outer(loc[,2],loc[,2],'-')))
    cov.mat <- matrix(exp(-(sqrt(diff.vector[,1]^2 + diff.vector[,2]^2)/r)^v), ncol=n) #+ diag(1e-6,n) 
    return(cov.mat)
}

alpha.func <- function(loc,par){
    beta.1 = par[1];beta.2 = par[2];beta.3 = par[3]
    alpha = beta.1 + loc[,1] * beta.2 + loc[,2] * beta.3
    return(alpha)
}

alpha.func <- function(coord,par=10){
    #alpha = 1 + 1.5*coord[,2] - par * exp(2*sin(2*coord[,2]))
    #alpha = par*exp(2*sin(2*coord[,2]))
    alpha = rep(par,nrow(coord))
}

## inference for simulated data ##  
fit.model <- function(data,loc,init,fixed,thres = 0.90,model="truncT",maxit=100,
                    ncores=NULL,method="Nelder-Mead",lb=NULL,ub=NULL,hessian=FALSE,bootstrap=FALSE,opt=FALSE){
    data.sum = apply(data,1,sum)
    idx.thres = which(data.sum>quantile(data.sum,thres))
    data = sweep(data[idx.thres,],1,data.sum[idx.thres],"/")
    #data = data[idx.thres,]
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    if(model == "logskew"){
    ## 5 parameters: 2 for the covariance function; 3 for the slant parameter
        object.func <- function(par,opt=TRUE,dat=data,ncore=ncores){
            par2 = init; par2[!fixed] = par
            par.1 = par2[1:2];par.2 = par2[-c(1:2)]
            cov.mat = cov.func(loc,par.1)
            alpha = alpha.func(loc,par.2)
            if(any(par < lb[!fixed]) | any(par > ub[!fixed])){return(Inf)}
            #print(par)
            para.temp = list(sigma=cov.mat,alpha=alpha)
            val = - intensity_logskew(dat,par=para.temp,log=TRUE,ncores=ncore)
            if(opt) return(mean(val)) else return(val)
        }
    }
    if(model == "truncT"){
    ## 3 parameters: 2 for the covariance function; 1 for the df parameter
        object.func <- function(par,opt=TRUE,dat=data,ncore=ncores){
            par2 = init; par2[!fixed] = par
            par.1 = par2[1:2];nu = par2[3]
            if(any(par < lb[!fixed]) | any(par > ub[!fixed])){return(Inf)}
            cov.mat = cov.func(loc,par.1)
            para.temp = list(sigma=cov.mat,nu=nu)
            val = -intensity_truncT(dat,par=para.temp,log=TRUE,ncores=ncore)
            if(opt) return(mean(val)) else return(val)
        }
    }
    if(opt){
        opt.result = optim(init[!fixed],object.func,method=method,control=list(maxit=maxit,trace=TRUE),hessian=hessian)
    }else{
        return(object.func(init[!fixed],opt=TRUE))
    }
    if(hessian){
        h = 1e-4
        par.mat.grad = matrix(opt.result$par,nrow=length(opt.result$par),ncol=length(opt.result$par),byrow=TRUE) + diag(h,length(opt.result$par))
        val.object = object.func(opt.result$par,opt=FALSE)
        val.object.grad = apply(par.mat.grad,1,function(x){(object.func(x,opt=FALSE) - val.object)/h})
        opt.result$K = var(val.object.grad)
        opt.result$hessian.inv = solve(opt.result$hessian)
        opt.result$sigma = opt.result$hessian.inv %*% opt.result$K %*% opt.result$hessian.inv
    }
    if(bootstrap){
        func <- function(id){
            ind = sample(1:ncol(data),ncol(data),replace=TRUE)
            data.boot = data[ind,]
            opt.boot = optim(opt.result$par,object.func,dat=data.boot,ncore=NULL,method=method,control=list(maxit=maxit,trace=FALSE),hessian=FALSE)
            val = opt.boot$par
            return(val)
        }
        boot.results <- matrix(unlist(mclapply(1:300,func,mc.cores=ncores,mc.preschedule=TRUE)),nrow=300,byrow=TRUE)
        opt.result$boot = boot.results
    }
    par2 = init; par2[!fixed] = opt.result$par
    opt.result$par = par2
    assign(".Random.seed", oldSeed, envir=globalenv())
    return(opt.result)
}

