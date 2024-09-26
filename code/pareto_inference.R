## weight functions that used in the gradient scoring method
weightFun <- function(x){
    val = x*(1- exp(1-rFun(x)))
    return(val)
}

dWeightFun <- function(x){
    r = rFun(x)
    val = (1 - exp(-(r - 1))) + x * exp(-(r - 1))*r^(1-xi)*(x)^(xi-1)
    return(val)
}

### score matching inference for the skewed Brown-Resnick and Truncated extremal-t Process ###
################################################################################################
scoreMatching <- function (par2, obs, loc, model="logskew", vario.func=NULL,cov.func=NULL, alpha.func=NULL,basis=NULL,idx.para=1:2, dof=2, weightFun = NULL, dWeightFun = NULL, nCores = 1L, ...){
    ellipsis <- list(...)
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    if ("weigthFun" %in% names(ellipsis) && is.null(weightFun)) {
        weightFun <- ellipsis[["weigthFun"]]
        ellipsis[["weigthFun"]] <- NULL
    }
    if ("dWeigthFun" %in% names(ellipsis) && is.null(dWeightFun)) {
        dWeightFun <- ellipsis[["dWeigthFun"]]
        ellipsis[["dWeigthFun"]] <- NULL
    }
    if (is.null(weightFun)) {
        stop("`weightFun` argument missing, with no default value.")
    }
    if (is.null(dWeightFun)) {
        stop("`dWeightFun` argument missing, with no default value.")
    }
    if (!inherits(obs, "list") || length(obs) < 1 || !inherits(obs[[1]], 
        c("numeric", "integer"))) {
        stop("`obs` must be a list of vectors.")
    }
    if (!inherits(weightFun, "function")) {
        stop("`weightFun` must be a function.")
    }
    if (!inherits(dWeightFun, "function")) {
        stop("`dWeightFun` must be a function.")
    }
    if (!is.numeric(nCores) || nCores < 1) {
        stop("`nCores` must a positive number of cores.")
    }
    if(model=="logskew"){
        n <- length(obs)
        SigmaS = vario.fun(loc, par2[idx.para])
        alpha = alpha.func(par2[-idx.para],b.mat=basis)
        delta = c(SigmaS %*% alpha)/sqrt(c(1+alpha %*% SigmaS %*% alpha))
        a = log(2) + pnorm(delta,log.p=TRUE)
        computeScores= function(i){
            obs.i = .subset2(obs, i)
            ind = !is.na(obs.i)
            obs.i = obs.i[ind]
            log.obs.i = log(obs.i) + a
            sigmaInv <- MASS::ginv(SigmaS[ind, ind])
            sigma <- diag(SigmaS[ind, ind])
            d = nrow(sigmaInv)
            alpha.i = c(1 - delta[ind] %*% sigmaInv %*% delta[ind])^(-1/2) * c(sigmaInv %*% delta[ind])
            q = rowSums(sigmaInv)
            sum.q = sum(q);sum.alpha = sum(alpha.i)
            q.mat = matrix(q,d,d,byrow=TRUE)

            beta = (1+sum.alpha^2/sum.q)^(-0.5)
            tau.tilde.beta = beta * (alpha.i - sum.alpha*q/sum.q)
            tau.tilde =  sum( tau.tilde.beta * (log.obs.i + sigma/2) ) + beta*sum.alpha/sum.q
            Phi.tau = pnorm(tau.tilde);phi.tau = dnorm(tau.tilde)
            A = sigmaInv - q %*% t(q)/sum.q
            mtp <- 2 * q/(sum(q)) + 2 + A %*% sigma
            gradient <- - A %*% log.obs.i * (1/obs.i) - 
                    1/2 * (1/obs.i) * mtp +  phi.tau/Phi.tau*tau.tilde.beta*(1/obs.i) 

            diagHessian <- -diag(A) * (1/obs.i^2) + 
                A %*% log.obs.i * (1/obs.i)^2 + 
                1/2 * (1/obs.i)^2 * mtp + (-phi.tau*tau.tilde*tau.tilde.beta^2+phi.tau*tau.tilde.beta-phi.tau^2*tau.tilde.beta^2/Phi.tau)/(obs.i^2*Phi.tau)

            weights <- do.call(what = "weightFun", args = c(ellipsis, 
                x = list(obs.i)))
            dWeights <- do.call(what = "dWeightFun", args = c(ellipsis, 
                x = list(obs.i)))

            sum(2 * (weights * dWeights) * gradient + weights^2 * 
                diagHessian + 1/2 * weights^2 * gradient^2)
        }
    }
    if(model=="truncT"){
        n <- length(obs)
        SigmaS = cov.func(loc, par2[idx.para])
        sigmaInv <- MASS::ginv(SigmaS)
        a_fun <- function(j,upper=rep(Inf,n-1)){
            sigma_j = (SigmaS[-j,-j] - SigmaS[-j,j,drop=F] %*% SigmaS[j,-j,drop=F]/SigmaS[j,j])/(dof + 1)/SigmaS[j,j]
            val = mvtnorm::pmvt(lower=-SigmaS[-j,j]/SigmaS[j,j],upper=rep(Inf,n-1),sigma=sigma_j,df=dof+1)[[1]]
            return(log(val))
        }
        
        logphi = log(mvtnorm::pmvnorm(lower=rep(0,n),upper=rep(Inf,n),sigma=SigmaS)[[1]])
        if(!is.null(ncores)) T_j = unlist(mclapply(1:n,a_fun,mc.cores=ncores,mc.set.seed = FALSE)) else T_j = unlist(lapply(1:n,a_fun))
        a = T_j - logphi+ (dof-2)/2 * log(2) + log(gamma((dof+1)/2)) - 1/2*log(pi)
        # a = par2[-idx.para]
        dim = nrow(SigmaS)
        if(any(is.na(unlist(obs)))){ stop("Data must be complete.") }
        computeScores <- function(i){
            obs.i = .subset2(obs, i)
            obs.i.a = (obs.i*a)^(1/dof)
            obs.i.a.quad = sum(t(obs.i.a) %*% sigmaInv %*% obs.i.a)
            
            gradient <- (1-dof)/dof/obs.i - (dof+dim)/dof/obs.i.a.quad*sigmaInv %*% obs.i.a * (obs.i.a^(1-dof))*a 

            diagHessian <- (dof-1)/dof/obs.i^2 - 
                a^2*(dof+dim)/dof/obs.i.a.quad*( 1/dof*diag(sigmaInv)*obs.i.a^(2-2*dof) + (1-dof)/dof*sigmaInv %*% obs.i.a * (obs.i.a^(1-2*dof)) ) + 
                2*a^2*(dof+dim)/(dof^2)/(obs.i.a.quad^2)*(sigmaInv %*% obs.i.a * (obs.i.a^(1-dof)))
            
            weights <- do.call(what = "weightFun", args = c(ellipsis, 
                x = list(obs.i)))
            dWeights <- do.call(what = "dWeightFun", args = c(ellipsis, 
                x = list(obs.i)))

            sum(2 * (weights * dWeights) * gradient + weights^2 * 
                diagHessian + 1/2 * weights^2 * gradient^2)
        }
    }
    if (nCores > 1) {
        scores <- parallel::mclapply(1:n, computeScores, mc.cores = nCores)
    }
    else {
        scores <- lapply(1:n, computeScores)
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
    return(sum(unlist(scores))/n)
}

fit.scoreMatching <- function(init, obs, loc,fixed=c(F,F,F,F,F), model="logskew", vario.func=NULL,cov.func=NULL,basis=NULL,thres=1, idx.para=1:2, alpha.func=NULL, dof=2, weightFun = NULL, dWeightFun = NULL, method="Nelder-Mead", maxit=1000, nCores = 1L,lb,ub, ...){
    t1 = proc.time()
    if (is.matrix(obs)) {
        obs <- split(obs, row(obs))
    }
    if (is.null(weightFun) | is.null(dWeightFun)) {
        weightFun <- function(x){
            x * (1 - exp(-(mean(x) - 1)))
        }
        dWeightFun <- function(x) {
            (1 - exp(-(mean(x) - 1))) + (x/length(x)) * exp(-(mean(x) - 
                1))
        }
    }
    fun <- function(par){
        par2 = init
        par2[!fixed] = par
        if(any(par < lb[!fixed]) | any(par > ub[!fixed])){return(Inf)}
        val = scoreMatching(par2, obs, loc, model, vario.func, cov.func, alpha.func, basis, idx.para, dof, weightFun=weightFun, dWeightFun=dWeightFun,  nCores, ...)
        return(val)
    }
    init2 = init[!fixed]
    result = optim(init2, fun, control = list(trace = TRUE, maxit = maxit), method = method, hessian = FALSE)
    # val.old = fun(init2)
    # n = sum(!fixed);m=length(init)  
    # fixed.old = fixed
    # fixed.idx = which(!fixed.old) 
    # ub[!is.finite(ub)] <- 1e+5
    # for(i in 1:maxit){
    #     for(j in 1:n){
    #         ind =  (j:j+2) %% n + 1
    #         fixed <- rep(T,m);fixed[fixed.idx[ind]] = F 
    #         result <- optim(init[fixed.idx[ind]], fun, control = list(trace = TRUE, maxit = 100), method = "L-BFGS-B", hessian = FALSE,lower=lb[fixed.idx][ind],upper=ub[fixed.idx][ind])
    #         init[fixed.idx[ind]] <- result$par
    #         print(paste0(i,": ",result$value,": ",result$par))
    #     }
    #     if(abs(result$value - val.old)<0.001){
    #         break
    #     }else{
    #         val.old = result$value
    #     }
    # }    
    t2 = proc.time() - t1
    init[!fixed] = result$par
    result$par = init
    result$time = t2
    return(result)
}



