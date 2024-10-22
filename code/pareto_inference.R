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
scoreMatching <- function (par2, obs, loc, model="logskew", vario.func=NULL,cov.func=NULL, alpha.func=NULL,basis=NULL,idx.para=1:2, dof=2, weightFun = NULL, dWeightFun = NULL, ncores = NULL,partial=FALSE, ...){
    ellipsis <- list(...)
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
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
    if(model=="logskew"){
        n <- nrow(loc)
        SigmaS = vario.func(loc, par2[idx.para])
        alpha = alpha.func(par2[-idx.para],b.mat=basis)
        delta = c(SigmaS %*% alpha)/sqrt(c(1+alpha %*% SigmaS %*% alpha))
        a = log(2) + pnorm(delta,log.p=TRUE)
        if(!partial){sigmaInv <- chol2inv(chol(SigmaS))}
        computeScores= function(i){
            if(partial){
                obs.i = .subset2(obs,i)
                ind = obs.i[[1]]
                obs.i = obs.i[[2]]
                sigmaInv <- chol2inv(chol(SigmaS[ind, ind]))
            }else {
                obs.i = .subset2(obs, i)
                ind = 1:length(obs.i)
            }
            log.obs.i = log(obs.i) + a[ind]
            sigma <- diag(SigmaS[ind, ind])
            d = nrow(sigmaInv)
            alpha.i = c(1 - delta[ind] %*% sigmaInv %*% delta[ind])^(-1/2) * c(sigmaInv %*% delta[ind])
            q = rowSums(sigmaInv)
            sum.q = sum(q)
            sum.alpha = sum(alpha.i)
            q.mat = matrix(q,d,d,byrow=TRUE)

            beta = (1+sum.alpha^2/sum.q)^(-0.5)
            tau.tilde.beta = beta * (alpha.i - sum.alpha*q/sum.q)
            tau.tilde =  sum( tau.tilde.beta * (log.obs.i + sigma/2) ) + beta*sum.alpha/sum.q
            Phi.tau = pnorm(tau.tilde,log=TRUE);phi.tau = dnorm(tau.tilde,log=TRUE)
            phi.diff = exp(phi.tau-Phi.tau)
            A = sigmaInv - q %*% t(q)/sum.q
            mtp <- 2 * q/(sum(q)) + 2 + A %*% sigma
            gradient <- - A %*% log.obs.i * (1/obs.i) - 
                    1/2 * (1/obs.i) * mtp +  phi.diff*tau.tilde.beta*(1/obs.i) 

            diagHessian <- -diag(A) * (1/obs.i^2) + 
                A %*% log.obs.i * (1/obs.i)^2 + 
                1/2 * (1/obs.i)^2 * mtp + (-phi.diff*tau.tilde*tau.tilde.beta^2+phi.diff*tau.tilde.beta-phi.diff^2*tau.tilde.beta^2)/(obs.i^2)

            weights <- do.call(what = "weightFun", args = c(ellipsis, 
                x = list(obs.i)))
            dWeights <- do.call(what = "dWeightFun", args = c(ellipsis, 
                x = list(obs.i)))

            sum(2 * (weights * dWeights) * gradient + weights^2 * 
                diagHessian + 1/2 * weights^2 * gradient^2)
        }
    }
    if(model=="truncT"){
        n <- nrow(loc)
        SigmaS = cov.func(loc, par2[idx.para])
        sigmaInv <- chol2inv(chol(SigmaS))
        a_fun <- function(j,upper=rep(Inf,n-1)){
            set.seed(747380)
            sigma_j = (SigmaS[-j,-j] - SigmaS[-j,j,drop=F] %*% SigmaS[j,-j,drop=F]/SigmaS[j,j])/(dof + 1)/SigmaS[j,j]
            val = mvtnorm::pmvt(lower=-SigmaS[-j,j]/SigmaS[j,j],upper=rep(Inf,n-1),sigma=sigma_j,df=dof+1)[[1]]
            assign(".Random.seed", oldSeed, envir=globalenv())
            return(log(val))
        }
        set.seed(747380)
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
    if (!is.null(ncores)) {
        scores <- parallel::mclapply(1:length(obs), computeScores, mc.cores = ncores,mc.set.seed = TRUE)
    }
    else {
        scores <- lapply(1:length(obs), computeScores)
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
    return(sum(unlist(scores))/length(obs))
}

fit.scoreMatching <- function(init, obs, loc,fixed=c(F,F,F,F,F), model="logskew", vario.func=NULL,cov.func=NULL,basis=NULL, idx.para=1:2, alpha.func=NULL, dof=2, weightFun = NULL, dWeightFun = NULL, method="Nelder-Mead", maxit=1000, ncores = NULL,lb,ub,trace=FALSE,step2=TRUE, partial=FALSE, ...){
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
    fixed2 = fixed
    fun <- function(par){
        par2 = init
        par2[!fixed2] = par
        if(any(par < lb[!fixed2]) | any(par > ub[!fixed2])){return(Inf)}
        val = scoreMatching(par2, obs, loc, model, vario.func, cov.func, alpha.func, basis, idx.para, dof, weightFun=weightFun, dWeightFun=dWeightFun,  ncores, partial=partial,...)
        return(val)
    }
    if(method=="L-BFGS-B"){
        opt.result = optim(init[!fixed],lower=lb[!fixed],upper=ub[!fixed],fun,method=method,control=list(maxit=maxit,trace=trace,factr=1e6))
    }else{
        opt.result = optim(init[!fixed],fun,method=method,control=list(maxit=maxit,trace=trace,reltol=1e-6))
    }
    if(model=="logskew" & any(!fixed[-idx.para]) & step2 & !is.null(ncores)){
        n.alpha = sum(!fixed[-idx.para])
        if(n.alpha==2){
            a = seq(0,2*pi,length.out=ncores)
            a = cbind(cos(a),sin(a))
        } else {
            a = matrix(rnorm(ncores*n.alpha),ncol=n.alpha)
            a = sweep(a,1,sqrt(rowSums(a^2)),FUN="/")
            a[,1] = pmax(a[,1],1)
        }
        init[!fixed] = opt.result$par
        fixed2[-idx.para] = fixed[-idx.para]
        fixed2[idx.para] = TRUE;
        init.list = split(a,row(a)) 
        if(method=="L-BFGS-B"){
            opt.result2 = mapply(optim,par=init.list,MoreArgs = list(fn=fun,lower=lb[!fixed2],upper=ub[!fixed2],method=method,control=list(maxit=maxit,trace=FALSE,factr=1e6),hessian=FALSE),SIMPLIFY=FALSE)
        }else{
            opt.result2 = mapply(optim,par=init.list,MoreArgs = list(fn=fun,method=method,control=list(maxit=maxit,trace=FALSE,reltol=1e-6)),SIMPLIFY=FALSE)
        }
        opt.values <- unlist(lapply(opt.result2,function(x){tryCatch(x$value,error=function(e){return(Inf)})}))
        opt.result = opt.result2[[which.min(opt.values)]]
        init[!fixed2] = opt.result$par
        fixed2 = fixed
        if(method=="L-BFGS-B"){
            opt.result = optim(init[!fixed2],lower=lb[!fixed2],upper=ub[!fixed2],fun,method=method,control=list(maxit=maxit,trace=trace,factr=1e6))
        }else{
            opt.result = optim(init[!fixed2],fun,method=method,control=list(maxit=maxit,trace=trace,reltol=1e-6))
        }
        #opt.result$others = opt.result2
    }
    if(model == "truncT" & !is.null(ncores) & step2){
        init[!fixed] <- opt.result$par
        fixed2 = fixed
        init.list = cbind(seq(0.5,log(init[1]),length.out=ncores),init[2])
        init.list = split(init.list,row(init.list))
        opt.result2 = mapply(optim,par=init.list,MoreArgs = list(fn=fun,method=method,control=list(maxit=maxit,trace=FALSE,reltol=1e-6)),SIMPLIFY=FALSE)
        opt.values <- unlist(lapply(opt.result2,function(x){tryCatch(x$par[1],error=function(e){return(Inf)})}))
        opt.result = opt.result2[[which.min(opt.values)]]
        init[!fixed2] = opt.result$par
    }
    t2 = proc.time() - t1
    init[!fixed2] = opt.result$par
    opt.result$par = init
    opt.result$time = t2
    return(opt.result)
}



