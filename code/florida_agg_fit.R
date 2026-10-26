# BR r-Pareto fits by the spectral likelihood with dry sites treated as missing, for candidate aggregations
args <- commandArgs(TRUE); ncores <- 10
for (arg in args) eval(parse(text = arg))
suppressMessages({library(parallel); library(matrixStats); library(numDeriv); library(mvtnorm)})
source("code/exponent_functions.R")
set.seed(12342)
pwm.gpd <- function(y){
    y <- sort(y); n <- length(y); p <- (seq_len(n)-0.35)/n
    a0 <- mean(y); a1 <- mean((1-p)*y)
    c(2-a0/(a0-2*a1), 2*a0*a1/(a0-2*a1))
}
to.pareto <- function(y){
    n <- length(y); Fe <- rank(y,ties.method="average")/(n+1)
    thr <- quantile(y[y>0],0.95,names=FALSE); ex <- y>thr
    if(sum(ex)>=30){
        pr <- pwm.gpd(y[ex]-thr)
        z <- pmax(1+pr[1]*(y[ex]-thr)/pr[2],1e-300)
        Fe[ex] <- 1-pmax(mean(ex)*z^(-1/pr[1]),0.1/(n+1))
    }
    1/(1-Fe)
}
lp.norm <- function(Z,p){ m = rowMaxs(Z); if(is.infinite(p)) m else m*rowSums((Z/m)^p)^(1/p) }
vario_fast <- function(loc,par){
    lambda = par[1];alpha = par[2];theta = par[3];a = par[4]
    V = matrix(c(cos(theta),a*sin(theta),-sin(theta),a*cos(theta)),2,2)
    L = loc %*% t(V); g0 = (sqrt(rowSums(L^2))/lambda)^alpha
    outer(g0,g0,"+") - (as.matrix(dist(L))/lambda)^alpha
}
lb = c(0.01,0.01,-pi/4,0.01); ub = c(Inf,1.99,pi/4,Inf)

configs <- list(c(3,1),c(4,3),c(6,3),c(4,6),c(6,6),c(8,6))
prep <- function(b,h){
    load(sprintf("data/florida_agg/hourly_b%d.RData",b))
    Yh <- if(h==1) Y else rowsum(Y,rep(seq_len(nrow(Y)/h),each=h),reorder=FALSE)
    wet <- Yh > 0; Z <- apply(Yh,2,to.pareto); Z[!wet] <- 0
    loc <- apply(coord.utm/60000,2,function(x){x-median(x)})
    list(b=b,h=h,Z=Z,wet=wet,loc=loc,D=ncol(Z))
}
jobs <- list()
for(cf in configs) for(p in c(3,10,Inf)) jobs[[length(jobs)+1]] <- list(b=cf[1],h=cf[2],p=p,q=0.86)

fit.job <- function(job,inner=1){
    set.seed(747380)
    d <- prep(job$b,job$h); D <- d$D; Z <- d$Z; loc <- d$loc
    r0 <- rowSums(Z)/D
    rp <- lp.norm(Z,job$p)/(if(is.infinite(job$p)) 1 else D^(1/job$p))
    rp[!is.finite(rp)] <- 0 # all-dry time steps
    u <- quantile(rp,job$q,names=FALSE,type=1)
    idx <- which(r0 > u); idx.p <- which(rp > u)
    X <- Z[idx,,drop=FALSE]/u; W <- d$wet[idx,,drop=FALSE]
    # group events sharing the same wet set so that each submatrix is factorised once
    key <- apply(W,1,function(w) paste(which(w),collapse=","))
    grp <- split(seq_along(idx),key)
    ll.grp <- function(gs,sigma){
        unlist(lapply(gs,function(g){
            U = which(W[g[1],])
            if(length(U)==1) -2*log(X[g,U]) else
                intensity_logskew_constraint(X[g,U,drop=FALSE],list(sigma[U,U],rep(0,length(U))),alpha.para=FALSE,log=TRUE)
        }))
    }
    # events split into chunks evaluated in parallel when inner > 1
    chunks <- split(grp,rep(seq_len(inner),length.out=length(grp)))
    ord.ev <- unlist(lapply(chunks,unlist))
    ll.vec <- function(par){
        sigma = vario_fast(loc,par); out = numeric(length(idx))
        v = if(inner==1) ll.grp(grp,sigma) else unlist(mclapply(chunks,ll.grp,sigma=sigma,mc.cores=inner))
        out[if(inner==1) unlist(grp) else ord.ev] = v
        out
    }
    nll <- function(par){ v = -sum(ll.vec(par)); if(!is.finite(v)) 1e10 else v }
    t0 <- proc.time()
    starts <- list(c(0.2,1,0,1),c(0.4,1.2,0.4,0.9),c(0.4,1.2,-0.4,1.1))
    fits <- lapply(starts,function(s) optim(s,nll,method="L-BFGS-B",lower=lb,upper=ub,control=list(maxit=1000,factr=1e7)))
    fit <- fits[[which.min(sapply(fits,function(f) f$value))]]
    time <- (proc.time()-t0)[3]
    # sandwich standard errors with scores summed within storms (gaps of at least 24 hours)
    cl <- cumsum(c(1,diff(idx) > 24/job$h))
    se <- tryCatch({
        S = numDeriv::jacobian(ll.vec,fit$par,method="simple")
        H = numDeriv::hessian(function(t) sum(ll.vec(t)),fit$par,method.args=list(r=4))
        Hi = solve(H); list(sand=sqrt(pmax(diag(Hi %*% crossprod(rowsum(S,cl)) %*% Hi),0)),hess=sqrt(pmax(diag(solve(-H)),0)))
    },error=function(e) list(sand=rep(NA,4),hess=rep(NA,4)))
    # model-implied retention 1/kappa(A_p) = D^(1/p-1)/E r'_p(W) with L1 angles W of the BR r_0-Pareto process
    sigma = vario_fast(loc,fit$par); m = 4000
    G = matrix(rnorm(m*D),m) %*% chol(sigma); j = sample(D,m,replace=TRUE)
    V = exp(G + t(sigma[,j]) - matrix(diag(sigma)/2,m,D,byrow=TRUE)); V = V/rowSums(V)
    ret.model = D^(1/job$p-1)/mean(lp.norm(V,job$p))
    # pairwise extremal coefficients: empirical (marginal level 0.98) against fitted
    E = (Z > 50)*1; ne = colSums(E); both = crossprod(E)
    emp = 2 - 2*both/outer(ne,ne,"+")
    fitc = 2*pnorm(sqrt(pmax(outer(diag(sigma),diag(sigma),"+")-2*sigma,0))/2)
    up = upper.tri(emp); dist.km = as.matrix(dist(loc))[up]*60
    rmse = sqrt(mean((emp[up]-fitc[up])^2)); rmse50 = sqrt(mean(((emp[up]-fitc[up])^2)[dist.km<50]))
    list(job=job,D=D,u=u,n.p=length(idx.p),n=length(idx),ncl=max(cl),npat=length(grp),miss=1-mean(W),
         par=fit$par,value=fit$value,conv=fit$convergence,se=se,time=time,
         ret.emp=length(idx)/length(idx.p),ret.model=ret.model,rmse=rmse,rmse50=rmse50,
         theta.bins=tapply(emp[up]-fitc[up],cut(dist.km,c(0,10,20,40,80,160,400)),mean))
}
heavy <- sapply(jobs,function(j) j$b==3 && j$h==1 && j$p==3)
res <- mclapply(jobs[!heavy],function(j) tryCatch(fit.job(j),error=function(e) list(job=j,err=conditionMessage(e))),mc.cores=ncores,mc.preschedule=FALSE)
save(res,file="data/florida_agg/fits_missing_BR.RData")
cat("light jobs done\n")
res[[length(res)+1]] <- fit.job(jobs[[which(heavy)]],inner=ncores)
save(res,file="data/florida_agg/fits_missing_BR.RData")
