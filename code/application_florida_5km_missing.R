rm(list=ls())
args <- commandArgs(TRUE)
computer = "local"
stage = "all" # "spectral", "score", "diag" or "all"
ncores = 10
for (arg in args) eval(parse(text = arg))
switch(computer,
    "ws" = {DataPath<-"~/Desktop/InteriorExtremes/"},
    "hpc" = {DataPath<-"/srv/scratch/z3536974/";.libPaths("../src")},
    "local" = {DataPath<-"./"}
)
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(matrixStats)
library(Matrix)
library(numDeriv)
source("code/simulation.R")
source("code/exponent_functions.R")
set.seed(12342)
## Florida rainfall at the 3x3 block (~5.7km) hourly resolution with dry hours treated as missing,
## u_p set to the 95% empirical quantile of r_p, and spectral likelihood fitted to {r_0 > u_p} (Theorem 3)

file.spectral = paste0(DataPath,"data/application_florida_5km_missing_spectral.RData")
file.score = paste0(DataPath,"data/application_florida_5km_missing_score.RData")
file.diag = paste0(DataPath,"data/application_florida_5km_missing_diag.RData")

### unconditional margins: empirical ranks for the body, GPD tail above the 95% wet quantile; dry hours are missing ###
load(paste0(DataPath,"data/florida_agg/hourly_b3.RData")) # Y, coord.utm, cnt from code/florida_agg_prep.R
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
        # floor the tail probability so that bounded GPD tails cannot map observations to infinity
        Fe[ex] <- 1-pmax(mean(ex)*z^(-1/pr[1]),0.1/(n+1))
    }
    1/(1-Fe)
}
wet <- Y > 0
Z <- apply(Y,2,to.pareto); Z[!wet] <- 0 # dry entries are missing; zeros only enter the risk functionals
rm(Y); invisible(gc())
D <- ncol(Z); Nt <- nrow(Z)

### model set-up ###
loc <- apply(coord.utm/60000,2,function(x){x-median(x)})
basis.centers <- as.matrix(expand.grid(quantile(loc[,1],c(0.2,0.8)),quantile(loc[,2],c(0.2,0.8))))
basis <- sapply(1:nrow(basis.centers),function(i){
    y=dnorm(sqrt((loc[,1]-basis.centers[i,1])^2 + (loc[,2]-basis.centers[i,2])^2),mean=0,sd=ncol(loc)*2)
    y=y-mean(y)
    y/sqrt(sum(y^2))
})
lb = c(0.01,0.01,-pi/4,0.01,rep(-Inf,4))
ub = c(Inf,1.99,pi/4,Inf,rep(Inf,4))

vario_fast <- function(loc,par,ncores=NULL){
    lambda = par[1];alpha = par[2];theta = par[3];a = par[4]
    V = matrix(c(cos(theta),a*sin(theta),-sin(theta),a*cos(theta)),2,2)
    L = loc %*% t(V); g0 = (sqrt(rowSums(L^2))/lambda)^alpha
    outer(g0,g0,"+") - (as.matrix(dist(L))/lambda)^alpha
}
stopifnot(max(abs(vario_fast(loc[1:20,],c(0.3,1.2,0.3,0.8)) - vario.func2(loc[1:20,],c(0.3,1.2,0.3,0.8)))) < 1e-10)
par.list <- function(par){
    sigma = vario_fast(loc,par[1:4]); alpha = c(basis %*% par[5:8])
    list(sigma,c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha)))
}
lp.norm <- function(Z,p){
    m = rowMaxs(Z)
    r = if(is.infinite(p)) m else m*rowSums((Z/m)^p)^(1/p)
    r[!is.finite(r)] = 0 # all-dry hours
    r
}

### two-step selection: r_p exceedances of u_p, then r_0 > u_p (normalised r_0 <= r_p) ###
r0 <- rowSums(Z)/D
sets <- list()
for(p in c(1,10,Inf)){ # r_1 = r_0: all r_1 exceedances are retained
    rp = lp.norm(Z,p)/(if(is.infinite(p)) 1 else D^(1/p))
    u = quantile(rp,0.95,names=FALSE,type=1)
    idx = which(r0 > u); idx.p = which(rp > u)
    stopifnot(all(rp[idx] > u))
    sets[[length(sets)+1]] <- list(p=p,q=0.95,u=u,u.unnorm=u*(if(is.infinite(p)) 1 else D^(1/p)),idx=idx,idx.p=idx.p,acc=length(idx)/length(idx.p))
}
names(sets) <- sapply(sets,function(s){paste0("p=",s$p,",q=",s$q)})
print(sapply(sets,function(s) c(n.p=length(s$idx.p),n=length(s$idx),acc=round(s$acc,4))))
# storm clusters: exceedances less than 24 hours apart
cluster.id <- function(idx){cumsum(c(1,diff(idx) > 24))}
cpu <- function(t) sum(t[c(1,2,4,5)],na.rm=TRUE)/60 # total CPU minutes including forked workers

# sandwich standard errors with scores summed within storm clusters
godambe <- function(fun.vec,par,free,cl,r=4){
    f <- function(t){pp = par; pp[free] = t; fun.vec(pp)}
    S = numDeriv::jacobian(f,par[free],method="simple")
    H = numDeriv::hessian(function(t) sum(f(t)),par[free],method.args=list(r=r))
    Hi = solve(H)
    list(se=sqrt(pmax(diag(Hi %*% crossprod(rowsum(S,cl)) %*% Hi),0)),se.hess=sqrt(pmax(diag(solve(-H)),0)))
}
# events sharing the same wet set are evaluated together; chunks of events are evaluated in parallel
event.chunks <- function(Wm,nc){
    key = apply(Wm,1,function(w) paste(which(w),collapse=","))
    grp = unname(split(seq_len(nrow(Wm)),key))
    grp = grp[order(-sapply(grp,function(g) sum(Wm[g[1],])))]
    split(grp,rep(seq_len(nc),length.out=length(grp)))
}
par.apply <- function(chunks,FUN,n){
    # intensity_logskew_constraint needs .Random.seed in forked workers
    v = mclapply(chunks,function(gs){set.seed(747380); lapply(gs,function(g) list(g=g,v=FUN(g)))},mc.cores=length(chunks))
    out = numeric(n)
    for(ch in v){
        if(inherits(ch,"try-error")) stop(ch)
        for(e in ch) out[e$g] = e$v
    }
    out
}

### spectral log-likelihood on the wet sets U, evaluated with fit.model(partial=TRUE) ###
# the normalising constant of {sum_U x_j > c} is |U|/c and does not depend on the parameters;
# the sub-vector slant does not sum to zero, which fit.model accounts for (intensity_logskew_constraint does not)
vario.fm <- function(loc,par,ncore=NULL) vario_fast(loc,par)
wet.list <- function(X,Wm) lapply(seq_len(nrow(X)),function(i){U = which(Wm[i,]); list(U,X[i,U])})
fit.fm <- function(dat,init,fixed,nc=ncores) fit.model(data=dat,loc=loc,init=init,fixed=fixed,model="logskew",maxit=1000,FUN=vario.fm,basis=basis,alpha.func=alpha.func,
    ncores=nc,method="L-BFGS-B",lb=lb,ub=ub,opt=TRUE,idx.para=1:4,partial=TRUE,step2=FALSE)
# per-event contributions, for the sandwich standard errors
ll.events <- function(dat,nc=ncores){
    function(par) unlist(mclapply(dat,function(d) -fit.model(data=list(d),loc=loc,init=par,model="logskew",FUN=vario.fm,basis=basis,
        alpha.func=alpha.func,idx.para=1:4,partial=TRUE,opt=FALSE)$value,mc.cores=nc))
}
# check against the complete-data intensity, for which the slant sums to zero
local({
    i = which(rowSums(wet)==D)[1:3]; par = c(0.2,1.2,0.3,0.95,0.1,0.05,-0.05,0.02)
    v1 = ll.events(wet.list(Z[i,,drop=FALSE],wet[i,,drop=FALSE]),1)(par)
    v2 = intensity_logskew_constraint(Z[i,,drop=FALSE],list(vario_fast(loc,par[1:4]),c(basis %*% par[5:8])),log=TRUE)
    stopifnot(max(abs(v1-v2)) < 1e-8*max(abs(v2)))
})

### spectral likelihood fits ###
# risk functionals already fitted are kept; only the missing ones are fitted
if(stage %in% c("spectral","all")){
    fit.spectral <- function(s){
        X = Z[s$idx,,drop=FALSE]/s$u; Wm = wet[s$idx,,drop=FALSE]; cl = cluster.id(s$idx)
        dat = wet.list(X,Wm); ll = ll.events(dat)
        fit <- function(init,free) tryCatch(fit.fm(dat,init,!(1:8 %in% free)),error=function(e){message("start failed: ",conditionMessage(e)); list(value=Inf)})
        t0 = proc.time()
        starts = list(c(0.2,1,0,1,0,0,0,0),c(0.19,1.2,0.4,0.94,0,0,0,0),c(0.22,1.2,-0.4,1.06,0,0,0,0))
        br.all = lapply(starts,fit,free=1:4)
        br = br.all[[which.min(sapply(br.all,function(x) x$value))]]
        br$time = cpu(proc.time()-t0)
        t0 = proc.time()
        sbr.all = lapply(list(c(0.05,0.05,0,0),c(-0.05,-0.05,0,0)),function(b0){init = br$par; init[5:8] = b0; fit(init,1:8)})
        sbr = sbr.all[[which.min(sapply(sbr.all,function(x) x$value))]]
        sbr$time = cpu(proc.time()-t0)
        g.br = tryCatch(godambe(ll,br$par,1:4,cl),error=function(e) list(se=rep(NA,4),se.hess=rep(NA,4)))
        g.sbr = tryCatch(godambe(ll,sbr$par,1:8,cl),error=function(e) list(se=rep(NA,8),se.hess=rep(NA,8)))
        br$se = g.br$se; br$se.hess = g.br$se.hess; sbr$se = g.sbr$se; sbr$se.hess = g.sbr$se.hess
        br$aic = 2*nrow(X)*br$value + 2*4; sbr$aic = 2*nrow(X)*sbr$value + 2*8
        list(n=nrow(X),ncl=max(cl),miss=1-mean(Wm),br=br,sbr=sbr,br.values=sapply(br.all,function(x) x$value))
    }
    res.spectral <- if(file.exists(file.spectral)) local({load(file.spectral); res.spectral}) else list()
    todo = setdiff(names(sets),names(res.spectral))
    if(length(todo)){
        for(nm in todo){res.spectral[[nm]] = fit.spectral(sets[[nm]]); cat(nm,"done\n")}
        res.spectral = res.spectral[names(sets)]
        save(res.spectral,sets,loc,basis,basis.centers,coord.utm,file=file.spectral)
    }
}

### gradient score matching of the BR model on the r_p exceedances, restricted to the wet sets ###
# same weight functions as weightFun and dWeightFun in pareto_inference.R, vanishing on {r_p(x) = 1}
make.score <- function(X,Wm,p,nc=ncores){
    m = rowMaxs(X); r = if(is.infinite(p)) m else m*rowSums((X/m)^p)^(1/p)
    e = exp(-(r-1)); chunks = event.chunks(Wm,nc)
    function(par){
        SigmaS = vario_fast(loc,par[1:4])
        par.apply(chunks,function(g){
            U = which(Wm[g[1],]); S = SigmaS[U,U,drop=FALSE]
            Si = tryCatch(chol2inv(chol(S)),error=function(e) NULL)
            if(is.null(Si)) return(rep(NA,length(g)))
            q = rowSums(Si); A = Si - tcrossprod(q)/sum(q)
            mtp = c(2*q/sum(q) + 2 + A %*% diag(S))
            XU = X[g,U,drop=FALSE]; AL = log(XU) %*% A; iX = 1/XU
            G = -AL*iX - 0.5*sweep(iX,2,mtp,"*")
            DH = -sweep(iX^2,2,diag(A),"*") + AL*iX^2 + 0.5*sweep(iX^2,2,mtp,"*")
            eg = e[g]; W = XU*(1-eg)
            dW = (1-eg) + XU*eg*(if(is.infinite(p)) (XU==r[g]) else (XU/r[g])^(p-1))
            rowSums(2*W*dW*G + W^2*DH + 0.5*W^2*G^2)
        },nrow(X))
    }
}
if(stage %in% c("score","all")){
    # check against scoreMatching on complete observations
    source("code/pareto_inference.R")
    for(p in c(1,Inf)){
        s = sets[[paste0("p=",p,",q=0.95")]]; i = s$idx.p[rowSums(wet[s$idx.p,])==D][1:10]
        obs = Z[i,,drop=FALSE]/s$u.unnorm
        rFun <- if(is.infinite(p)) function(x) max(x) else function(x) max(x)*sum((x/max(x))^p)^(1/p)
        wF <- function(x){x*(1-exp(1-rFun(x)))}
        dwF <- function(x){r = rFun(x); (1-exp(-(r-1))) + x*exp(-(r-1))*(if(is.infinite(p)) (x==r) else (x/r)^(p-1))}
        v.repo = scoreMatching(c(0.21,1.3,0.4,0.93,0,0,0,0),split(obs,row(obs)),loc,"logskew",vario_fast,NULL,alpha.func,basis,1:4,weightFun=wF,dWeightFun=dwF,ncores=NULL)
        v.new = mean(make.score(obs,wet[i,,drop=FALSE],p,1)(c(0.21,1.3,0.4,0.93)))
        stopifnot(abs(v.repo - v.new) < 1e-8*abs(v.repo))
    }
    fit.score <- function(s){
        obs = Z[s$idx.p,,drop=FALSE]/s$u.unnorm; Wm = wet[s$idx.p,,drop=FALSE]; cl = cluster.id(s$idx.p)
        sc = make.score(obs,Wm,s$p)
        t0 = proc.time()
        # started at the complete-data score matching estimate (8.3 km, 1.6, 0.79, 0.95) to shorten the optimisation path
        fit = optim(c(8.3/60,1.6,0.7,0.95),function(par){v = mean(sc(par)); if(is.finite(v)) v else 1e10},method="L-BFGS-B",lower=lb[1:4],upper=ub[1:4],control=list(maxit=1000,factr=1e6))
        fit$time = cpu(proc.time()-t0)
        fit$par = c(fit$par,0,0,0,0)
        g = tryCatch(godambe(function(par) sc(par),fit$par,1:4,cl),error=function(e) list(se=rep(NA,4),se.hess=rep(NA,4)))
        fit$se = g$se; fit$se.hess = g$se.hess
        list(n=nrow(obs),ncl=max(cl),miss=1-mean(Wm),fit=fit)
    }
    # each risk functional is saved as soon as it is fitted
    res.score <- list()
    for(nm in names(sets)){
        f.nm = sub("\\.RData$",paste0("_",gsub("[^0-9a-zA-Z]+","_",nm),".RData"),file.score)
        if(file.exists(f.nm)){load(f.nm)} else {r.nm = fit.score(sets[[nm]]); save(r.nm,file=f.nm)}
        res.score[[nm]] = r.nm; cat(nm,"done\n")
    }
    save(res.score,file=file.score)
}

### diagnostics: retained proportions and bivariate extremal coefficients ###
if(stage %in% c("diag","all") & !file.exists(file.diag)){
    load(file.spectral); load(file.score)
    acceptance <- do.call(rbind,lapply(names(sets),function(nm){
        s = sets[[nm]]
        data.frame(set=nm,p=s$p,q=s$q,n.p=length(s$idx.p),n=length(s$idx),empirical=s$acc)
    }))
    print(acceptance)
    idx.ref <- sapply(list(c(0.25,0.25),c(0.25,0.75)),function(pr){
        target = c(quantile(loc[,1],pr[1]),quantile(loc[,2],pr[2]))
        which.min(colSums((t(loc)-target)^2))
    })
    emp.extcoef <- sapply(idx.ref,function(j){
        u = 100
        sapply(1:D,function(i) 2 - 2*sum(Z[,i]>u & Z[,j]>u)/(sum(Z[,i]>u)+sum(Z[,j]>u)))
    })
    fitted.extcoef <- function(par){
        pl = par.list(par)
        sapply(idx.ref,function(j){
            sapply(1:D,function(i){
                if(i==j) return(1)
                pair = c(i,j)
                V_bi_logskew(c(1,1),list(pl[[1]][pair,pair],pl[[2]][pair]),alpha.para=FALSE)
            })
        })
    }
    set.seed(747380)
    extcoef <- lapply(names(sets),function(nm){
        list(BR=fitted.extcoef(res.spectral[[nm]]$br$par),sBR=fitted.extcoef(res.spectral[[nm]]$sbr$par),BR.sm=fitted.extcoef(res.score[[nm]]$fit$par))
    })
    names(extcoef) <- names(sets)
    save(acceptance,idx.ref,emp.extcoef,extcoef,file=file.diag)
}
