rm(list=ls())
args <- commandArgs(TRUE)
computer = "local"
stage = "all" # "prep", "spectral", "score", "diag" or "all"
ncores = 10
for (arg in args) eval(parse(text = arg))
switch(computer,
    "ws" = {DataPath<-"~/Desktop/InteriorExtremes/"},
    "hpc" = {DataPath<-"/srv/scratch/z3536974/";.libPaths("../src")},
    "local" = {DataPath<-"~/Documents/Github/InteriorExtremes/"}
)
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(partitions)
library(Rfast)
library(matrixStats)
library(Matrix)
library(numDeriv)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
source("code/pareto_inference.R")
set.seed(12342)

file.prep = paste0(DataPath,"data/application_florida_5km.RData")
file.spectral = paste0(DataPath,"data/application_florida_5km_spectral.RData")
file.score = paste0(DataPath,"data/application_florida_5km_score.RData")
file.diag = paste0(DataPath,"data/application_florida_5km_diag.RData")

### upscale the radar data to 3x3 pixel blocks (~5.7km) and hourly totals ###
if(!file.exists(file.prep)){
    library(sf)
    load("data/application_florida.RData")
    ids <- as.numeric(sub("P","",names(tmp)[-1]))
    grid.sf <- read_sf("data/Florida/DOPGrid/DOPGrid.shp")
    xy <- suppressWarnings(st_coordinates(st_centroid(grid.sf)))[match(ids,grid.sf$PIXEL),1:2]
    X <- as.matrix(tmp[,-1]); rm(tmp); invisible(gc())
    # the radar grid is rotated by ~0.43 degrees, so rows and columns are indexed in the rotated frame
    phi <- atan2(14,1858)
    xr <- xy %*% t(matrix(c(cos(phi),sin(phi),-sin(phi),cos(phi)),2))
    col <- round((xr[,1]-min(xr[,1]))/1888); row <- round((xr[,2]-min(xr[,2]))/1991)
    b <- 3
    key <- paste(col%/%b,row%/%b); tb <- table(key); keep.key <- names(tb)[tb >= b*b/2]
    keep <- key %in% keep.key; f <- factor(key[keep],levels=keep.key); cnt <- as.numeric(table(f))
    M <- sparseMatrix(i=which(keep),j=as.integer(f),x=1/cnt[as.integer(f)],dims=c(ncol(X),length(keep.key)))
    Y <- as.matrix(X %*% M); rm(X); invisible(gc())
    Y <- matrix(colSums(matrix(Y,nrow=4)),ncol=ncol(Y)) # hourly totals
    coord.utm <- as.matrix(t(M) %*% xy)

    # unconditional margins: empirical ranks for the body, GPD tail (probability weighted moments) above the 95% wet quantile
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
    Z <- apply(Y,2,to.pareto)
    wet <- Y > 0
    save(Z,wet,coord.utm,cnt,file=file.prep)
}else{
    load(file.prep)
}
D <- ncol(Z); Nt <- nrow(Z)
if(stage=="prep") stop("data preparation finished")

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

# vectorised version of vario.func2
vario_fast <- function(loc,par,ncores=NULL){
    lambda = par[1];alpha = par[2];theta = par[3];a = par[4]
    V = matrix(c(cos(theta),a*sin(theta),-sin(theta),a*cos(theta)),2,2)
    L = loc %*% t(V); g0 = (sqrt(rowSums(L^2))/lambda)^alpha
    outer(g0,g0,"+") - (as.matrix(dist(L))/lambda)^alpha
}
stopifnot(max(abs(vario_fast(loc[1:20,],c(0.3,1.2,0.3,0.8)) - vario.func2(loc[1:20,],c(0.3,1.2,0.3,0.8)))) < 1e-10)

lp.norm <- function(Z,p){
    m = rowMaxs(Z)
    if(is.infinite(p)) m else m*rowSums((Z/m)^p)^(1/p)
}

### Holder selection: {||z||_1 > D^(1-1/p) u_p} is a subset of {||z||_p > u_p} ###
L1 <- rowSums(Z)
sets <- list()
for(np in c(5000,8000)){
    for(p in c(3,10,Inf)){
        rp = lp.norm(Z,p)
        u = sort(rp,decreasing=TRUE)[np+1]
        idx = which(L1 > D^(1-1/p)*u*(1+1e-12))
        stopifnot(all(rp[idx] > u))
        sets[[length(sets)+1]] <- list(p=p,q=round(1-np/Nt,3),u=u,idx=idx,idx.p=which(rp > u),acc=length(idx)/np)
    }
}
names(sets) <- sapply(sets,function(s){paste0("p=",s$p,",q=",s$q)})
# storm clusters: exceedances less than 24 hours apart
cluster.id <- function(idx){cumsum(c(1,diff(idx) > 24))}

# sandwich standard errors with scores summed within storm clusters
godambe <- function(fun.vec,par,free,cl,r=4){
    f <- function(t){pp = par; pp[free] = t; fun.vec(pp)}
    S = numDeriv::jacobian(f,par[free],method="simple")
    H = numDeriv::hessian(function(t) sum(f(t)),par[free],method.args=list(r=r))
    Hi = solve(H)
    sqrt(pmax(diag(Hi %*% crossprod(rowsum(S,cl)) %*% Hi),0))
}

### spectral likelihood fits using fit.model ###
if(stage %in% c("spectral","all") & !file.exists(file.spectral)){
    fit.spectral <- function(s){
        set.seed(747380) # intensity_logskew_constraint needs .Random.seed in forked workers
        X = Z[s$idx,,drop=FALSE]; cl = cluster.id(s$idx)
        fit <- function(init,fixed) fit.model(data=X,loc=loc,init=init,fixed=fixed,model="logskew",maxit=1000,FUN=vario_fast,basis=basis,alpha.func=alpha.func,ncores=NULL,method="L-BFGS-B",lb=lb,ub=ub,opt=TRUE,idx.para=1:4,partial=FALSE,step2=FALSE)
        ll.obs <- function(par) intensity_logskew_constraint(X,list(vario_fast(loc,par[1:4]),c(basis %*% par[5:8])),log=TRUE)
        starts = list(c(0.2,1,0,1,0,0,0,0),c(0.19,1.2,0.4,0.94,0,0,0,0),c(0.22,1.2,-0.4,1.06,0,0,0,0))
        br.all = lapply(starts,fit,fixed=c(F,F,F,F,T,T,T,T))
        br = br.all[[which.min(sapply(br.all,function(x) x$value))]]
        sbr.all = lapply(list(c(0.05,0.05,0,0),c(-0.05,-0.05,0,0)),function(b0){init = br$par; init[5:8] = b0; fit(init,rep(F,8))})
        sbr = sbr.all[[which.min(sapply(sbr.all,function(x) x$value))]]
        br$se = tryCatch(godambe(ll.obs,br$par,1:4,cl),error=function(e) rep(NA,4))
        sbr$se = tryCatch(godambe(ll.obs,sbr$par,1:8,cl),error=function(e) rep(NA,8))
        list(n=nrow(X),ncl=max(cl),br=br,sbr=sbr,br.values=sapply(br.all,function(x) x$value))
    }
    res.spectral <- mclapply(sets,fit.spectral,mc.cores=min(ncores,length(sets)),mc.set.seed=TRUE)
    save(res.spectral,sets,loc,basis,basis.centers,file=file.spectral)
}

### gradient score matching fits of the BR model using the r_p exceedances directly ###
# vectorised version of scoreMatching(model="logskew") for the BR model (alpha = 0, complete vectors):
# the matrices that scoreMatching recomputes for every observation do not depend on the observation
score_br <- function(par,X,p){
    SigmaS = vario_fast(loc,par[1:4]); sigmaInv = chol2inv(chol(SigmaS))
    q = rowSums(sigmaInv); A = sigmaInv - tcrossprod(q)/sum(q)
    mtp = c(2*q/sum(q) + 2 + A %*% diag(SigmaS))
    LX = log(X); AL = LX %*% A; iX = 1/X
    G = -AL*iX - 0.5*sweep(iX,2,mtp,"*")
    DH = -sweep(iX^2,2,diag(A),"*") + AL*iX^2 + 0.5*sweep(iX^2,2,mtp,"*")
    # weight functions vanishing on {r_p(x) = 1}, as weightFun and dWeightFun in pareto_inference.R
    m = rowMaxs(X); r = if(is.infinite(p)) m else m*rowSums((X/m)^p)^(1/p)
    e = exp(-(r-1)); W = X*(1-e)
    dW = (1-e) + X*e*(if(is.infinite(p)) (X==r) else (X/r)^(p-1))
    rowSums(2*W*dW*G + W^2*DH + 0.5*W^2*G^2)
}
if(stage %in% c("score","all") & !file.exists(file.score)){
    # check against scoreMatching on a few observations
    for(p in c(3,Inf)){
        s = sets[[paste0("p=",p,",q=0.912")]]; obs = Z[s$idx.p[1:20],,drop=FALSE]/s$u
        rFun <- if(is.infinite(p)) function(x) max(x) else function(x) max(x)*sum((x/max(x))^p)^(1/p)
        wF <- function(x){x*(1-exp(1-rFun(x)))}
        dwF <- function(x){r = rFun(x); (1-exp(-(r-1))) + x*exp(-(r-1))*(if(is.infinite(p)) (x==r) else (x/r)^(p-1))}
        v.repo = scoreMatching(c(0.21,1.3,0.4,0.93,0,0,0,0),split(obs,row(obs)),loc,"logskew",vario_fast,NULL,alpha.func,basis,1:4,weightFun=wF,dWeightFun=dwF,ncores=NULL)
        stopifnot(abs(v.repo - mean(score_br(c(0.21,1.3,0.4,0.93),obs,p))) < 1e-8*abs(v.repo))
    }
    fit.score <- function(s){
        obs = Z[s$idx.p,,drop=FALSE]/s$u; cl = cluster.id(s$idx.p)
        t0 = proc.time()
        # same optimiser settings as fit.scoreMatching
        fit = optim(c(0.2,1,0.3,0.95),function(par) mean(score_br(par,obs,s$p)),method="L-BFGS-B",lower=lb[1:4],upper=ub[1:4],control=list(maxit=1000,factr=1e6))
        fit$time = proc.time() - t0
        fit$par = c(fit$par,0,0,0,0)
        fit$se = tryCatch(godambe(function(par) score_br(par,obs,s$p),fit$par,1:4,cl),error=function(e) rep(NA,4))
        list(n=nrow(obs),ncl=max(cl),fit=fit)
    }
    res.score <- mclapply(sets,fit.score,mc.cores=min(ncores,length(sets)),mc.set.seed=TRUE)
    save(res.score,file=file.score)
}

### diagnostics: acceptance probabilities and bivariate extremal coefficients ###
if(stage %in% c("diag","all") & !file.exists(file.diag)){
    load(file.spectral); load(file.score)
    par.list <- function(par){
        sigma = vario_fast(loc,par[1:4]); alpha = c(basis %*% par[5:8])
        list(sigma,c(sigma %*% alpha)/sqrt(c(1+alpha %*% sigma %*% alpha)))
    }
    # model-implied acceptance probability D^(1/p-1)/E||w||_p, where w are the L1 angles of the r_0-Pareto process
    acc.model <- function(par,p,m=4000){
        set.seed(747380)
        W = simu_Pareto_logskew(m,par.list(par),riskr=sum,ncores=ncores)
        W = W/rowSums(W)
        D^(1/p-1)/mean(lp.norm(W,p))
    }
    acceptance <- do.call(rbind,lapply(names(sets),function(nm){
        s = sets[[nm]]; r = res.spectral[[nm]]
        data.frame(set=nm,p=s$p,q=s$q,n.p=length(s$idx.p),n=length(s$idx),empirical=s$acc,BR=acc.model(r$br$par,s$p),sBR=acc.model(r$sbr$par,s$p))
    }))
    print(acceptance)

    # reference blocks in the south-west and north-west of the domain
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
