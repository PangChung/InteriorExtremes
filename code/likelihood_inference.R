###########################################################
###########################################################
## BROWN-RESNICK MAX-STABLE MODEL with General Variogram###
###########################################################
###########################################################
###########################
## Variogram functions ####
###########################
## semivariogram function but returns a covariance matrix
# vario.func <- function(loc,par){ 
#     oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
#     set.seed(747380)    
#     alpha = par[1];lambda = par[2];a = par[3]; theta = par[4]
#     #Sigma <- matrix(c(par[3],-par[4],-par[4],1),2,2)
#     A = matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),2,2)
#     Sigma <- A%*%diag(c(1,a),2)%*%t(A)
#     loc = matrix(loc,ncol=2)
#     n = nrow(loc)
#     if(n==1){
#         val=2*(sqrt(t(loc[1,])%*%Sigma%*%loc[1,])/lambda)^alpha
#         return(val)
#     }
#     fun <- function(idx){
#         loc.temp <- loc[idx,]
#         if(idx[1]==idx[2]){
#             h <- loc.temp[1,]
#             val = 2*(sqrt(t(h)%*%Sigma%*%h)/lambda)^alpha
#         }else{
#             h <- loc.temp[1,]-loc.temp[2,]
#             val <- (sqrt(t(loc.temp[1,])%*%Sigma%*%(loc.temp[1,]))/lambda)^alpha+
#             (sqrt(t(loc.temp[2,])%*%Sigma%*%(loc.temp[2,]))/lambda)^alpha-
#             (sqrt(t(h)%*%Sigma%*%h)/lambda)^alpha 
#       }
#       return(val)
#     }
#     idx <- cbind(rep(1:n, times = n:1),unlist(lapply(1:n, function(x){x:n})))
#     val <- apply(idx,1,fun)
#     val.mat <- matrix(1,n,n)
#     val.mat[idx]<-val;
#     val.mat[idx[,c(2,1)]]<-val
#     assign(".Random.seed", oldSeed, envir=globalenv())
#     return(val.mat + .Machine$double.eps * diag(n))
# }

### Exponent function V for the Brown-Resnick model
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# sigma: covariance matrix of dimension DxD
V <- function(data,sigma){
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)  
    if(!is.matrix(data)) data <- matrix(data,nrow=1)
    D <- ncol(data)
    if(D==1){
    return(1/data)
  } else{
    fun.i <- function(i){
      if(nrow(data)==1){
            eval.i <- t(log(data[,-i]/data[,i])) + diag(sigma)[-i]/2 + sigma[i,i]/2 - sigma[i,-i]
      }else{
            eval.i <- t(t(log(data[,-i]/data[,i]))) + diag(sigma)[-i]/2 + sigma[i,i]/2 - sigma[i,-i]
      }
      Id <- diag(1,D-1)
      if(i==1){ ### see Wadsworth and Tawn (2014) for the definition of the matrix T.i...
        T.i <- cbind(-1,Id)
      } else if(i==D){
        T.i <- cbind(Id,-1)
      } else{
        T.i <- cbind(Id[,1:(i-1)],-1,Id[,i:(D-1)]) 
      }
      sigma.i <- T.i%*%sigma%*%t(T.i)
      return(apply(eval.i,1,function(x){return(mvtnorm::pmvnorm(upper=x,sigma=sigma.i))})/data[,i])
    }
    if(nrow(data)==1){
        assign(".Random.seed", oldSeed, envir=globalenv())
        return( sum(sapply(1:D,fun.i)))
    }
    assign(".Random.seed", oldSeed, envir=globalenv())
    return( rowSums(sapply(1:D,fun.i)) )
  }
}

### Negative partial derivatives of the exponent function V for the Brown-Resnick model
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# sigma: covariance matrix of dimension DxD
# I: vector of indices with respect to which the partial derivatives are computed; if I==c(1:D), the function returns the full mixed derivative.
nVI <- function(data,sigma,I,logval=FALSE){
    oldSeed <- get(".Random.seed", mode="numeric", envir=globalenv())
    set.seed(747380)
    if(!is.matrix(data)) data <- matrix(data,nrow=1)
    D <- ncol(data)
    nI <- length(I)
    if(nI==0){ ## If the index set is empty
        return(-V(data,sigma))
    }else if(nI==D){ ## full derivative
        if(D==1){
        return(1/data^2)
    } else{
      sigma.DD <- diag(sigma)
      sigma.chol <- chol(sigma)
      sigma.inv <- chol2inv(sigma.chol)
      log.det.sigma <- 2*sum(log(diag(sigma.chol)))
      q <- rowSums(sigma.inv)
      q.sum <- sum(q)
      A <- (sigma.inv - q%*%t(q)/q.sum)
      sigma.q <- t(sigma.DD)%*%q
      
      log.data <- log(data)
      
      log.Part1 <- 0
      log.Part2 <- ((D-1)/2)*log(2*pi) + (1/2)*log.det.sigma + (1/2)*log(q.sum) + rowSums(log.data)
      log.Part3 <- c(-(1/2)*( 1/4*t(sigma.DD)%*%sigma.inv%*%sigma.DD -1/4*(sigma.q)^2/q.sum + sigma.q/q.sum - 1/q.sum))
      log.Part4 <- c(-(1/2)*(apply(log.data,1,function(x){return(t(x)%*%A%*%x)}) + log.data%*%(q%*%(2-sigma.q)/q.sum + sigma.inv%*%sigma.DD)))
      
      res <- log.Part1-log.Part2+log.Part3+log.Part4
      assign(".Random.seed", oldSeed, envir=globalenv())
      if(logval) return(drop(res)) else return(drop(exp(res)))
      return( drop(res) )
    }
  } else{ ## Partial derivative
    if(D==1){ ## If there is only one variable
      return(1/data^2)
    } else{
        sigma.DD <- diag(sigma)
        sigma.II <- sigma.DD[I]
        sigma.inv <- chol2inv(chol(sigma))
        q <- rowSums(sigma.inv)
        A <- (sigma.inv - q%*%t(q)/sum(q))
        
        
        sigma.I <- as.matrix(sigma[I,I])
        sigma.I.chol = chol(sigma.I)
        sigma.I.inv <- chol2inv(sigma.I.chol)
        logdet.sigma.I = 2*sum(log(diag(chol(sigma.I))))
        q.I <- rowSums(sigma.I.inv)
        A.I <- (sigma.I.inv - q.I%*%t(q.I)/sum(q.I))
        
        K10 <- matrix(0,nrow=D,ncol=nI)
        K10[I,] <- diag(nI)
        K01 <- matrix(0,nrow=D,ncol=D-nI)
        K01[-I,] <- diag(D-nI)
        
        log.data <- log(data)
        log.data.I <- matrix(log.data[,I],ncol=nI)
        
        gamma.inv <- t(K01)%*%A%*%K01
        tryCatch({
            gamma <- chol2inv(chol(gamma.inv))
            mu <- -gamma%*%(t(K01)%*%A%*%K10%*%t(log.data.I) + c(t(K01)%*%( (q - 1/2*q%*%t(q)%*%sigma.DD )/sum(q)+1/2*sigma.inv%*%sigma.DD )) )
            eval.x <- log(data[,-I])-t(mu)
            sigma.q.II <- t(sigma.II)%*%q.I
            q.I.sum <- sum(q.I)
            
            log.Part1 <- apply(eval.x,1,function(x){return(max(mvtnorm::pmvnorm(upper=x,sigma=gamma),0))})
            log.Part2 <- ((nI-1)/2)*log(2*pi) + (1/2)*logdet.sigma.I + (1/2)*log(q.I.sum) + rowSums(log.data.I)
            log.Part3 <- c(-(1/2)*( 1/4*t(sigma.II)%*%sigma.I.inv%*%sigma.II -1/4*(sigma.q.II)^2/q.I.sum + sigma.q.II/q.I.sum - 1/q.I.sum))
            log.Part4 <- c(-(1/2)*(apply(log.data.I,1,function(x){return(t(x)%*%A.I%*%x)}) + log.data.I%*%(q.I%*%(2-sigma.q.II)/q.I.sum + sigma.I.inv%*%sigma.II)))
            #res <- drop(log.Part1*exp(-log.Part2+log.Part3+log.Part4))
            res <- drop(log(log.Part1) - log.Part2+log.Part3+log.Part4)
            assign(".Random.seed", oldSeed, envir=globalenv())},error=function(e){browser()})
        if(logval) return(drop(res)) else return(drop(exp(res)))
    }
  }
}

### BIVARIATE Exponent function V for the Brown-Resnick model
# data: matrix of dimension nx2, containing n 2-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# sigma: covariance matrix of dimension 2x2
V.biv <- function(data,sigma){
  vario <- (sigma[1,1]+sigma[2,2]-2*sigma[1,2])
  a <- sqrt(vario)
  P12 <- a/2-log(data[,1]/data[,2])/a
  P21 <- a/2-log(data[,2]/data[,1])/a
  return(pnorm(P12)/data[,1]+pnorm(P21)/data[,2])
}

### Negative partial derivatives of the exponent function V for the Brown-Resnick model
# data: matrix of dimension nx2, containing n 2-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# sigma: covariance matrix of dimension 2x2
# I: vector of indices with respect to which the partial derivatives are computed; if I==c(1:2), the function returns the full mixed derivative.
nVI.biv <- function(data,sigma,I){
  nI <- length(I)
  if(nI==0){
    return(-V.biv(data,sigma))
  } else if(nI==2){
    vario <- (sigma[1,1]+sigma[2,2]-2*sigma[1,2])
    a <- sqrt(vario)
    P12 <- a/2-log(data[,1]/data[,2])/a
    P21 <- a/2-log(data[,2]/data[,1])/a
    return( (dnorm(P12)*(1-P12/a)/data[,1]+dnorm(P21)*(1-P21/a)/data[,2])/(a*data[,1]*data[,2]) )
  } else if(nI==1){
    vario <- (sigma[1,1]+sigma[2,2]-2*sigma[1,2])
    a <- sqrt(vario)
    P12 <- a/2-log(data[,1]/data[,2])/a
    P21 <- a/2-log(data[,2]/data[,1])/a
    if(I==1){
      return( (pnorm(P12)+dnorm(P12)/a)/data[,1]^2 - dnorm(P21)/(a*data[,1]*data[,2]) )
    } else if(I==2){
      return( (pnorm(P21)+dnorm(P21)/a)/data[,2]^2 - dnorm(P12)/(a*data[,1]*data[,2]) )
    }
  }
}

### Negative log likelihood function for Brown-Resnick data with unit Frechet margins
# par: parameter vector
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# loc: matrix of dimension DxD, containing all pairwise distances between locations or Locations
# FUN: the variogram function that returns the covraiance matrix
nloglik <- function(par,data,model="BR"){
    #fix random seed (and save the current random seed to restore it at the end)
    if(!is.matrix(data)){data <- matrix(data,nrow=1)}
    D <- ncol(data)
    if(D>2){
        all_combn <- lapply(1:D,FUN=Rfast::comb_n,n=D,simplify=FALSE) 
        all_nVI <- list() ## will contain all the terms nVI (total number is equal to 2^D-1), used later to assemble the log-likelihood...
        if(model == "BR"){
            sigma = par[[1]]
            all_nVI <- lapply(all_combn,FUN = function(idx){sapply(idx,nVI,data=data,sigma=sigma)})
            Vdata = V(data,sigma)
        }
        if(model == "logskew"){
            all_nVI <- lapply(all_combn,FUN = function(idx){sapply(idx,partialV_logskew,x=data,par=par,alpha.para=FALSE,log=FALSE)})
            Vdata = V_logskew(data,par,alpha.para=FALSE)
        }
        get.nVI <- function(I){
            nI <- length(I)
            return(all_nVI[[nI]][,which(sapply(all_combn[[nI]],function(x){return(all(I%in%x))}))])
        }
        parts <- listParts(D) ## using package `partitions'
        contribution.partition <- function(partition){
            return( matrixStats::rowProds(as.matrix(as.data.frame(lapply(partition,FUN=get.nVI) ))))
        }
        # tmp = log(rowSums(as.matrix(as.data.frame(lapply(parts,contribution.partition)))))
        # if(any(is.nan(tmp))){print(data);browser()}
        res <- log(rowSums(as.matrix(as.data.frame(lapply(parts,contribution.partition))))) - Vdata
    }else{
        if(model == "BR"){
            sigma = par[[1]]
            res <- log(nVI.biv(data,sigma,1:2)+nVI.biv(data,sigma,1)*nVI(data,sigma,2)) - V(data,sigma)
        }
        if(model == "logskew"){
            res <- log(partialV_logskew(data,1:2,par,alpha.para=FALSE)+partialV_logskew(data,1,par,alpha.para=FALSE)*partialV_logskew(data,2,par,alpha.para=FALSE)) - V_bi_logskew(data,par)
        }
    }
    if(any(!is.finite(res))){return(rep(Inf,nrow(data)))}
    return(-res)
}

###########################
## Composite likelihood ###
###########################
### Negative log composite-likelihood function for max-stable models with unit Frechet margins ###
# par: parameter vector (par[1]=range, par[2]=smoothness)
# data: matrix of dimension nxD, containing n D-dimensional random vectors (each row = 1 vector) on the unit Frechet scale
# index: q-by-Q matrix of q-dimensional margins to be used in the composite likelihood. Here Q refers to the number of composite likelihood contributions (with 1<=Q<=choose(D,q))
nlogcomplik <- function(par,data,index,ncores,model){
    if(model == "logskew"){
        par <- alpha2delta(par)  
    }
    nlogcomplik.contribution <- function(ind){
      par.index <- par
      if(model == "BR"){par.index[[1]] = par[[1]][ind,ind]}
      if(model == "logskew"){par.index[[1]] = par[[1]][ind,ind];par.index[[2]] = par.index[[2]][ind]} 
      val <- nloglik(par=par.index,data[,ind],model)
    }
    if(!is.null(ncores)) res <- mclapply(as.list(as.data.frame(index)),nlogcomplik.contribution,mc.cores = ncores,mc.set.seed = F,mc.preschedule = TRUE)
    else res = lapply(as.list(as.data.frame(index)),nlogcomplik.contribution)
    res <- mean(unlist(res))
    return(res)
}

### Function that returns the MCLE (maximum composite likelihood estimator) for the Brown-Resnick model with unit Frechet margins
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# init: initial parameter vector (init[1]=initial range, init[2]=initial smoothness)
# fixed: vector of booleans indicating whether or not the parameters are fixed to initial values
# loc: coordinates
# sigmaFUN: function returns covariance matrix
# index: q-by-Q matrix of q-dimensional margins to be used in the composite likelihood. Here Q refers to the number of composite likelihood contributions (with 1<=Q<=choose(D,q)).
MCLE <- function(data,init,fixed,loc,FUN,index,ncores,maxit=200,model="BR",hessian=FALSE,lb=-Inf,ub=Inf,alpha.func=NULL,trace =FALSE,basis=NULL,idx.para=1:2,...){
    t <- proc.time()
    
    object.func <- function(par2,opt=TRUE){
        par1 <- init
        par1[!fixed] <- par2
        if( any(par1 < lb) | any( par1 > ub)  ){return(Inf)}
        sigma = FUN(loc,par1[idx.para])
        if(model=="BR"){par.list <- list(sigma=sigma)}
        if(model=="logskew"){b.mat <- basis 
                            par.list <- list(sigma=sigma,alpha=alpha.func(par=par1[-idx.para],b.mat=b.mat))
                            par.list <- alpha2delta(par.list)}
        val = nlogcomplik(par.list,data=data,index,ncores,model=model)
        if(opt){ 
            val = mean(val,na.rm=TRUE)
            # print(c(par2,val))
        }
        return(val)
    }
    if(sum(!fixed)==1){
        opt <- optim(par=init[!fixed],fn=object.func,lower=lb[!fixed],upper=ub[!fixed],method="Brent",control=list(maxit=maxit,trace=trace),hessian=hessian)
    }else{
        opt <- optim(par=init[!fixed],fn=object.func,method="Nelder-Mead",control=list(maxit=maxit,trace=trace),hessian=hessian)
    }
    if(hessian){
        h = 1e-4
        par.mat.grad = matrix(opt$par,nrow=length(opt$par),ncol=length(opt$par),byrow=TRUE) + diag(h,length(opt$par))
        val.object = object.func(opt$par,opt=FALSE)
        val.object.grad = apply(par.mat.grad,1,function(x){(object.func(x,opt=FALSE) - val.object)/h})
        opt$K = var(val.object.grad)
        opt$hessian.inv = solve(opt$hessian)
        opt$sigma = opt$hessian.inv %*% opt$K %*% opt$hessian.inv
    }
    init[!fixed] = opt$par
    time <- proc.time()-t
    opt$par = init
    opt$time <- time[3]
  return(opt)
}

##############################
## Vecchia's approximation ###
##############################
### Negative log likelihood function for Max-stable data with unit Frechet margins, based on Vecchia's approximation
# par: parameter vector (par[1]=range, par[2]=smoothness)
# data: matrix of dimension nxD, containing n D-dimensional random vectors (each row = 1 vector) on the unit Frechet scale
# loc: coordinates
# FUN: function returns covariance matrix
# vecchia.seq: vector of length D (with integers from {1,...,D}), indicating the sequence of variables to be considered for the Vecchia approximation
# neighbours: an q-by-D matrix with the corresponding the neighbors of each observation in the Vecchia sequence (where q is the number of neighbours, i.e., the size of the conditioning set)
nlogVecchialik <- function(par,data,vecchia.seq,neighbours,ncores,model="BR"){
    if(model == "logskew"){
        par <- alpha2delta(par)  
    }
    logVecchialik.contribution <- function(i){
        par.index <- par
        if(i==1 & !any(!is.na(neighbours[,i]))){
            par.index[[1]] = par[[1]][vecchia.seq[1],vecchia.seq[1]]
            if(model == "logskew"){par.index[[1]] = par[[1]][vecchia.seq[1]]}  
            contribution <- nloglik(par.index,data[,vecchia.seq[1],drop=FALSE],model) #density of 1st variable in the sequence (unit Fréchet)
        }else{
            ind.i <- vecchia.seq[i] #index of ith-variable in the Vecchia sequence
            ind.neighbours <- na.omit(neighbours[,i])
            ind <- c(ind.i,ind.neighbours)
            par.index[[1]] = par[[1]][ind,ind]
            if(model == "logskew"){par.index[[2]] = par.index[[2]][ind]}  
            num <- nloglik(par.index,data[,ind,drop=FALSE],model) #joint density of ith-variable and its conditioning set
            par.index[[1]] = par[[1]][ind.neighbours,ind.neighbours]
            if(model == "logskew"){par.index[[2]] = par[[2]][ind.neighbours]}  
            denom <- nloglik(par.index,data[,ind.neighbours,drop=FALSE],model) #joint density of conditioning set only
            contribution <- num-denom
        }
    return(contribution)
    }
    if(!is.null(ncores)){
        res <- rowSums(matrix(unlist(mclapply(1:length(vecchia.seq),FUN=logVecchialik.contribution,mc.cores=ncores,mc.set.seed = F,mc.preschedule = TRUE)),ncol=length(vecchia.seq),byrow=FALSE),na.rm=TRUE)
    }else{
        res <- rowSums(matrix(unlist(lapply(1:length(vecchia.seq),FUN=logVecchialik.contribution)),ncol=length(vecchia.seq),byrow=FALSE),na.rm=TRUE)
    }
    return(res)
}

### Function that returns the MVLE (maximum Vecchia likelihood estimator) for the Brown-Resnick model with unit Frechet margins
# data: matrix of dimension nxD, containing n D-dimensional random Brown-Resnick vectors (each row = 1 vector) on the unit Frechet scale
# init: initial parameter vector (init[1]=initial range, init[2]=initial smoothness)
# fixed: vector of booleans indicating whether or not the parameters are fixed to initial values
# loc: coordinates
# sigma.FUN: function returns covariance matrix
# vecchia.seq: vector of length D (with integers from {1,...,D}), indicating the sequence of variables to be considered for the Vecchia approximation
# neighbours: an q-by-D matrix with the corresponding the neighbors of each observation in the Vecchia sequence (where q is the number of neighbours, i.e., the size of the conditioning set)
MVLE <- function(data,init,fixed,loc,FUN,vecchia.seq,neighbours,ncores,model="BR",maxit=1000,hessian=FALSE,alpha.func=NULL,lb=-Inf,ub=Inf,trace=FALSE,basis=NULL,idx.para=1:2,...){
    t <- proc.time()
    object.func <- function(par2,opt=TRUE){
        par1 <- init
        par1[!fixed] <- par2
        if( any(par1 < lb) | any( par1 > ub)  ){return(Inf)}
        sigma = FUN(loc,par1[idx.para])
        if(model=="BR"){par.list=list(sigma=sigma)}
        if(model=="logskew"){b.mat <- basis / sqrt(diag(sigma))
                            par.list <- list(sigma=sigma,alpha=alpha.func(par=par1[-idx.para],b.mat=b.mat))
                            par.list <- alpha2delta(par.list)}
        val = nlogVecchialik(par.list,data,vecchia.seq,neighbours,ncores,model)
        if(opt){
            val = mean(val,na.rm=TRUE)
            #print(c(par2, val))
        } 
        return(val)
    }
    if(sum(!fixed)==1){
        opt <- optim(par=init[!fixed],fn=object.func,method="Brent",lower=lb[!fixed],upper=ub[!fixed],control=list(maxit=maxit,trace=trace),hessian=hessian)
    }else{
        opt <- optim(par=init[!fixed],fn=object.func,method="Nelder-Mead",control=list(maxit=maxit,trace=trace),hessian=hessian)
        #opt <- optim(par=init[!fixed],fn=object.func,method="L-BFGS-B",lower=lb,upper=ub,control=list(maxit=maxit,trace=TRUE),hessian=hessian)
    }
    if(hessian){
        h = 1e-4
        par.mat.grad = matrix(opt$par,nrow=length(opt$par),ncol=length(opt$par),byrow=TRUE) + diag(h,length(opt$par))
        val.object = object.func(opt$par,opt=FALSE)
        val.object.grad = apply(par.mat.grad,1,function(x){(object.func(x,opt=FALSE) - val.object)/h})
        opt$K = var(val.object.grad)
        opt$hessian.inv = solve(opt$hessian)
        opt$sigma = opt$hessian.inv %*% opt$K %*% opt$hessian.inv
    }
    init[!fixed] = opt$par
    time <- proc.time()-t
    opt$par = init
    opt$time <- time[3]
  return( opt )
}

neighbours <- function(ind,vecchia.seq,q,loc){
    if(ind==1){
      ind.neighbours <- rep(NA,q)
    }
    if(ind>=2){
      ind.ind <- vecchia.seq[ind] #index of ith-variable in the Vecchia sequence
      ind.past <- vecchia.seq[1:(ind-1)] #index of the "past" observations in the Vecchia sequence  
      d.past <- loc[ind.past,vecchia.seq[ind]] #distance of the ith-variable to the "past" observations
      ind.neighbours <- ind.past[order(d.past)[1:q]] #choose "neighbours" as the closest observations in the "past"  
    }
    return(ind.neighbours)
}
