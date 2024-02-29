## internal functions: do not use any of them directly!
simu_px_brownresnick <- function(no.simu=1, idx,  N, trend, chol.mat) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  res <- t(chol.mat)%*%matrix(rnorm(N*no.simu), ncol=no.simu)
  if (!is.matrix(trend)) {
    res <- exp(t(res - trend))
  } else {
    res <- exp(t(res - trend[,idx]))   
  }
  return(res/res[cbind(1:no.simu,idx)])
}

simu_px_extremalt <- function(no.simu=1, idx, N, dof, mu, chol.mat) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  if (!is.list(chol.mat)) {  
    res <- t(chol.mat)%*%matrix(rnorm(N*no.simu), ncol=no.simu)
    res <- sqrt((dof+1)/rgamma(no.simu, shape=0.5*(dof+1), scale=2))*(t(res) - res[idx,])  
    res <- pmax(matrix(mu, nrow=no.simu, ncol=N, byrow=TRUE)+res,0)^dof
  } else {
    res <- matrix(NA, nrow=no.simu, ncol=N) 
    for (i in 1:no.simu) {
      res[i,] <- as.vector(t(chol.mat[[idx[i]]])%*%matrix(rnorm(N), ncol=1))
      res[i,] <- sqrt((dof+1)/rgamma(1, shape=0.5*(dof+1), scale=2))*(res[i,] - res[i,idx[i]])
      res[i,] <- pmax(mu[[idx[i]]]+res[i,],0)^dof
    }  
  }
  return(res)
}       
       
simu_px_logistic <- function(no.simu=1, idx, N, theta) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  res       <- matrix(1/gamma(1-theta)*(-log(runif(no.simu*N)))^(-theta), nrow=no.simu, ncol=N) 
  res[cbind(1:no.simu,idx)] <- 1/gamma(1-theta)*rgamma(no.simu,shape=1-theta)^(-theta) 
  return(res/res[cbind(1:no.simu,idx)])
}

simu_px_neglogistic <- function(no.simu=1, idx, N, theta) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  res       <- matrix(rweibull(no.simu*N, shape=theta, scale=1/gamma(1+1/theta)), nrow=no.simu, ncol=N)
  res[cbind(1:no.simu,idx)] <- 1/gamma(1+1/theta)*rgamma(no.simu,shape=1+1/theta)^(1/theta) 
  return(res/res[cbind(1:no.simu,idx)])
}

simu_px_dirichlet <- function(no.simu, idx, N, weights, alpha, norm.alpha) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  if (length(idx)==1) {
    k <- sample(1:length(weights), no.simu, replace=TRUE, prob=N*weights*norm.alpha[idx,])
  } else {
    k <- sapply(1:no.simu, function(i) sample(1:length(weights), 1, prob=N*weights*norm.alpha[idx[i],]))
  }
  shape.mat <- alpha[,k,drop=FALSE]
  shape.mat[cbind(idx,1:no.simu)] <- shape.mat[cbind(idx,1:no.simu)]+1
  res <- t(matrix(rgamma(N*no.simu, shape=shape.mat), nrow=N, ncol=no.simu))
  return(res/res[cbind(1:no.simu,idx)])
}


## main functions

simu_extrfcts <- function(model, N, loc=1, scale=1, shape=1, no.simu=1, 
                          coord, cov.mat, corr, dof, theta, weights, alpha) {
                            
  stopifnot(model %in% c("brownresnick", "extremalt", "logistic", "neglogistic", "dirichlet"))
  
  if (model %in% c("brownresnick", "extremalt")) {
    stopifnot(!missing(coord))
    if (!is.matrix(coord)) coord <- matrix(coord, ncol=1)   
    if (!missing(N)) stopifnot(N==nrow(coord))   
    N <- nrow(coord)   
  }
  
  stopifnot((N==round(N)) & (N>=1))
  stopifnot((no.simu==round(no.simu)) & (no.simu>=1))
  
  if (length(loc)  ==1) loc   <- rep(loc  , times=N)
  if (length(scale)==1) scale <- rep(scale, times=N)
  if (length(shape)==1) shape <- rep(shape, times=N)
  stopifnot(all(scale>1e-12))
  
  if (model=="brownresnick") {
    # stopifnot(is.function(vario))
    # cov.mat <- sapply(1:N, function(i) sapply(1:N, function(j) 
    #                  vario(coord[i,]) + vario(coord[j,]) - vario(coord[c(i,j),])))
    cov.mat <- cov.mat + 1e-6 
    #add constant random effect to avoid numerical problems            
    chol.mat <- chol(cov.mat)
  } else if (model=="extremalt") {
    stopifnot(is.function(corr))
    stopifnot(dof>1e-12)
    diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),
                         as.vector(outer(coord[,2],coord[,2],'-')))  
    cov.mat.tmp <- matrix(apply(diff.vector, 1, function(x) corr(x)), ncol=N)       
  } else if (model=="logistic") {
    stopifnot(1e-12 < theta & theta < 1 - 1e-12)    
  } else if (model=="neglogistic") {
    stopifnot(theta > 1e-12) 
  } else if (model=="dirichlet") {
    m <- length(weights)
    stopifnot(all(weights>=0))
    stopifnot(abs(sum(weights)-1)<1e-12)
    stopifnot(length(alpha)==N*m)
    stopifnot(all(alpha>1e-12))
    if (N > 1 & m > 1) {
      stopifnot(is.matrix(alpha))
    } else {
      if (!is.matrix(alpha)) dim(alpha) <- c(N,m)
    }
    stopifnot(all(dim(alpha)==c(N,m)))
    norm.alpha <- apply(alpha, 2, function(alpha) alpha/sum(alpha))
    stopifnot(all(abs(norm.alpha %*% weights - rep(1/N, times=N))<1e-12))  
  }
   
  res <- matrix(0, nrow=no.simu, ncol=N)
  counter <- rep(0, times=no.simu)
   
  for (k in 1:N) {
    poisson <- rexp(no.simu)
    if (model == "brownresnick") {
      #trend <- sapply(1:N, function(j) vario(coord[c(j,k),]))
        trend <- diag(cov.mat) - cov.mat[k,]
    } else if (model == "extremalt") {
      cov.vec  <- apply(coord, 1, function(x) corr(x-coord[k,]))
      cov.mat  <- (cov.mat.tmp - outer(cov.vec, cov.vec, '*'))/(dof+1) + 1e-6
      chol.mat <- chol(cov.mat)
      mu <- apply(coord, 1, function(x) corr(coord[k,]-x))    
    }
    while (any(1/poisson > res[,k])) {
      ind <- (1/poisson > res[,k])
      n.ind <- sum(ind)
      idx <- (1:no.simu)[ind]
      counter[ind] <- counter[ind] + 1
      proc <- switch(model,
                "brownresnick" = simu_px_brownresnick(no.simu=n.ind, idx=k, N=N, trend=trend, chol.mat=chol.mat),
                "extremalt"    = simu_px_extremalt(no.simu=n.ind, idx=k, N=N, dof=dof, mu=mu, chol.mat=chol.mat),
                "logistic"     = simu_px_logistic(no.simu=n.ind, idx=k, N=N, theta=theta),
                "neglogistic"  = simu_px_neglogistic(no.simu=n.ind, idx=k, N=N, theta=theta),
                "dirichlet"    = simu_px_dirichlet(no.simu=n.ind, idx=k, N=N, weights=weights, alpha=alpha, norm.alpha=norm.alpha)
              )
      stopifnot(dim(proc)==c(n.ind, N))
      if (k==1) {
        ind.upd <- rep(TRUE, times=n.ind)
      } else {
        ind.upd <- sapply(1:n.ind, function(i) 
                                   all(1/poisson[idx[i]]*proc[i,1:(k-1)] <= res[idx[i],1:(k-1)]))
      }
      if (any(ind.upd)) {
        idx.upd <- idx[ind.upd]
        res[idx.upd,] <- pmax(res[idx.upd,], 1/poisson[idx.upd]*proc[ind.upd,])
      }
      poisson[ind] <- poisson[ind] + rexp(n.ind)
    } 
  }
  res <- sapply(1:N, function(i) {
           if (abs(shape[i]<1e-12)) {   
             return(log(res[,i])*scale[i] + loc[i])
           } else {
             return(1/shape[i]*(res[,i]^shape[i]-1)*scale[i] + loc[i])
           }
         })   
   
  return(list(res=res, counter=counter))  
}


simu_specfcts <- function(model, N, loc=1, scale=1, shape=1, no.simu=1, 
                          coord, vario, corr, dof, theta, weights, alpha) {
 stopifnot(model %in% c("brownresnick", "extremalt", "logistic", "neglogistic", "dirichlet"))
  
  if (model %in% c("brownresnick", "extremalt")) {
    stopifnot(!missing(coord))
    if (!is.matrix(coord)) coord <- matrix(coord, ncol=1)   
    if (!missing(N)) stopifnot(N==nrow(coord))   
    N <- nrow(coord)   
  }
  
  stopifnot((N==round(N)) & (N>=1))
  stopifnot((no.simu==round(no.simu)) & (no.simu>=1))
  
  if (length(loc)  ==1) loc   <- rep(loc  , times=N)
  if (length(scale)==1) scale <- rep(scale, times=N)
  if (length(shape)==1) shape <- rep(shape, times=N)
  stopifnot(all(scale>1e-12))
  
  if (model=="brownresnick") {
    stopifnot(is.function(vario))
    cov.mat <- sapply(1:N, function(i) sapply(1:N, function(j) 
                     vario(coord[i,]) + vario(coord[j,]) - vario(coord[c(i,j),])))
    cov.mat <- cov.mat + 1e-6 
    #add constant random effect to avoid numerical problems            
    chol.mat <- chol(cov.mat)
    trend <- sapply(1:N, function(k) sapply(1:N, function(j) vario(coord[c(j,k),])))    
  } else if (model=="extremalt") {
    stopifnot(is.function(corr))
    stopifnot(dof>1e-12)
    diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),
                         as.vector(outer(coord[,2],coord[,2],'-')))  
    cov.mat.tmp <- matrix(apply(diff.vector, 1, function(x) corr(x)), ncol=N)   
    chol.mat <- vector("list", length=N)
    mu <- vector("list", length=N)
    for (k in 1:N) {
      cov.vec  <- apply(coord, 1, function(x) corr(x-coord[k,]))
      cov.mat  <- (cov.mat.tmp - outer(cov.vec, cov.vec, '*'))/(dof+1) + 1e-6
      chol.mat[[k]] <- chol(cov.mat)
      mu[[k]]       <- apply(coord, 1, function(x) corr(coord[k,]-x))   
    }    
  } else if (model=="logistic") {
    stopifnot(1e-12 < theta & theta < 1 - 1e-12)    
  } else if (model=="neglogistic") {
    stopifnot(theta > 1e-12) 
  } else if (model=="dirichlet") {
    m <- length(weights)
    stopifnot(all(weights>=0))
    stopifnot(abs(sum(weights)-1)<1e-12)
    stopifnot(length(alpha)==N*m)
    stopifnot(all(alpha>1e-12))
    if (N > 1 & m > 1) {
      stopifnot(is.matrix(alpha))
    } else {
      if (!is.matrix(alpha)) dim(alpha) <- c(N,m)
    }
    stopifnot(all(dim(alpha)==c(N,m)))
    norm.alpha <- apply(alpha, 2, function(alpha) alpha/sum(alpha))
    stopifnot(all(abs(norm.alpha %*% weights - rep(1/N, times=N))<1e-12))  
  }
   
  res <- matrix(0, nrow=no.simu, ncol=N)
  counter <- rep(0, times=no.simu)
 
  poisson <- rexp(no.simu)
  ind <- rep(TRUE, times=no.simu)
  while (any(ind)) {
    n.ind <- sum(ind)
    counter[ind] <- counter[ind] + 1
    shift <- sample(1:N, n.ind, replace=TRUE)
    proc <- switch(model,
                "brownresnick" = simu_px_brownresnick(no.simu=n.ind, idx=shift, N=N, trend=trend, chol.mat=chol.mat),
                "extremalt"    = simu_px_extremalt(no.simu=n.ind, idx=shift, N=N, dof=dof, mu=mu, chol.mat=chol.mat),
                "logistic"     = simu_px_logistic(no.simu=n.ind, idx=shift, N=N, theta=theta),
                "neglogistic"  = simu_px_neglogistic(no.simu=n.ind, idx=shift, N=N, theta=theta),
                "dirichlet"    = simu_px_dirichlet(no.simu=n.ind, idx=shift, N=N, weights=weights, alpha=alpha, norm.alpha=norm.alpha)
              )
    stopifnot(dim(proc)==c(n.ind, N))
    proc <- N*proc/rowSums(proc)     
    res[ind,] <- pmax(res[ind,], proc/poisson[ind])
    poisson[ind] <- poisson[ind] + rexp(n.ind)
    ind <- (N/poisson > apply(res, 1, min))
  }
  res <- sapply(1:N, function(i) {
           if (abs(shape[i]<1e-12)) {   
             return(log(res[,i])*scale[i] + loc[i])
           } else {
             return(1/shape[i]*(res[,i]^shape[i]-1)*scale[i] + loc[i])
           }
         })   
   
  return(list(res=res, counter=counter)) 
}