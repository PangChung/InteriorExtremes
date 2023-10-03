library(tmvnsim)
library(parallel)
library(evd)
library(doParallel)
source("code/simulation.R")
source("code/exponent_functions.R")
### testing the simulator ###
d <- 10
coord = as.matrix(expand.grid(1:d,1:d)/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),
                         as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
corr <- function(x,r=0.5,v=1) exp(- (sum(x^2)/r)^v)                          
cov.mat <- matrix(apply(diff.vector, 1, corr), ncol=nrow(coord)) + diag(1e-6,nrow(coord))       
chol(cov.mat)
nu = 10
par1 <- list(nu=nu,sigma=cov.mat)

# Simulate a truncated extremal-t max-stable process
system.time(Z.trunc <- simu_truncT(m=10000,par=par1,ncores=10))
hist(pgev(Z.trunc[,1],1,1,1),50,prob=TRUE)
image(1:10,1:10,z=matrix(log(Z.trunc[2,]),nrow=10),col=rev(heat.colors(10)) )

# Simulate a log-skew normal based max-stable process
alpha = rep(0.5,nrow(coord))
par2 <- list(alpha=alpha,sigma=cov.mat)
system.time(Z.logskew <- simu_logskew(m=10000,par=par2,ncores=10))
hist(pgev(Z.logskew[,1],1,1,1),50,prob=TRUE)
image(1:10,1:10,z=matrix(log(Z.logskew[2,]),nrow=10),col=rev(heat.colors(10)) )

# calculate empirical extremal coefficients
empirical_extcoef <- function(p,data){
    return(min(2,max(1,1/mean(1/pmax(data[,p[1]],data[,p[2]])))))
}

truc_extcoef <- function(idx,par,model="logskew1"){
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
        x = rep(Inf,nrow(par[[2]]))
        x[idx] = 1
        val = V_truncT(x,par,parallel=FALSE)
    }

    if(model == "truncT2"){
        x = rep(1,length(idx))
        val = V_truncT(x,list(nu = par[[1]],sigma=par[[2]][idx,idx]),parallel=FALSE)
    }
    return(val)

}

all.pairs <- combn(1:ncol(Z.trunc),2)
ec.trunc <- apply(all.pairs,2,empirical_extcoef,data=Z.trunc)
plot(x=diff.mat[t(all.pairs)],y=ec.trunc,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",pch=20)
tc.truncT1 <- apply(all.pairs,2,truc_extcoef,par=par1,model="truncT1")
tc.truncT2 <- apply(all.pairs,2,truc_extcoef,par=par1,model="truncT2")
points(x=diff.mat[t(all.pairs)],y=tc.truncT1,type="p",cex=0.5,col="red",pch=20)
points(x=diff.mat[t(all.pairs)],y=tc.truncT2,type="p",cex=0.5,col="blue",pch=20)

ec.logskew <- apply(all.pairs,2,empirical_extcoef,data=Z.logskew)
tc.logskew1 <- apply(all.pairs,2,truc_extcoef,par=par2,model="logskew1")
tc.logskew2 <- apply(all.pairs,2,truc_extcoef,par=par2,model="logskew2")
plot(x=diff.mat[t(all.pairs)],y=ec.logskew,type="p",cex=0.5,ylim=c(1,2),xlab="Distance",ylab="Extremal Coefficient",pch=20)
points(x=diff.mat[t(all.pairs)],y=tc.logskew1,type="p",cex=0.5,col="red",pch=20)
points(x=diff.mat[t(all.pairs)],y=tc.logskew2,type="p",cex=0.5,col="blue",pch=20)
abline(h=c(1,2),col="grey",lty=2)


