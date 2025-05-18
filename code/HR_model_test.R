d = 3
loc = matrix(rnorm(d*2),ncol=2)*10

alpha = rnorm(d)
data = exp(rnorm(d))

cov.mat = matrix(runif(d^2),ncol=d);cov.mat = (cov.mat + t(cov.mat))/2
cov.mat = cov.mat + diag(d)*d

delta=alpha2delta(list(cov.mat,alpha))[[2]]
par = list(cov.mat,delta)

## HR model without skewness ##
intensity_logskew(data,list(cov.mat=cov.mat,delta=rep(0,d)),alpha.para=FALSE,log=FALSE)
for(i in 1:d){
    print(intensity_HR(data,par,i))
}

for(i in 1:d){
    print(intensity_skewedHR(data,list(cov.mat=cov.mat,delta=rep(0,d)),i))
}

V_logskew(data,list(cov.mat,alpha=rep(0,d)),alpha.para=TRUE)
for(i in 1:d){
    print(V_HR(data,list(cov.mat,alpha=rep(0,d)),i))
}

for(i in 1:d){
    print(V_skewedHR(data,list(cov.mat,alpha=rep(0,d)),i))
}

## skewed-HR model ##
intensity_logskew(data,par,alpha.para=FALSE,log=FALSE)
for(i in 1:d){
    print(intensity_skewedHR(data,par,i))
}

V_logskew(data,par,alpha.para=FALSE)
for(i in 1:d){
    print(V_skewedHR(data,par,i))
}


## 
d = 10
coord = as.matrix(expand.grid(1:d,1:d))
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
# loading library and setting path
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(partitions)
library(Rfast)
library(matrixStats)
library(splines)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
source("code/pareto_inference.R")


Lpnorm <- function(x,xi=1){
    val = sum((x)^xi)^{1/xi}
    return(val)
}
xi_vec = c(2,3,5,10)
d_vec = c(2,4,8,10)
for(i in 1:length(xi_vec)){
    for(j in 1:length(d_vec)){
        idx = coord[,1] <= d_vec[j] & coord[,2] <= d_vec[j]
        sigma = vario.func(coord[idx,],c(2,1))
        D = nrow(sigma)
        samples.skew.normal <- simu_Pareto_logskew(m=10^5,par=list(sigma,rep(0,D)),sum,ncores=NULL)
        xi = xi_vec[i]
        Lp = apply(samples.skew.normal,1,Lpnorm,xi=xi)
        a = format(mean(Lp>D^(xi-2))*100,scientific=TRUE,digits=2)
        b = format(mean(Lp>1)*100,scientific=TRUE,digits=2)
        val[i,j] = paste0(a,"/",b)
    }
}

kableExtra::kable(val,format="latex",digits=2)
save(val,file="code/simulation_study_speed.RData")
