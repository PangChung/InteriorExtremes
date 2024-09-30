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
