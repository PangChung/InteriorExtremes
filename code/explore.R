args <- commandArgs(TRUE)
computer = "local"
id = 1
d <- 15 ## 10 * 10 grid on [0,1]^2
m <- 1000 ## number of samples
model = "logskew"; # "logskew" or "truncT"
#model = "truncT"; # "logskew" or "truncT"
for (arg in args) eval(parse(text = arg))
switch(computer,
    "ws" = {DataPath<-"~/Desktop/InteriorExtremes/"},
    "hpc" = {DataPath<-"/srv/scratch/z3536974/";.libPaths("../src")},
    "local" = {DataPath<-"~/Documents/Github/InteriorExtremes/"}
)
# settings 
coord = as.matrix(expand.grid(1:d,1:d))
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
para.range = c(4,8) #c(0.5,1,2) ## range for the correlation function ##      
para.nu = c(1,1.5) #c(0.5,1,1.5) ## smoothness parameter for the correlation function ##
para.alpha = rbind(c(0,0),c(-1,-2),c(-1,1)) ## slant parameter for skewed norm model ##
para.deg = 2 ## degree of the freedom for the truncated t model ##
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
thres = c(0.95,0.9)
# loading library and setting path
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(gridExtra)
library(partitions)
library(ggplot2)
library(Rfast)
library(matrixStats)
library(splines)
library(numDeriv)
library(cubature)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
ncores=detectCores()
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
set.seed(init.seed)

## compute the basis ###
centers <- rbind(c(0.25,0.25),c(0.5,0.5),c(0.75,0.75))*d
idx.centers <- apply(centers,1,function(x){which.min(apply(coord,1,function(y){sum((x-y)^2)}))})
basis <- sapply(idx.centers,function(x){y=dnorm(diff.mat[x,],mean=0,sd=d*2);y=y-mean(y);y/max(abs(y))})
# basis <- basis * 100
summary(basis)

basis <- cbind(bs(coord[,1],degree = 1),coord)
basis <- apply(basis,2,function(x){x-mean(x)})
basis[,1] = basis[,1]*20
summary(basis)


## contour of the bi-variate extremal coefficient ##
coord = as.matrix(expand.grid(1:d,1:d))
sigma = cov.func(coord,c(16,1))
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))

centers <- rbind(c(0.25,0.25),c(0.5,0.5),c(0.75,0.75))*d
idx.centers <- apply(centers,1,function(x){which.min(apply(coord,1,function(y){sum((x-y)^2)}))})
basis <- sapply(idx.centers,function(x){y=dnorm(diff.mat[x,],mean=0,sd=32);y=y-mean(y);y/sqrt(sum(y^2))})

basis <- cbind(bs(coord[,1],degree = 1),coord)
basis <- apply(basis,2,function(x){x-mean(x)})

idx = floor(matrix(seq(1,nrow(coord),length.out=6),ncol=2,3))
basis <- sapply(1:(ncol(para.alpha)+1),function(x){y <- rep(0,nrow(coord));y[idx[x,]] <- c(-2,2);y})

all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
idx.center = c(18,18)
idx.center = which.min(abs(coord[,1] - idx.center[1]) + abs(coord[,2] - idx.center[2]))
ind.idx.center = all.pairs[1,] == idx.center |  all.pairs[2,] == idx.center
bi.extcoef <- list()
p1.list <- list()
for(i in 1:nrow(para.alpha)){
    alpha = alpha.func(par=para.alpha[i,])
    bi.extcoef[[i]] = mcmapply(true_extcoef,idx=all.pairs.list[ind.idx.center],MoreArgs=list(par=alpha2delta(list(sigma,alpha))),mc.cores=ncores,mc.set.seed = FALSE)
    
    data <- data.frame( x = coord[-idx.center,1],
                    y = coord[-idx.center,2],
                    z = bi.extcoef[[i]])
    p1.list[[i]] <- ggplot(data, aes(x = x, y = y, z = z)) +
    geom_contour(aes(colour = after_stat(level))) +
    scale_colour_gradient(low = "blue", high = "red") +
    theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
    labs(title = paste("Bivariate Extremal Coef.:",paste(para.alpha[i,],collapse = " ")), x = "X", y = "Y", colour = "Z")
}

grid.arrange(grobs=p1.list,nrow=1)


sigma = cov.func(coord,c(4,1))
delta = alpha2delta(list(sigma,alpha.func(para.alpha[2,])))[[2]]
alpha = delta2alpha(list(sigma,delta))[[2]]
range(delta)
alpha = delta2alpha(list(sigma,delta*d))[[2]]
alpha = delta2alpha(list(sigma,rep(1/sqrt(15*10),d*d)))[[2]]
range(alpha)
# # plot contours #

# 
# pdf(file="figures/simulation_samples_extcoef_contours.pdf",width=10,height = 8,onefile = TRUE)
# do.call(grid.arrange, c(p1.list[1:9], ncol = 3,nrow=3))
# do.call(grid.arrange, c(p1.list[1:9+9], ncol = 3,nrow=3))
# do.call(grid.arrange, c(p1.list[1:9+18], ncol = 3,nrow=3))
# dev.off()


idx = floor(matrix(seq(1,nrow(coord),length.out=6),ncol=2,3))
basis <- sapply(1:(ncol(para.alpha)+1),function(x){y <- rep(0,nrow(coord));y[idx[x,]] <- c(-2,2);y})
summary(basis)
#basis[,1] <- rep(0,nrow(basis))
# 1: fixing the first one, 2: unfixing the first one
#pdf("figures/delta_basis_1.pdf",width=5*3,height = 5,onefile = TRUE)
pdf("figures/delta_basis_4.pdf",width=5*3,height = 5,onefile = TRUE)
for(idx in 1:nrow(para.alpha)){
    beta1 = alpha2delta(list(cov.func(coord,c(4,1)),alpha.func(para.alpha[idx,])))[[2]]
    df = data.frame(x = coord[,1], y = coord[,2], z = beta1)
    p1 <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_tile() + ggtitle(paste(c("alpha:",para.alpha[idx,]),collapse = " ")) + 
    scale_fill_gradient2(low = "blue",mid="white" ,high = "red",limits=c(-1,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5))
    #p1
    beta2 = alpha2delta(list(cov.func(coord,c(4,1)),alpha.func(2*para.alpha[idx,])))[[2]]
    df = data.frame(x = coord[,1], y = coord[,2], z = beta2)
    p2 <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_tile() + ggtitle(paste(c("alpha:",2*para.alpha[idx,]),collapse = " ")) + 
    scale_fill_gradient2(low = "blue",mid="white",high = "red",limits=c(-1,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5))
    #p2
    df = data.frame(x = coord[,1], y = coord[,2], z = beta1-beta2)
    p3 <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_tile() + ggtitle(paste("Max difference: ",max(abs(df$z)))) +
    scale_fill_gradient2(low = "blue",mid="white",high = "red",limits=c(-1,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5))
    #p3
    grid.arrange(p1,p2,p3,ncol=3)
}
dev.off()

## explore the basis to use ##
alpha.1 = seq(-5,5,0.1)
alpha.addon = 4
alpha.grid = as.matrix(expand.grid(alpha.1,alpha.1))
alpha.grid.list <- split(alpha.grid,row(alpha.grid))

delta.grid <- lapply(alpha.grid.list,function(x){alpha2delta(list(cov.func(coord[1:3,],c(2,1)),c(x,-sum(x))+c(0.5,0,-0.5)*alpha.addon))[[2]]})
delta.grid2 <- lapply(alpha.grid.list,function(x){alpha2delta(list(cov.func(coord[1:3,],c(2,1)),c(x,-sum(x))+c(0.5,0,-0.5)*alpha.addon*2))[[2]]})
delta.diff <- lapply(1:length(delta.grid),function(x){sum(abs(delta.grid[[x]]-delta.grid2[[x]]))})

biv.ext <- lapply(alpha.grid.list,function(x){V_logskew(c(1,1),par=list(matrix(c(1,0.5,0.5,1),2,2),x),alpha.para=TRUE)})

#data = data.frame(x=alpha.grid[,1],y=alpha.grid[,2],z=unlist(delta.diff))
z=matrix(matrix(unlist(delta.grid),ncol=3,byrow=TRUE)[,3],length(alpha.1),length(alpha.1),byrow = TRUE)
z=matrix(unlist(biv.ext),length(alpha.1),length(alpha.1),byrow = TRUE)
p <- plotly::plot_ly(z=~z, type = "surface")
p 

########
diff.delta <- function(par,sigma){
    omega = diag(sqrt(diag(sigma)))
    omega.inv = diag(diag(omega)^(-1))
    sigma.bar = omega.inv %*% sigma %*% omega.inv
    fun <- function(i){
        alpha.1 = alpha.func(par=par[i,])
        alpha.2 = alpha.func(par=par[i,]*2)
        delta.1 = c(sigma.bar %*% alpha.1)/sqrt(c(1+t(alpha.1) %*% sigma.bar %*% alpha.1))
        delta.2 = c(sigma.bar %*% alpha.2)/sqrt(c(1+alpha.2 %*% sigma.bar %*% alpha.2))
        max(abs(delta.1-delta.2))
    }
    return(unlist(mclapply(1:nrow(par),fun,mc.cores=ncores)))
}

centers <- rbind(c(0.25,0.25),c(0.5,0.5),c(0.75,0.75))*d
idx.centers <- apply(centers,1,function(x){which.min(apply(coord,1,function(y){sum((x-y)^2)}))})
basis <- sapply(idx.centers,function(x){y=dnorm(diff.mat[x,],mean=0,sd=d*2);y=y-mean(y);y/sqrt(sum(y^2))})
#basis <- sapply(idx.centers,function(x){y=exp(-diff.mat[x,]/d/2);y=y-mean(y);y/max(abs(y))})

basis <- cbind(bs(coord[,1],degree = 1),coord)
basis <- apply(basis,2,function(x){x-mean(x)})

idx = floor(matrix(seq(1,nrow(coord),length.out=6),ncol=2,3))
basis <- sapply(1:(ncol(para.alpha)+1),function(x){y <- rep(0,nrow(coord));y[idx[x,]] <- c(-2,2);y})

alpha.1 = seq(-5,5,0.1)
sigma = cov.func(coord,c(4,1.5))
alpha.grid = as.matrix(expand.grid(alpha.1,alpha.1))
values <- diff.delta(alpha.grid,sigma)
summary(values)

para.norm = as.matrix(expand.grid(para.range,para.nu))
p.list = list()
for(i in 1:nrow(para.norm)){
    sigma = cov.func(coord,para.norm[i,])
    alpha.grid = as.matrix(expand.grid(alpha.1,alpha.1))
    values <- diff.delta(alpha.grid,sigma)
    data = data.frame(x=alpha.grid[,1],y=alpha.grid[,2],z=values)[values>0.1,]    
    p.list[[i]] <- ggplot(data) + geom_point(aes(x=x, y=y, color=z)) + scale_color_gradient(low = "blue", high = "red") + ggtitle(paste("sigma:",para.norm[i,])) + xlim(-5,5) + ylim(-5,5)
}

grid.arrange(grobs=p.list,ncol=2)

## explore the likelihood surface ##
samples.skew.normal <- simu_logskew(m=m,par=alpha2delta(list(cov.func(coord,c(4,1)),alpha.func(c(2,1)))),ncores=ncores)

alpha.1 = seq(-5,5,0.1)
length(alpha.1)
alpha.grid = as.matrix(expand.grid(alpha.1,alpha.1))
alpha.grid.list <- split(alpha.grid,row(alpha.grid))
t0 <- proc.time()
fit.values <- unlist(mclapply(alpha.grid.list,function(x){mean(fit.model(data=samples.skew.normal,loc=coord,init=c(4,1,x),fixed=c(F,F,F,F),thres=0.95,model="logskew",ncores=NULL,lb=lb,ub=ub,bootstrap=FALSE,hessian=FALSE,opt=FALSE,alpha.func=alpha.func,FUN=cov.func))},mc.cores=ncores,mc.set.seed = FALSE))
print(t0 <- proc.time()- t0)
alpha.grid.list[[which.min(unlist(fit.values))]]
para.alpha[idx,]
min(unlist(fit.values))

z=matrix(unlist(fit.values),nrow=length(alpha.1),ncol=length(alpha.1))
p <- plotly::plot_ly(z=~z, type = "surface")
p 


n=3#ceiling(100^(1/ncol(para.alpha)))
alpha.vec = c(-1,0,1)#seq(-2,2,length.out=n)
alpha.vec = matrix(alpha.vec,ncol=ncol(para.alpha),nrow=length(alpha.vec))
alphas.grid = as.matrix(do.call(expand.grid,split(alpha.vec,col(alpha.vec))))
alphas.grid.list <- split(alphas.grid,row(alphas.grid))
fit.values <- unlist(mclapply(alphas.grid.list,function(x){mean(fit.model(data=samples.skew.normal,loc=coord,init=c(1,1,x),fixed=c(F,F,F,F),thres=0.9,model="logskew",ncores=NULL,lb=lb,ub=ub,bootstrap=FALSE,hessian=FALSE,opt=FALSE))},mc.cores=ncores,mc.set.seed = FALSE))
init = c(1,1,alphas.grid.list[[which.min(unlist(fit.values))]])
print(init)
fit.result <- fit.model(data=samples.skew.normal,loc=coord,init=init,fixed=c(F,F,F,F),thres=0.95,model="logskew",ncores=ncores,lb=lb,ub=ub,bootstrap=FALSE,hessian=FALSE,opt=TRUE,method="Nelder-Mead",trace=FALSE,maxit=1000)
fit.result$par
para.alpha[idx,]


## play 
files.list = list.files("data/samples/",pattern="simulation_logskew_comp_*",full.names = TRUE)
samples.skew.normal <- NULL
for(i in 1:length(files.list)){
    load(files.list[i],e<-new.env())
    samples.skew.normal <- lapply(1:length(e$samples.skew.normal),function(x){rbind(samples.skew.normal[[x]],e$samples.skew.normal[[x]])})
    print(i)
}


ub=c(Inf,Inf,rep(Inf,ncol(para.alpha)))
i=6
data = samples.skew.normal[[i]][1:1000,]
data.avg = rowMeans(data)
idx = data.avg > 100 & data.avg < 1000
if(sum(idx)<2){idx = c(which(idx),which.max(data.avg))}
data = data[idx,,drop=FALSE]
fit.logskew.angular[[i]] <- fit.model(data=data,init=e$par.skew.normal[i,],fixed=c(F,F,T,T),loc=diff.mat,thres=nrow(data),FUN=cov.func,alpha.func=alpha.func,model="logskew",lb=lb,ub=ub,ncores=ncores,maxit=1000,trace=TRUE,method="Nelder-Mead",opt=TRUE,hessian=FALSE,basis=basis)
fit.logskew.angular[[i]]$par - e$par.skew.normal[i,]   
sum(idx)




source("code/exponent_functions.R")
set.seed(1)
coord = cbind(1,c(1,5))
par.logskew = alpha2delta(list(vario.func(coord,c(3,1)),c(-2,2)))
par.logskew = alpha2delta(list(matrix(c(1,0.5,0.5,2),2,2),c(-2,2)))
par.logskew = alpha2delta(list(matrix(c(1,0.5,0.5,1),2,2),c(0,0)))
par.logskew = alpha2delta(list(vario.func(coord,c(2,1)),c(-1,0)))

x = c(2,5)
func <- function(x.i){
    x.new <- x;x.new[2] <- x.i
    -V_logskew(x.new,par.logskew,alpha.para=FALSE)
}
grad(func,x[2])
partialV_logskew(x,idx=2,par.logskew,alpha.para=FALSE)

func <- function(x.i){
    x.new <- x;x.new[1] <- x.i
    -V_logskew(x.new,par.logskew,alpha.para=FALSE)
}
grad(func,x[1])
partialV_logskew(x,idx=1,par.logskew,alpha.para=FALSE)

func <- function(x.i){
    exp(-nloglik(par.logskew,x.i,model="logskew"))
}

hcubature(func,lowerLimit = rep(0,2),upperLimit=rep(Inf,2))

