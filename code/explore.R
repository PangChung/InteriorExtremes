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
coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
para.range = c(1,2) #c(0.5,1,2) ## range for the correlation function ##      
para.nu = 1 #c(0.5,1,1.5) ## smoothness parameter for the correlation function ##
para.alpha = rbind(c(0,0,0),c(-1,-2,-3),c(-2,-1,4),c(2,1,4)) ## slant parameter for skewed norm model ##
#para.alpha = rbind(c(0,0),c(-1,-2),c(-2,-1),c(2,1),c(-1,2))
para.deg = c(2,3) ## degree of the freedom for the truncated t model ##
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
thres = c(0.95,0.9)
lb=c(0.01,0.01,rep(-Inf,nrow(para.alpha)))
ub=c(10,2.0,rep(Inf,ncol(para.alpha)))
# loading library and setting path
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
#library(gridExtra)
library(partitions)
#library(ggplot2)
library(Rfast)
library(matrixStats)
library(splines)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
ncores=detectCores()
file2save = paste0(DataPath,"data/simulation_study_",model,"_",id,"_",m,".RData")
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
set.seed(init.seed)

##compute the basis ###
centers <- rbind(c(0.5,0.5),c(0.25,0.75),c(0.25,0.25),c(0.75,0.75))
#centers <- rbind(c(0.25,0.25),c(0.5,0.5),c(0.75,0.75))
idx.centers <- apply(centers,1,function(x){which.min(apply(coord,1,function(y){sum((x-y)^2)}))})
basis <- sapply(idx.centers,function(x){ y=dnorm(diff.mat[x,],mean=0,sd=0.125);y=y-mean(y) })

## plot the basis functions
alphas = apply(para.alpha,1,alpha.func)

idx=3
beta = alpha2delta(list(cov.func(coord,c(0.5,1)),alphas[,idx]))[[2]]
df = data.frame(x = coord[,1], y = coord[,2], z = beta)
library(ggplot2)
p <- ggplot(df, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()
p

samples.skew.normal <- simu_logskew(m=m,par=alpha2delta(list(cov.func(coord,c(0.5,1)),alphas[,idx])),ncores=ncores)

alpha.1 = seq(-4,4,0.3)
length(alpha.1)
alpha.grid = as.matrix(expand.grid(alpha.1,alpha.1,alpha.1))
alpha.grid.list <- split(alpha.grid,row(alpha.grid))


t0 <- proc.time()
fit.values <- unlist(mclapply(alpha.grid.list,function(x){mean(fit.model(data=samples.skew.normal,loc=coord,init=c(0.5,1,x),fixed=c(F,F,F,F),thres=0.9,model="logskew",ncores=NULL,lb=lb,ub=ub,bootstrap=FALSE,hessian=FALSE,opt=FALSE))},mc.cores=ncores,mc.set.seed = FALSE))
print(t0 <- proc.time()- t0)

init = c(1,1,0.5,-0.5,0.5)
fit.values <- fit.model(data=samples.skew.normal,loc=coord,init=init,fixed=c(F,F,F,F,F),thres=0.9,model="logskew",ncores=ncores,lb=lb,ub=ub,bootstrap=FALSE,hessian=FALSE,opt=TRUE,method="Nelder-Mead",trace=TRUE,maxit=500)
# Library
library(plotly)

# Data: volcano is provided by plotly

# Plot
data = data.frame(x=alpha.grid[,1],y=alpha.grid[,2],z=unlist(fit.values))
z=matrix(unlist(fit.values),nrow=length(alpha.1),ncol=length(alpha.1))
p <- plot_ly(z=~z, type = "surface")
p 


alpha.grid.list[[which.min(unlist(fit.values))]]
para.alpha[idx,]
