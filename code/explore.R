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
para.alpha = rbind(c(0,0),c(-1,-2),c(2,1)) ## slant parameter for skewed norm model ##
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
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
ncores=detectCores()
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
set.seed(init.seed)

## compute the basis ###
centers <- rbind(c(0.5,0.5),c(0.25,0.25),c(0.75,0.75))*d
idx.centers <- apply(centers,1,function(x){which.min(apply(coord,1,function(y){sum((x-y)^2)}))})
basis <- sapply(idx.centers,function(x){y=dnorm(diff.mat[x,],mean=0,sd=d);y=y-mean(y)})
basis <- basis * 1000
summary(basis)

basis <- cbind(bs(coord[,1],degree = 1),coord)
basis <- apply(basis,2,function(x){x-mean(x)})

basis <- sapply(1:(ncol(para.alpha)+1),function(x){y <- rep(0,nrow(coord));y[sample(1:nrow(coord),2)] <- c(-1,1);y})

#basis[,1] <- rep(0,nrow(basis))
# 1: fixing the first one, 2: unfixing the first one
#pdf("figures/delta_basis_1.pdf",width=5*3,height = 5,onefile = TRUE)
pdf("figures/delta_basis_5.pdf",width=5*3,height = 5,onefile = TRUE)
for(idx in 1:nrow(para.alpha)){
    beta1 = alpha2delta(list(cov.func(coord,c(0.5,1)),alpha.func(para.alpha[idx,])))[[2]]
    df = data.frame(x = coord[,1], y = coord[,2], z = beta1)
    p1 <- ggplot(df, aes(x = x, y = y, fill = z)) +
    geom_tile() + ggtitle(paste(c("alpha:",para.alpha[idx,]),collapse = " ")) + 
    scale_fill_gradient2(low = "blue",mid="white" ,high = "red",limits=c(-1,1)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5))
    #p1
    beta2 = alpha2delta(list(cov.func(coord,c(0.5,1)),alpha.func(2*para.alpha[idx,])))[[2]]
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

#data = data.frame(x=alpha.grid[,1],y=alpha.grid[,2],z=unlist(delta.diff))
z=matrix(unlist(delta.diff),length(alpha.1),length(alpha.1),byrow = TRUE)
p <- plotly::plot_ly(z=~z, type = "surface")
p 



## explore the likelihood surface ##
idx = 3
samples.skew.normal <- simu_logskew(m=m,par=alpha2delta(list(cov.func(coord,c(0.5,1)),alphas[,idx])),ncores=ncores)

alpha.1 = seq(-4,4,0.1)
length(alpha.1)
alpha.grid = as.matrix(expand.grid(alpha.1,alpha.1))
alpha.grid.list <- split(alpha.grid,row(alpha.grid))
t0 <- proc.time()
fit.values <- unlist(mclapply(alpha.grid.list,function(x){mean(fit.model(data=samples.skew.normal,loc=coord,init=c(1,1,x),fixed=c(F,F,F,F),thres=0.9,model="logskew",ncores=NULL,lb=lb,ub=ub,bootstrap=FALSE,hessian=FALSE,opt=FALSE))},mc.cores=ncores,mc.set.seed = FALSE))
print(t0 <- proc.time()- t0)
alpha.grid.list[[which.min(unlist(fit.values))]]
para.alpha[idx,]
min(unlist(fit.values))
# Data: volcano is provided by plotly
# Plot
data = data.frame(x=alpha.grid[,1],y=alpha.grid[,2],z=unlist(fit.values))
z=matrix(unlist(fit.values),nrow=length(alpha.1),ncol=length(alpha.1))
p <- plot_ly(z=~z, type = "surface")
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



# # plots using ggplot2 #
# # Load the ggplot2 package
# library(ggplot2)
# library(gridExtra)
# # Create some data
# set.seed(123)
# p.list <- list()
# for(idx.case in 1:27){
# data <- data.frame(x = coord[,1],
#                    y = coord[,2],
#                    z = samples.skew.normal[[idx.case]][1,])

# # Plot a gridded image
# p.list[[idx.case]] <- ggplot(data, aes(x = x, y = y, fill = z)) +
#   geom_tile() +
#   scale_fill_gradient(low = "blue", high = "red") +
#   labs(title = paste(par.skew.normal[idx.case,],collapse = " "), x = "X", y = "Y", fill = "Z") + theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot")
# }

# pdf(file="figures/simulation_samples.pdf",width=10,height = 10,onefile = TRUE)
# do.call(grid.arrange, c(p.list[1:9], ncol = 3,nrow=3))
# do.call(grid.arrange, c(p.list[1:9+9], ncol = 3,nrow=3))
# do.call(grid.arrange, c(p.list[1:9+18], ncol = 3,nrow=3))
# dev.off()

# p1.list <- list()
# idx.center = 313
# ind.idx.center = all.pairs[1,] == idx.center |  all.pairs[2,] == idx.center
# for(idx.case in 1:27){
# tc.logskew.idx.center <- tc.logskew[[idx.case]][ind.idx.center]
# ind.idx = apply(all.pairs[,ind.idx.center],2,function(x) x[x!=idx.center])
# data <- data.frame( x = coord[-idx.center,1],
#                     y = coord[-idx.center,2],
#                     z = tc.logskew.idx.center)

# # plot contours #

# p1.list[[idx.case]] <- ggplot(data, aes(x = x, y = y, z = z)) +
#     geom_contour(aes(colour = ..level..)) +
#     scale_colour_gradient(low = "blue", high = "red") +
#     theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
#     labs(title = paste("Bivariate Extremal Coef.:",paste(par.skew.normal[idx.case,],collapse = " ")), x = "X", y = "Y", colour = "Z")

# }

# pdf(file="figures/simulation_samples_extcoef_contours.pdf",width=10,height = 8,onefile = TRUE)
# do.call(grid.arrange, c(p1.list[1:9], ncol = 3,nrow=3))
# do.call(grid.arrange, c(p1.list[1:9+9], ncol = 3,nrow=3))
# do.call(grid.arrange, c(p1.list[1:9+18], ncol = 3,nrow=3))
# dev.off()

# save(p1.list,p.list,file="data/simulation_samples_plots.RData")

# # system("say \'your program has finished\' ")

# idx.case = 11
# fit.logskew.angular[[idx.case]]$par
# par.skew.normal[idx.case,]

