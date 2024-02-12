args <- commandArgs(TRUE)
# settings 
id = 1
ncores=detectCores()
d <- 25 ## 10 * 10 grid on [0,1]^2
m <- 200 ## number of samples
m <- 1000 ## number of samples
set.seed(1342342)
coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
para.range = c(1,2)#c(0.5,1,2) ## range for the correlation function ##
para.nu = c(0.5,1)#c(0.5,1,1.5) ## smoothness parameter for the correlation function ##
para.alpha = rbind(c(0,0,0),c(-1,2,3),c(-2,-1,4)) ## slant parameter for skewed norm model ##
para.deg = 2 ## degree of the freedom for the truncated t model ##
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
thres = 0.9
computer = "local"
# loading library and setting path
for (arg in args) eval(parse(text = arg))
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
switch(computer,
    "ws" = {DataPath<-"~/Desktop/InteriorExtremes/"},
    "hpc" = {DataPath<-"/srv/scratch/z3536974/";.libPaths("../src")},
    "local" = {DataPath<-"~/Documents/Github/InteriorExtremes/"}
)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")

id = 1
ncores=detectCores()
d <- 25 ## 10 * 10 grid on [0,1]^2
m <- 1000 ## number of samples
set.seed(1342342)
coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
para.range = c(1,2)#c(0.5,1,2) ## range for the correlation function ##
para.nu = c(0.5,1)#c(0.5,1,1.5) ## smoothness parameter for the correlation function ##
para.alpha = rbind(c(0,0,0),c(-1,2,3),c(-2,-1,4)) ## slant parameter for skewed norm model ##
para.deg = 2 ## degree of the freedom for the truncated t model ##
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
for (arg in args) eval(parse(text = arg))
file2save = paste0("/srv/scratch/z3536974/data/simulation_study_",id,"_",thres*100,"_",m,".RData")
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
set.seed(init.seed)

########################################################################
### simulation study for the log-skew normal based max-stable process ##
########################################################################
par.skew.normal <- as.matrix(expand.grid(para.range,para.nu,1:3))
par.skew.normal <- cbind(par.skew.normal[,-3],para.alpha[par.skew.normal[,3],]);colnames(par.skew.normal) <- NULL
samples.skew.normal <- list()
par.skew.list <- list()
ec.logskew <- list()
tc.logskew <- list()
fit.logskew.angular <- list()
for(i in 1:nrow(par.skew.normal)){
    par.skew.list[[i]] <- list(sigma=cov.func(coord,par.skew.normal[i,1:2]),alpha=alpha.func(coord,par.skew.normal[i,-c(1:2)]))
    samples.skew.normal[[i]] <- simu_logskew(m=m,par=alpha2delta(par.skew.list[[i]]),ncores=ncores)
    # ec.logskew[[i]] <- lapply(all.pairs,empirical_extcoef,data=samples.skew.normal[[i]])
    # tc.logskew[[i]] <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=alpha2delta(par.skew.list[[i]]),model="logskew1"),mc.cores=ncores,mc.set.seed=FALSE)
    tryCatch(fit.logskew.angular[[i]] <- fit.model(data=samples.skew.normal[[i]],loc=coord,init=par.skew.normal[i,],fixed=c(F,F,F,F,F),thres=thres,model="logskew",ncores=ncores,maxit=100,lb=c(0.01,0.01,-Inf,-Inf,-Inf),ub=c(10,2.0,Inf,Inf,Inf),bootstrap=FALSE,hessian=TRUE,opt=TRUE),
             error=function(e){print(e);fit.logskew.angular[[i]] <- fit.model(data=samples.skew.normal[[i]],loc=coord,init=c(0.5,1.5,-0.1,0.1,0.1),fixed=c(F,F,F,F,F),thres=thres,model="logskew",ncores=ncores,maxit=100,lb=c(0.01,0.01,-Inf,-Inf,-Inf),ub=c(10,2.0,Inf,Inf,Inf),bootstrap=FALSE,hessian=TRUE,opt=TRUE)})    
}

save(fit.logskew.angular,par.skew.normal,file=file2save)

## trainning the model with the log-skew normal based max-stable process ##
# idx.case = 11
# fit.logskew.comp <- MCLE(data=samples.skew.normal[[idx.case]],init=c(0.5,1,0,0,0),fixed=c(F,F,F,F,F),loc=coord,FUN=cov.func,index=all.pairs[,sample(1:ncol(all.pairs),1000,replace=FALSE)],ncores=ncores,maxit=200,model="logskew",lb=c(0.1,0.1,-Inf,-Inf,-Inf),ub=c(10,2.5,Inf,Inf,Inf), alpha.func=alpha.func,hessian=TRUE)

# vecchia.seq <- sample(1:nrow(coord),size=nrow(coord),replace=FALSE)
# neighbours.mat <- sapply(1:nrow(coord),FUN=neighbours,vecchia.seq=vecchia.seq,
# 					q=3,loc=diff.mat)
# fit.logskew.vecchia <- MVLE(data=samples.skew.normal[[idx.case]],init=c(0.5,1,0,0,0),fixed=c(F,F,F,F,F),loc=coord,FUN=cov.func,vecchia.seq=vecchia.seq,neighbours = neighbours.mat,alpha.func=alpha.func,maxit=200,model="logskew",lb=c(0.1,0.1,-Inf,-Inf,-Inf),ub=c(10,2,Inf,Inf,Inf),ncores=ncores)

# # trainning the model with angular density 
# fit.logskew.angular <- list()
# for(i in 1:nrow(par.skew.normal)){
#     tryCatch(fit.logskew.angular[[i]] <- fit.model(data=samples.skew.normal[[i]],loc=coord,init=c(1,1.5,-0.1,0.1,0.1),fixed=c(F,F,F,F,F),thres=0.9,model="logskew",ncores=ncores,maxit=100,lb=c(0.01,0.01,-Inf,-Inf,-Inf),ub=c(10,2.0,Inf,Inf,Inf),bootstrap=FALSE,hessian=TRUE,opt=TRUE),
#              error=function(e){print(e);fit.logskew.angular[[i]] <- fit.model(data=samples.skew.normal[[i]],loc=coord,init=c(0.5,1.5,-0.1,0.1,0.1),fixed=c(F,F,F,F,F),thres=0.9,model="logskew",ncores=ncores,maxit=100,lb=c(0.01,0.01,-Inf,-Inf,-Inf),ub=c(10,2.0,Inf,Inf,Inf),bootstrap=FALSE,hessian=TRUE,opt=TRUE)})
# }

# save(samples.skew.normal,par.skew.list,ec.logskew,tc.logskew,fit.logskew.angular,par.skew.normal,neighbours.mat,vecchia.seq,file="data/simulation_study_logskew.RData")

# ## high dimensional case ## 
# d <- 25 ## 10 * 10 grid on [0,1]^2
# m <- 1000 ## number of samples
# set.seed(1342342)
# coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
# diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
# diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
# samples.skew.normal <- list()
# par.skew.list <- list()
# ec.logskew <- list()
# tc.logskew <- list()
# all.pairs = combn(1:nrow(coord),2)
# all.pairs.list = split(all.pairs,col(all.pairs))

# for(i in 1:nrow(par.skew.normal)){
#     par.skew.list[[i]] <- list(sigma=cov.func(coord,par.skew.normal[i,1:2]),alpha=alpha.func(coord,par.skew.normal[i,-c(1:2)]))
#     system.time(samples.skew.normal[[i]] <- simu_logskew(m=m,par=alpha2delta(par.skew.list[[i]]),ncores=ncores))
#     ec.logskew[[i]] <- mclapply(all.pairs.list,empirical_extcoef,data=samples.skew.normal[[i]],mc.cores=ncores,mc.set.seed=FALSE)
#     tc.logskew[[i]] <- mcmapply(true_extcoef,all.pairs.list,MoreArgs=list(par=alpha2delta(par.skew.list[[i]]),model="logskew1"),mc.cores=ncores,mc.set.seed=FALSE)
#     print(i)
# }

# save(samples.skew.normal,par.skew.list,ec.logskew,tc.logskew,fit.logskew.angular,par.skew.normal,neighbours.mat,vecchia.seq,file="data/simulation_study_logskew_high.RData")

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

