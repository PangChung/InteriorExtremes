#files.list <- list()
    #files.list[[i]] <- list.files(path = "data/simulation_25_25_1000/", pattern = paste0("simulation_study_\\d+_",thres.list[[i]],".RData"), full.names = TRUE, recursive = FALSE)
source("code/exponent_functions.R")
idx.file = 9
files.list <- list.files(path=paste0("data/simulation_",idx.file,"_1000"),pattern="simulation_study_logskew_\\d+_1000.RData",full.names=TRUE,recursive=FALSE)
thres.list = c(0.95,0.9)

extract_results <- function(files){
    fit.results <- list()
    for(i in 1:length(files)){
        load(files[[i]],e<-new.env())
        fit.results[[i]] <- e$fit.logskew.angular
    }
    par.skew <- e$par.skew.normal
    n1 = nrow(par.skew);n2 = length(e$fit.logskew.angular[[1]])
    est.mat.list <- lapply(1:n2,function(x){  lapply(1:n1,function(x1){list()})})
    for(i in 1:n1){
        for(j in 1:n2){
            est.mat.list[[j]][[i]] <- matrix(unlist(lapply(fit.results,function(x){x[[i]][[j]]$par})),ncol=ncol(par.skew),byrow=TRUE)
        }
    }
    return(list(est.mat.list,par.skew))
}

est.mat.list <- extract_results(files.list)
par.skew.normal = est.mat.list[[2]];est.mat.list = est.mat.list[[1]]
par.skew.normal = as.data.frame(par.skew.normal)
save(est.mat.list,files.list,file=paste0("data/simulation_study_logskew_results_",idx.file,"_1000.RData"))

library(ggplot2)
library(gridExtra)
library(tidyr)

variable.names <- c(expression(lambda), expression(nu), expression(alpha[1]), expression(alpha[2]), expression(alpha[3]))

n1 = length(est.mat.list[[1]]);n2 = length(est.mat.list)
p.list = create_lists(c(n2,n1))
for(idx.thres in 1:n2){
    for(idx.case in 1:n1){
        data = as.data.frame(est.mat.list[[idx.thres]][[idx.case]])
        data.true <- pivot_longer(par.skew.normal[idx.case,], everything(), names_to = "Variable", values_to = "Value")
        data_long <- pivot_longer(data, everything(), names_to = "Variable", values_to = "Value")
        
        p<- ggplot(data_long, aes(x = Variable, y = Value)) +
        geom_boxplot() + scale_x_discrete(labels=variable.names) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("Threshold: ",thres.list[idx.thres],"%"," with 1000 replicates")) + geom_point(data=data.true,aes(x=Variable, y = Value), color = "red") + ylim(max(-5,min(data_long$Value)),min(5,max(data_long$Value)))
        p.list[[idx.thres]][[idx.case]] <- p
    }
}

pdf(file=paste0("figures/simulation_est_boxplots_",idx.file,"_1000.pdf"),width=4*n2,height = 5,onefile = TRUE)
for(idx.case in 1:n1){
    do.call(grid.arrange, c(lapply(p.list,function(x){x[[idx.case]]}), ncol = n2,nrow=1))
}
dev.off()


extract_results_truncT <- function(files){
    fit.results <- list()
    for(i in 1:length(files)){
        load(files[[i]],e<-new.env())
        fit.results[[i]] <- e$fit.truncT.angular
    }
    par.truncT <- e$par.truncT
    n1 = nrow(par.truncT);n2 = length(e$fit.truncT.angular[[7]])
    est.mat.list <- lapply(1:n2,function(x){  lapply(1:n1,function(x1){list()})})
    for(i in 6:n1){
        for(j in 1:n2){
            est.mat.list[[j]][[i]] <- matrix(unlist(lapply(fit.results,function(x){x[[i]][[j]]$par})),ncol=3,byrow=TRUE)
        }
    }
    return(list(est.mat.list,par.truncT))
}

est.mat.list <- extract_results_truncT(files.list)
par.truncT = est.mat.list[[2]];est.mat.list = est.mat.list[[1]]
par.truncT = as.data.frame(par.truncT)
variable.names <- c(expression(lambda), expression(nu), expression(deg))

n1 = length(est.mat.list[[1]]);n2 = length(est.mat.list)
p.list = create_lists(c(n2,n1))
for(idx.thres in 1:n2){
    for(idx.case in 6:n1){
        data = as.data.frame(est.mat.list[[idx.thres]][[idx.case]])
        data.true <- pivot_longer(par.truncT[idx.case,], everything(), names_to = "Variable", values_to = "Value")
        data.true$Variable <- paste0("V",1:3)
        data_long <- pivot_longer(data, everything(), names_to = "Variable", values_to = "Value")
        
        p<- ggplot(data_long, aes(x = Variable, y = Value)) +
        geom_boxplot() + scale_x_discrete(labels=variable.names) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("Threshold: ",thres.list[idx.thres],"%"," with 1000 replicates")) + geom_point(data=data.true,aes(x=Variable, y = Value), color = "red") + ylim(max(-5,min(data_long$Value)),min(5,max(data_long$Value)))
        p.list[[idx.thres]][[idx.case]] <- p
    }
}

pdf(file="figures/simulation_est_boxplots_truncT_2_1000.pdf",width=4*n2,height = 5,onefile = TRUE)
for(idx.case in 6:n1){
    do.call(grid.arrange, c(lapply(p.list,function(x){x[[idx.case]]}), ncol = n2,nrow=1))
}
dev.off()

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


