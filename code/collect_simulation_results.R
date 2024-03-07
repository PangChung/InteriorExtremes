#files.list <- list()
    #files.list[[i]] <- list.files(path = "data/simulation_25_25_1000/", pattern = paste0("simulation_study_\\d+_",thres.list[[i]],".RData"), full.names = TRUE, recursive = FALSE)
rm(list=ls())
source("code/exponent_functions.R")
library(ggplot2)
library(gridExtra)
library(tidyr)

idx.file = 13
files.list <- list.files(path=paste0("data/simulation_",idx.file,"_1000"),pattern="simulation_study_logskew_\\d+_1000_1.RData",full.names=TRUE,recursive=FALSE)
thres.list = c(0.95,0.9)
load(files.list[[1]],e<-new.env())
par.skew.normal = e$par.skew.normal
rm(e)
extract_results <- function(files){
    fit.results <- list()
    for(i in 1:length(files)){
        load(files[[i]],e<-new.env())
        #fit.logskew.angular2 <- lapply(e$fit.logskew.angular2,function(x){values = lapply(x,function(x1){unlist(lapply(x1, function(x2){x2$value}))});return(lapply(1:length(values),function(i){x[[i]][[which.min(values[[i]])]]}))})
        #fit.results[[i]] <- fit.logskew.angular2
        fit.results[[i]] <- e$fit.logskew.angular
    }
    par.skew <- e$par.skew.normal
    n1 = nrow(par.skew);n2 = length(e$fit.logskew.angular[[1]])
    est.mat.list <- lapply(1:n2,function(x){  lapply(1:n1,function(x1){list()})})
    for(i in 1:n1){
        for(j in 1:n2){
            est.mat.list[[j]][[i]] <- matrix(unlist(lapply(fit.results,function(x){x[[i]][[j]]$par[1:4]})),ncol=ncol(par.skew),byrow=TRUE)
        }
    }
    return(list(est.mat.list,par.skew))
}

mse.max = matrix(NA,nrow=length(files.list),ncol=24)
for(k in 1:length(files.list)){
    load(files.list[k],e<-new.env())    
    fit.result <- lapply(1:nrow(par.skew.normal),function(i){values = lapply(1:2,function(j){matrix(unlist(lapply(e$fit.logskew.angular2[[i]][[j]], function(x2){x2$par[1:4]-par.skew.normal[i,]})),ncol=4,byrow=TRUE)})})
    fit.result <- unlist(lapply(fit.result,function(x){lapply(x,function(x1){mse=apply(abs(x1[,3:4]),1,mean);min(mse)})}))
    mse.max[k,] <- fit.result
    print(k)
}

boxplot(mse.max)
summary(mse.max)
idx= 1
max(mse.max[,idx])
#apply(mse.max,2,function(x) which(x>1))
error.idx = which.max(mse.max[,idx])
load(files.list[error.idx],e<-new.env())
e$fit.logskew.angular2[[(idx-1) %/% 2 + 1]][[(idx-1) %% 2 + 1]]
e$fit.logskew.angular[[(idx-1) %/% 2 + 1]][[(idx-1) %% 2 + 1]]
par.skew.normal[(idx-1) %/% 2 + 1,]

est.mat.list <- extract_results(files.list)
par.skew.normal = est.mat.list[[2]];est.mat.list = est.mat.list[[1]]
par.skew.normal = as.data.frame(par.skew.normal)
save.image(file=paste0("data/simulation_study_logskew_results_",idx.file,"_1000.RData"))


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
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("Threshold: ",thres.list[idx.thres],"%"," with 1000 replicates")) + geom_point(data=data.true,aes(x=Variable, y = Value), color = "red") + ylim(max(-10,min(data_long$Value)),min(10,max(data_long$Value)))
        p.list[[idx.thres]][[idx.case]] <- p
    }
}

pdf(file=paste0("figures/simulation_est_boxplots_",idx.file,"_1000.pdf"),width=4*n2,height = 5,onefile = TRUE)
for(idx.case in 1:n1){
    do.call(grid.arrange, c(lapply(p.list,function(x){x[[idx.case]]}), ncol = n2,nrow=1))
}
dev.off()


idx.file = "8_8"
files.list <- list.files(path=paste0("data/simulation2_",idx.file,"_5000"),pattern="simulation_study_truncT_\\d+",full.names=TRUE,recursive=FALSE)
thres.list = c(0.95,0.9)

extract_results_truncT <- function(files){
    fit.results <- list()
    for(i in 1:length(files)){
        load(files[[i]],e<-new.env())
        fit.results[[i]] <- e$fit.truncT.angular
    }
    par.truncT <- e$par.truncT
    n1 = nrow(par.truncT);n2 = length(e$fit.truncT.angular[[1]])
    est.mat.list <- lapply(1:n2,function(x){  lapply(1:n1,function(x1){list()})})
    for(i in 1:n1){
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
    for(idx.case in 1:n1){
        data = as.data.frame(est.mat.list[[idx.thres]][[idx.case]])
        data.true <- pivot_longer(par.truncT[idx.case,], everything(), names_to = "Variable", values_to = "Value")
        data.true$Variable <- paste0("V",1:3)
        data_long <- pivot_longer(data, everything(), names_to = "Variable", values_to = "Value")
        
        p<- ggplot(data_long, aes(x = Variable, y = Value)) +
        geom_boxplot() + scale_x_discrete(labels=variable.names) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("Threshold: ",thres.list[idx.thres],"%"," with 5000 replicates")) + geom_point(data=data.true,aes(x=Variable, y = Value), color = "red") + ylim(0,3)
        p.list[[idx.thres]][[idx.case]] <- p
    }
}

pdf(file="figures/simulation_est_boxplots_truncT2_3_5000.pdf",width=4*n2,height = 5,onefile = TRUE)
for(idx.case in 1:n1){
    do.call(grid.arrange, c(lapply(p.list,function(x){x[[idx.case]]}), ncol = n2,nrow=1))
}
dev.off()


