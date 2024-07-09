#files.list <- list()
    #files.list[[i]] <- list.files(path = "data/simulation_25_25_1000/", pattern = paste0("simulation_study_\\d+_",thres.list[[i]],".RData"), full.names = TRUE, recursive = FALSE)
rm(list=ls())
source("code/exponent_functions.R")
library(ggplot2)
library(gridExtra)
library(tidyr)

files.list1 <- list.files(path="data/simulation_comp/",pattern=paste0("simulation_comp2_logskew_\\d+_500_1.RData"),full.names=TRUE,recursive=FALSE)
files.list2 <- list.files(path="data/simulation_comp2/",pattern=paste0("simulation_comp2_2_logskew_\\d+_500_1.RData"),full.names=TRUE,recursive=FALSE)
load(files.list2[[1]],e<-new.env())
par.skew.normal = e$par.skew.normal[,-3]
n1 = nrow(par.skew.normal);n2 = 1 #length(e$fit.logskew.angular[[1]])
extract_results <- function(files,comp=FALSE){
    fit.results <- list()
    for(i in 1:length(files)){
        load(files[[i]],e<-new.env())
        if(!comp){fit.results[[i]] <- e$fit.logskew.angular}else{fit.results[[i]] <- e$fit.logskew.comp}
    }
    est.mat.list <- create_lists(c(n2,n1))
    time.used.list <- create_lists(c(n2,n1))
    for(i in 1:n1){
        for(j in 1:n2){
            est.mat.list[[j]][[i]] <- matrix(unlist(lapply(fit.results,function(x){x[[i]]$par})),ncol=length(fit.results[[1]][[1]]$par),byrow=TRUE)       
            time.used.list[[j]][[i]] <- unlist(lapply(fit.results,function(x){x[[i]]$time[1]}))
        }
    }
    return(list(est.mat.list,time.used.list))
}

results.list1.angular <- extract_results(files.list1,comp=F)
results.list1.comp <- extract_results(files.list1,comp=T)

results.list2.angular <- extract_results(files.list2,comp=F)
results.list2.comp <- extract_results(files.list2,comp=T)

data.est  <- data.frame(results.list1.angular[[1]][[1]][[2]],method=1,step=1,case=1)
data.est <- rbind(data.est,data.frame(results.list1.angular[[1]][[1]][[2]],method=1,step=1,case=2))
data.est <- rbind(data.est,data.frame(results.list1.comp[[1]][[1]][[1]],method=2,step=1,case=1))
data.est <- rbind(data.est,data.frame(results.list1.comp[[1]][[1]][[2]],method=2,step=1,case=2))

data.est <- rbind(data.est,data.frame(results.list2.angular[[1]][[1]][[1]],method=1,step=2,case=1))
data.est <- rbind(data.est,data.frame(results.list2.angular[[1]][[1]][[2]],method=1,step=2,case=2))
data.est <- rbind(data.est,data.frame(results.list2.comp[[1]][[1]][[1]],method=2,step=2,case=1))
data.est <- rbind(data.est,data.frame(results.list2.comp[[1]][[1]][[2]],method=2,step=2,case=2))
data.est <- data.est[data.est[,3]>0,]
data.est[,3:5] <- data.est[,3:5]/data.est[,3]
names(data.est) <- c("lambda","nu","b0","b1","b2","method","step","case")

save.image("data/simulation_comp_results.RData")

par.skew.normal = as.data.frame(par.skew.normal)

save.image(file=paste0("data/simulation_study_logskew_results_",idx.file,"_",basis.idx,".RData"))

## for skewed BR ##
n1 = length(est.mat.list[[1]]);n2 = length(est.mat.list)
p.list = create_lists(c(n2,n1))
for(idx.thres in 1:n2){
    for(idx.case in 1:n1){
        variable.names = c(bquote(lambda==.(par.skew.normal[idx.case,1])),bquote(nu==.(par.skew.normal[idx.case,2])),bquote(alpha[1]==.(par.skew.normal[idx.case,3])),bquote(alpha[2]==.(par.skew.normal[idx.case,4])))
        data = matrix(unlist(apply(est.mat.list[[idx.thres]][[idx.case]],1,function(x){x- par.skew.normal[idx.case,]})),ncol=4,byrow = TRUE)
        data = as.data.frame(data)
        data_long <- pivot_longer(data, everything(), names_to = "Variable", values_to = "Value")
        p<- ggplot(data_long, aes(x = Variable, y = Value)) + 
        geom_boxplot() + scale_x_discrete(labels=variable.names)  + #ylim(-5,5) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5,size=10),plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("Threshold: 50"," with 2000 replicates")) + 
        geom_hline(yintercept = 0, linetype="dashed", color = "red")
        p.list[[idx.thres]][[idx.case]] <- p
    }
}

pdf(file=paste0("figures/simulation_est_boxplots_",idx.file,"_",basis.idx,".pdf"),width=4*4,height = 5*3,onefile = TRUE)
grid.arrange(grobs=p.list[[1]],ncol=4)
dev.off()

## for BR ##
idx.file = "vario_50";basis.idx="BR"
load(paste0("data/simulation_study_logskew_results_",idx.file,"_",basis.idx,".RData"))
n1 = length(est.mat.list[[1]]);n2 = length(est.mat.list)
p.list = create_lists(c(n2,n1))
for(idx.thres in 1:n2){
    for(idx.case in 1:n1){
        variable.names = c(bquote(lambda==.(par.skew.normal[idx.case,1])),bquote(nu==.(par.skew.normal[idx.case,2])))
        data = matrix(unlist(apply(est.mat.list[[idx.thres]][[idx.case]],1,function(x){x- par.skew.normal[idx.case,]})),ncol=4,byrow = TRUE)
        data = as.data.frame(data[,1:2])
        data_long <- pivot_longer(data, everything(), names_to = "Variable", values_to = "Value")
        p<- ggplot(data_long, aes(x = Variable, y = Value)) +
        geom_boxplot() + scale_x_discrete(labels=variable.names) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5,size=10),plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("Threshold: 50"," with 500 replicates")) + 
        geom_hline(yintercept = 0, linetype="dashed", color = "red")
        p.list[[idx.thres]][[idx.case]] <- p
    }
}

pdf(file=paste0("figures/simulation_est_boxplots_",idx.file,"_",basis.idx,".pdf"),width=5*3,height = 5*2,onefile = TRUE)
grid.arrange(grobs=p.list[[1]],ncol=3,nrow=2)
dev.off()


idx.file = "final"
files.list <- list.files(path=paste0("data/simulation_",idx.file),pattern="simulation_study_truncT_\\d+_1000_1.RData",full.names=TRUE,recursive=FALSE)
thres.list = c(0.98,0.95,0.9)

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
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("Threshold: ",thres.list[idx.thres],"%"," with 1000 replicates")) + geom_point(data=data.true,aes(x=Variable, y = Value), color = "red") + ylim(0,10)
        p.list[[idx.thres]][[idx.case]] <- p
    }
}

pdf(file="figures/simulation_est_boxplots_truncT_final_2000.pdf",width=4*n2,height = 5,onefile = TRUE)
for(idx.case in 1:n1){
    do.call(grid.arrange, c(lapply(p.list,function(x){x[[idx.case]]}), ncol = n2,nrow=1))
}
dev.off()

save(p.list,est.mat.list,par.truncT,variable.names,file=paste0("data/simulation_study_truncT_results_",idx.file,"_2000.RData"))

