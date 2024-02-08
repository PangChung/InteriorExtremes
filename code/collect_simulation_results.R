thres.list = c(80,85,90,95,99)
files.list <- list()
files.list1000 <- list()
for(i in 1:length(thres.list)){
    files.list[[i]] <- list.files(path = "data", pattern = paste0("simulation_study_\\d+_",thres.list[[i]],".RData"), full.names = TRUE, recursive = FALSE)
    files.list1000[[i]] <- list.files(path = "data", pattern = paste0("simulation_study_\\d+_",thres.list[[i]],"_1000.RData"), full.names = TRUE, recursive = FALSE)
}

extract_results <- function(files){
    fit.results <- list()
    for(i in 1:length(files)){
        load(files[[i]],e<-new.env())
        fit.results[[i]] <- e$fit.logskew.angular
    }
    est.mat.list <- list()
    for(i in 1:27){
        est.mat.list[[i]] <- matrix(unlist(lapply(fit.results,function(x){x[[i]]$par})),ncol=5,byrow=TRUE)
    }
    return(est.mat.list)
}



est.mat.list <- lapply(files.list,extract_results)
est.mat.list1000 <- lapply(files.list1000,extract_results)

save(est.mat.list,est.mat.list1000,files.list,files.list1000,file="data/simulation_study_logskew_results.RData")

library(ggplot2)
library(gridExtra)
library(tidyr)
para.range = c(0.5,1,2) ## range for the correlation function ##
para.nu = c(0.5,1,1.5) ## smoothness parameter for the correlation function ##
para.alpha = rbind(c(0,0,0),c(-1,2,3),c(-2,-1,4)) ## slant parameter for skewed norm model ##
par.skew.normal <- as.matrix(expand.grid(para.range,para.nu,1:3))
par.skew.normal <- cbind(par.skew.normal[,-3],para.alpha[par.skew.normal[,3],]);colnames(par.skew.normal) <- NULL
par.skew.normal <- as.data.frame(par.skew.normal)

variable.names <- c(expression(lambda), expression(nu), expression(alpha[1]), expression(alpha[2]), expression(alpha[3]))

p.list <- list(list(),list(),list(),list(),list())
for(idx.thres in 1:5){
    for(idx.case in 1:27){
        data = as.data.frame(est.mat.list[[idx.thres]][[idx.case]])
        data.true <- pivot_longer(par.skew.normal[idx.case,], everything(), names_to = "Variable", values_to = "Value")
        data_long <- pivot_longer(data, everything(), names_to = "Variable", values_to = "Value")
        p<- ggplot(data_long, aes(x = Variable, y = Value)) +
        geom_boxplot() + scale_x_discrete(labels=variable.names) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("Threshold: ",thres.list[idx.thres],"%"," with 10000 replicates")) + geom_point(data=data.true,aes(x=Variable, y = Value), color = "red") + ylim(max(-10,min(data_long$Value)),min(10,max(data_long$Value)))
        p.list[[idx.thres]][[idx.case]] <- p
    }
}

p.list1000 <- list(list(),list(),list(),list(),list())
for(idx.thres in 1:5){
    for(idx.case in 1:27){
        data = as.data.frame(est.mat.list1000[[idx.thres]][[idx.case]])
        data.true <- pivot_longer(par.skew.normal[idx.case,], everything(), names_to = "Variable", values_to = "Value")
        data_long <- pivot_longer(data, everything(), names_to = "Variable", values_to = "Value")
        
        p<- ggplot(data_long, aes(x = Variable, y = Value)) +
        geom_boxplot() + scale_x_discrete(labels=variable.names) +
        theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + ggtitle(paste0("Threshold: ",thres.list[idx.thres],"%"," with 1000 replicates")) + geom_point(data=data.true,aes(x=Variable, y = Value), color = "red") + ylim(max(-10,min(data_long$Value)),min(10,max(data_long$Value)))
        p.list1000[[idx.thres]][[idx.case]] <- p
    }
}


pdf(file="figures/simulation_est_boxplots.pdf",width=4*5,height = 5*2,onefile = TRUE)
for(idx.case in 1:27){
    do.call(grid.arrange, c(lapply(c(p.list,p.list1000),function(x){x[[idx.case]]}), ncol = 5,nrow=2))
}
dev.off()

