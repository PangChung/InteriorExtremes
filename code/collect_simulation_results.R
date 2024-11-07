#files.list <- list()
    #files.list[[i]] <- list.files(path = "data/simulation_25_25_1000/", pattern = paste0("simulation_study_\\d+_",thres.list[[i]],".RData"), full.names = TRUE, recursive = FALSE)
rm(list=ls())
source("code/exponent_functions.R")
library(ggplot2)
library(gridExtra)
library(tidyr)

## for Pareto process ## 
source("code/exponent_functions.R")
library(ggplot2)
library(gridExtra)
library(tidyr)
idx.file = "_5"
path = paste0("data/simulation_pareto",idx.file,"/")
files.pareto.logskew.1 <- list.files(path=path,pattern="simulation_pareto_logskew_\\d+_1.RData",full.names=TRUE,recursive=FALSE)
files.pareto.logskew.3 <- list.files(path=path,pattern="simulation_pareto_logskew_\\d+_3.RData",full.names=TRUE,recursive=FALSE)
files.pareto.truncT.1 <- list.files(path=path,pattern="simulation_pareto_truncT_\\d+_1.RData",full.names=TRUE,recursive=FALSE)
files.pareto.truncT.3 <- list.files(path=path,pattern="simulation_pareto_truncT_\\d+_3.RData",full.names=TRUE,recursive=FALSE)

collect_results <- function(files){
    load(files,e<-new.env())
    #if(is.null(e$fit.logskew)){fit.results <- e$fit.truncT;n<-ncol(e$par.truncT);n<-c(n-1,n)}else{fit.results <- e$fit.logskew;n<-ncol(e$par.skew.normal);n<-c(n,n)}
    if(is.null(e$fit.logskew)){fit.results <- e$fit.truncT;n<-ncol(e$par.truncT);n<-c(n,n)}else{fit.results <- e$fit.logskew;n<-ncol(e$par.skew.normal);n<-c(n,n)}
    est.1 <- do.call(rbind,lapply(fit.results,function(x){if(is.numeric(x[[1]][[1]])) x[[1]]$par else rep(NA,n[1])}))
    val.1 <- unlist(lapply(fit.results,function(x){if(is.numeric(x[[1]][[1]])) x[[1]]$val else NA}))
    time.1 <- unlist(lapply(fit.results,function(x){if(is.numeric(x[[1]][[1]])) x[[1]]$time[[3]] else NA}))
    est.2 <- do.call(rbind,lapply(fit.results,function(x){if(is.numeric(x[[1]][[1]])) x[[2]]$par else rep(NA,n[2])}))
    val.2 <- unlist(lapply(fit.results,function(x){if(is.numeric(x[[1]][[1]])) x[[2]]$val else NA}))
    time.2 <- unlist(lapply(fit.results,function(x){if(is.numeric(x[[1]][[1]])) x[[2]]$time[[3]] else NA}))
    m = length(val.1)
    if(ncol(est.1)!=ncol(est.2)){if(ncol(est.1)<ncol(est.2)) est.1 <- cbind(est.1,est.2[,ncol(est.1)+1]) else est.2 <- cbind(est.2,est.1[,ncol(est.2)+1])}
    result <- data.frame(case=rep(1:m,2),method=rep(1:2,each=m),time=c(time.1,time.2),val=c(val.1,val.2))
    result <- cbind(result,rbind(est.1,est.2))
    return(result)
}

fit.pareto.logskew.1 <- do.call(rbind,lapply(files.pareto.logskew.1,collect_results))
fit.pareto.logskew.3 <- do.call(rbind,lapply(files.pareto.logskew.3,collect_results))
fit.pareto.truncT.1 <- do.call(rbind,lapply(files.pareto.truncT.1,collect_results))
fit.pareto.truncT.3 <- do.call(rbind,lapply(files.pareto.truncT.3,collect_results))

load(files.pareto.logskew.1[1],e<-new.env());par.logskew <- e$par.skew.normal
load(files.pareto.truncT.1[1],e<-new.env());par.truncT <- e$par.truncT
names(fit.pareto.logskew.1) <- names(fit.pareto.logskew.3) <- c("case","method","time","val","log(hat(lambda))","hat(vartheta)","hat(b)[0]","hat(b)[1]","hat(b)[2]")

fit.pareto.logskew.1[,5] = log(fit.pareto.logskew.1[,5])
fit.pareto.logskew.3[,5] = log(fit.pareto.logskew.3[,5])
fit.pareto.truncT.1[,5] = log(fit.pareto.truncT.1[,5])
fit.pareto.truncT.3[,5] = log(fit.pareto.truncT.3[,5])

names(fit.pareto.truncT.1) <- names(fit.pareto.truncT.3) <- c("case","method","time","val","log(hat(lambda))","hat(vartheta)")#,"deg")

par.logskew <- as.data.frame(par.logskew);par.truncT <- as.data.frame(par.truncT)
names(par.logskew) = names(fit.pareto.logskew.1)[5:9];names(par.truncT) = names(fit.pareto.truncT.1)[5:7]
par.logskew$case = 1:nrow(par.logskew);par.truncT$case = 1:nrow(par.truncT)

levels = c("log(hat(lambda))","hat(vartheta)","hat(b)[1]","hat(b)[2]")
data_true = par.logskew
data_true[,1] = log(data_true[,1])

data_true <- pivot_longer(data_true, cols=levels, names_to = "Variable", values_to = "Value")
data_true$facet = factor(paste0(data_true$Variable),levels=levels)
data_true <- rbind(data_true,data_true)
data_true$method = rep(1:2,each=nrow(data_true)/2)

p.list.logskew <- list()
thres.b = Inf
for(i in 1:2){
    data = subset(fit.pareto.logskew.1,method==i & abs(`hat(b)[1]`) < thres.b & abs(`hat(b)[2]`) < thres.b) #& !(case %in% c(5,7,8,9,11,12)))
    data.true.sub <- subset(data_true,method==i) #& !(case %in% c(5,7,8,9,11,12)))
    print(dim(data));print(dim(subset(fit.pareto.logskew.1,method==i & !is.na(`hat(b)[1]`))))
    data_long <- pivot_longer(data, cols=levels, names_to = "Variable", values_to = "Value")
    data_long$facet = factor(paste0(data_long$Variable),level=levels)
    p <- ggplot(data_long, aes(x = factor(case), y = Value)) + #,fill=factor(method,labels=c("Score","Spectral")))) +
    geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
    geom_point(data=data.true.sub,aes(x=factor(case),y=Value),color="black",size=1,position=position_dodge(width = 1)) +
    facet_wrap(~ facet, scales = "free",nrow=2,ncol=2,labeller = label_parsed) +
    labs(x = "Cases",
            y = "Value",fill="Method") + 
    theme(axis.text = element_text(size = 14),
            axis.title.x = element_text(size = 16),
            strip.text = element_text(size = 16),
            axis.title.y = element_text(size = 16),
            plot.title = element_text(hjust = 0.5, size = 16),
            legend.title = element_text(size = 16)) 
    p.list.logskew[[i]] <- p
}

for(i in 1:2){
        data = subset(fit.pareto.logskew.3,method==i & abs(`hat(b)[1]`) < thres.b) #& abs(`hat(b)[2]`) < thres.b & !(case %in% c(5,7,8,9,11,12)))
        data.true.sub <- subset(data_true,method==i )#& !(case %in% c(5,7,8,9,11,12)))
        print(dim(data));print(dim(subset(fit.pareto.logskew.3,method==i & !is.na(`hat(b)[1]`))))
        data_long <- pivot_longer(data, cols=levels, names_to = "Variable", values_to = "Value")
        data_long$facet = factor(paste0(data_long$Variable),level=levels)
        p <- ggplot(data_long, aes(x = factor(case), y = Value)) + #,fill=factor(method,labels=c("Score","Spectral")))) +
        geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
        geom_point(data=data.true.sub,aes(x=factor(case),y=Value),color="black",size=1,position=position_dodge(width = 1)) +
        facet_wrap(~ facet, scales = "free",nrow=2,ncol=2,labeller = label_parsed) +
        labs(x = "Cases",
                y = "Value",fill="Method") + 
        theme(axis.text = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                strip.text = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                plot.title = element_text(hjust = 0.5, size = 16),
                legend.title = element_text(size = 16)) 
    p.list.logskew[[i+2]] <- p
}

levels = c("log(hat(lambda))","hat(vartheta)")
data_true = par.truncT[,-3]
data_true[,1] = log(data_true[,1])
data_true <- pivot_longer(data_true, cols=levels, names_to = "Variable", values_to = "Value")
data_true$facet = factor(paste0(data_true$Variable),levels=levels)
data_true <- rbind(data_true,data_true)
data_true$method = rep(1:2,each=nrow(data_true)/2)

p.list.truncT <- list()
thres.lambda = Inf
for(i in 1:2){
    data = subset(fit.pareto.truncT.1, method==i & `log(hat(lambda))` < thres.lambda) # & case %in% c(3,4,7,8))
    data.true.sub <- subset(data_true,method==i) #& case %in% c(3,4,7,8))
    data_long <- pivot_longer(data, cols=levels, names_to = "Variable", values_to = "Value")
    data_long$facet = factor(paste0(data_long$Variable),level=levels)
    p <- ggplot(data_long, aes(x = factor(case), y = Value))+ #,fill=factor(method,labels=c("Score","Spectral")))) +
        geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
        geom_point(data=data.true.sub,aes(x=factor(case),y=Value),color="black",size=1,position=position_dodge(width = 1)) +
        facet_wrap(~ facet, scales = "free",nrow=2,ncol=2,labeller = label_parsed) +
        labs(x = "Cases",
                y = "Value",fill="Method") + 
        theme(axis.text = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                strip.text = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                plot.title = element_text(hjust = 0.5, size = 16),
                legend.title = element_text(size = 16))
    p.list.truncT[[i]] <- p
}

for(i in 1:2){
    data = subset(fit.pareto.truncT.3, method==i & `log(hat(lambda))` < thres.lambda) #& case %in% c(3,4,7,8))
    data.true.sub <- subset(data_true,method==i) #& case %in% c(3,4,7,8))
    data_long <- pivot_longer(data, cols=levels, names_to = "Variable", values_to = "Value")
    data_long$facet = factor(paste0(data_long$Variable),level=levels)
    p <- ggplot(data_long, aes(x = factor(case), y = Value))+ #,fill=factor(method,labels=c("Score","Spectral")))) +
        geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
        geom_point(data=data.true.sub,aes(x=factor(case),y=Value),color="black",size=1,position=position_dodge(width = 1)) +
        facet_wrap(~ facet, scales = "free",nrow=2,ncol=2,labeller = label_parsed) +
        labs(x = "Cases",
                y = "Value",fill="Method") + 
        theme(axis.text = element_text(size = 14),
                axis.title.x = element_text(size = 16),
                strip.text = element_text(size = 16),
                axis.title.y = element_text(size = 16),
                plot.title = element_text(hjust = 0.5, size = 16),
                legend.title = element_text(size = 16))
    p.list.truncT[[i+2]] <- p
}

pdf(file=paste0("figures/simulation_pareto_violin_",idx.file,"_",thres.b,"_",thres.lambda,".pdf"),width=14,height = 5,onefile = TRUE)
grid.arrange(grobs=p.list.logskew[1:2],ncol=2,nrow=1)
grid.arrange(grobs=p.list.logskew[3:4],ncol=2,nrow=1)
grid.arrange(grobs=p.list.truncT[1:2],ncol=2,nrow=1)
grid.arrange(grobs=p.list.truncT[3:4],ncol=2,nrow=1)
dev.off()


## for max-stable ##
files.logskew <- list.files(path="data/simulation_logskew/",pattern="simulation_study_comp_\\d+.RData",full.names=TRUE,recursive=FALSE)

collect_results <- function(files){
    load(files,e<-new.env())
    #if(is.null(e$fit.logskew)){fit.results <- e$fit.truncT;n<-ncol(e$par.truncT);n<-c(n-1,n)}else{fit.results <- e$fit.logskew;n<-ncol(e$par.skew.normal);n<-c(n,n)}
    if(is.null(e$fit.logskew)){fit.results <- e$fit.truncT;n<-ncol(e$par.truncT);n<-c(n,n)}else{fit.results <- e$fit.logskew;n<-ncol(e$par.skew.normal);n<-c(n,n)}
    est.1 <- do.call(rbind,lapply(fit.results,function(x){if(is.numeric(x[[1]][[1]])) x$par else rep(NA,n[1])}))
    val.1 <- unlist(lapply(fit.results,function(x){if(is.numeric(x[[1]][[1]])) x$val else NA}))
    time.1 <- unlist(lapply(fit.results,function(x){if(is.numeric(x[[1]][[1]])) x$time[[3]] else NA}))
    result <- data.frame(case=1:length(time.1),time=time.1,val=val.1)
    result <- cbind(result,est.1)
    return(result)
}
fit.logskew <- do.call(rbind,lapply(files.logskew,collect_results))
load(files.logskew[1],e<-new.env());par.logskew <- e$par.skew.normal
names(fit.logskew) <- c("case","time","val","hat(lambda)","hat(vartheta)","hat(b)[0]","hat(b)[1]","hat(b)[2]")

par.logskew <- as.data.frame(par.logskew)
names(par.logskew) = names(fit.logskew)[-c(1:3)]
par.logskew$case = 1:nrow(par.logskew)

levels = c("hat(lambda)","hat(vartheta)","hat(b)[1]","hat(b)[2]")
data_true = par.logskew

data_true <- pivot_longer(data_true, cols=levels, names_to = "Variable", values_to = "Value")
data_true$facet = factor(paste0(data_true$Variable),levels=levels)
data_true <- rbind(data_true,data_true)

p.list.logskew <- list()

data = fit.logskew
data_long <- pivot_longer(data, cols=levels, names_to = "Variable", values_to = "Value")
data_long$facet = factor(paste0(data_long$Variable),level=levels)
p <- ggplot(data_long, aes(x = factor(case), y = Value)) + #,fill=factor(method,labels=c("Score","Spectral")))) +
geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
geom_point(data=data_true,aes(x=factor(case),y=Value),color="black",size=1,position=position_dodge(width = 1)) +
facet_wrap(~ facet, scales = "free",nrow=2,ncol=2,labeller = label_parsed) +
labs(x = "Cases",
        y = "Value") + 
theme(axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        strip.text = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16)) 
p
