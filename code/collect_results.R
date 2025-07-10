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
names(fit.pareto.logskew.1) <- names(fit.pareto.logskew.3) <- c("case","method","time","val","hat(lambda)","hat(vartheta)","hat(b)[0]","hat(b)[1]","hat(b)[2]")

# fit.pareto.logskew.1[,5] = log(fit.pareto.logskew.1[,5])
# fit.pareto.logskew.3[,5] = log(fit.pareto.logskew.3[,5])
fit.pareto.truncT.1[,5] = log(fit.pareto.truncT.1[,5])
fit.pareto.truncT.3[,5] = log(fit.pareto.truncT.3[,5])

names(fit.pareto.truncT.1) <- names(fit.pareto.truncT.3) <- c("case","method","time","val","log(hat(lambda))","hat(vartheta)")#,"deg")

par.logskew <- as.data.frame(par.logskew);par.truncT <- as.data.frame(par.truncT)
names(par.logskew) = names(fit.pareto.logskew.1)[5:9];names(par.truncT) = names(fit.pareto.truncT.1)[5:7]
par.logskew$case = 1:nrow(par.logskew);par.truncT$case = 1:nrow(par.truncT)

levels = c("hat(lambda)","hat(vartheta)","hat(b)[1]","hat(b)[2]")
data_true = par.logskew

data_true <- pivot_longer(data_true, cols=levels, names_to = "Variable", values_to = "Value")
data_true$facet = factor(paste0(data_true$Variable),levels=levels)
data_true <- rbind(data_true,data_true)
data_true$method = rep(1:2,each=nrow(data_true)/2)

custom_scales <- list(
  `hat(lambda)` =  c(4, 14),
  `hat(vartheta)` = c(0.8, 1.7),
  `hat(b)[1]` =  c(-8, 3),
  `hat(b)[2]` = c(-8, 3)
)

p.list.logskew <- list()
thres.b = Inf
data = fit.pareto.logskew.1
data.true.sub <- data_true 
data_long <- pivot_longer(data, cols=levels, names_to = "Variable", values_to = "Value")
data_long$facet = factor(paste0(data_long$Variable),level=levels)

for(i in 1:4){
    p <- ggplot(subset(data_long,facet==levels[i]), aes(x = factor(case), y = Value,fill=factor(method,labels=c("Score","Spectral")))) +
    geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
    geom_point(data=subset(data.true.sub,facet==levels[i]),aes(x=factor(case),y=Value),color="black",size=1,position=position_dodge(width = 1)) +
    labs(x = "Cases",y = parse(text=levels[i]),fill="Method") + 
    theme(axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        strip.text = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title =element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16),legend.position = "none") + coord_cartesian(ylim = custom_scales[[i]])
    p.list.logskew[[i]] <- p
}

data = fit.pareto.logskew.3
data.true.sub <- data_true 
data_long <- pivot_longer(data, cols=levels, names_to = "Variable", values_to = "Value")
data_long$facet = factor(paste0(data_long$Variable),level=levels)
for(i in 1:4){
    p <- ggplot(subset(data_long,facet==levels[i]), aes(x = factor(case), y = Value,fill=factor(method,labels=c("Score","Spectral")))) +
    geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
    geom_point(data=subset(data.true.sub,facet==levels[i]),aes(x=factor(case),y=Value),color="black",size=1,position=position_dodge(width = 1)) +
    labs(x = "Cases",y = parse(text=levels[i]),fill="Method") +
    theme(axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        strip.text = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16),legend.position = "none") + coord_cartesian(ylim = custom_scales[[i]])
    p.list.logskew[[i+4]] <- p
}

# r-Pareto truncT
levels = c("log(hat(lambda))","hat(vartheta)")
data_true = par.truncT[,-3]
data_true[,1] = log(data_true[,1])
data_true <- pivot_longer(data_true, cols=levels, names_to = "Variable", values_to = "Value")
data_true$facet = factor(paste0(data_true$Variable),levels=levels)
data_true <- rbind(data_true,data_true)
data_true$method = rep(1:2,each=nrow(data_true)/2)

p.list.truncT <- list()
thres.lambda = Inf
data = fit.pareto.truncT.1
data.true.sub <- data_true
data_long <- pivot_longer(data, cols=levels, names_to = "Variable", values_to = "Value")
data_long$facet = factor(paste0(data_long$Variable),level=levels)
for(i in 1:2){
    p <- ggplot(subset(data_long,facet==levels[i]), aes(x = factor(case), y = Value,fill=factor(method,labels=c("Score","Spectral")))) +
    geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
    geom_point(data=subset(data.true.sub,facet==levels[i]),aes(x=factor(case),y=Value),color="black",size=1,position=position_dodge(width = 1)) +
    labs(x = "Cases",y = parse(text=levels[i]),fill="Method") +
    theme(axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        strip.text = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16),legend.position = "none") 
    p.list.truncT[[i]] <- p
}


data = fit.pareto.truncT.3
data.true.sub <- data_true
data_long <- pivot_longer(data, cols=levels, names_to = "Variable", values_to = "Value")
data_long$facet = factor(paste0(data_long$Variable),level=levels)
for(i in 1:2){
    p <- ggplot(subset(data_long,facet==levels[i]), aes(x = factor(case), y = Value,fill=factor(method,labels=c("Score","Spectral")))) +
    geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
    geom_point(data=subset(data.true.sub,facet==levels[i]),aes(x=factor(case),y=Value),color="black",size=1,position=position_dodge(width = 1)) +
    labs(x = "Cases",y = parse(text=levels[i]),fill="Method") + 
    theme(axis.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        strip.text = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16),legend.position = "none")
    p.list.truncT[[i+2]] <- p
}

pdf(file=paste0("figures/simulation_pareto_violin_logskew",idx.file,".pdf"),width=14,height = 5,onefile = TRUE)
for(i in 1:length(p.list.logskew)) print(p.list.logskew[[i]])
dev.off()

pdf(file=paste0("figures/simulation_pareto_violin_truncT",idx.file,".pdf"),width=14,height = 5,onefile = TRUE)
for(i in 1:length(p.list.truncT)) print(p.list.truncT[[i]])
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

## max-stable skewed BR ##
files.logskew <- list.files(path="data/simulation_logskew/",pattern="simulation_comp2_2_logskew_\\d+_500_1.RData",full.names=TRUE,recursive=FALSE)

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
data_true <- subset(data_true,case==1)
p.list.logskew <- list()

data = subset(fit.logskew,case==1)
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


## collect results: applications ##
# Red Sea #
data.RedSea <- data.frame(
    lambda = numeric(),
    vartheta = numeric(),
    theta = numeric(),
    a = numeric(),
    b1 = numeric(),
    b2 = numeric(),
    b3 = numeric(),
    b4 = numeric(),
    b5 = numeric(),
    case = character(),
    method = character(),
    time = numeric(),
    type = character(),
    lik = numeric(),
    stringsAsFactors = FALSE
)
count = 1
cases = c("BR","sBR","BR","BR")
method = c("spectral","spectral","composite","vecchia")
methods = c("Nelder-Mead","L-BFGS-B","Nelder-Mead","Nelder-Mead")
for(i.case in 1:4){
    file.origin = paste0("data/application2/application_RedSea_results_",i.case,"_",methods[i.case],".RData")
    load(file.origin,e <- new.env())
    if(i.case %in% c(1,2)){ 
        data.RedSea[count,1:9] = e$fit.result$par;data.RedSea$time[count] = e$fit.result$time[[3]]
    }else{
        data.RedSea[count,1:9] = c(e$fit.result$par[1:4],rep(0,5)); data.RedSea$time[count] = e$fit.result$time 
    }
    data.RedSea$case[count] = cases[i.case]
    data.RedSea$method[count] = method[i.case]
    data.RedSea$type[count] = "full"
    data.RedSea$lik[count] = e$fit.result$value
    count = count + 1 
    for(j in 1:300){
        file = paste0("data/application2/application_RedSea_results_",i.case,"_",methods[i.case],"_",j,".RData")
        if(file.exists(file)){
            load(file,e <- new.env())
            if(i.case %in% c(1,2)){ 
                data.RedSea[count,1:9] = e$fit.result$par;data.RedSea$time[count] = e$fit.result$time[[3]]
            }else{
                data.RedSea[count,1:9] = c(e$fit.result$par[1:4],rep(0,5)); data.RedSea$time[count] = e$fit.result$time 
            }
            data.RedSea$case[count] = cases[i.case]
            data.RedSea$method[count] = method[i.case]
            data.RedSea$type[count] = "partial"
            data.RedSea$lik[count] = e$fit.result$value
            count = count + 1
        }
    }
}

# Florida Tampa Bay #
data.florida <- data.frame(
    lambda = numeric(),
    vartheta = numeric(),
    theta = numeric(),
    a = numeric(),
    b1 = numeric(),
    b2 = numeric(),
    b3 = numeric(),
    b4 = numeric(),
    case = character(),
    method = character(),
    time = numeric(),
    type = character(),
    lik = numeric(),
    stringsAsFactors = FALSE
)
cases = c("BR-SUM","BR-MAX","sBR-SUM","sBR-MAX","BR-SUM","BR-MAX")
method = c("spectral","spectral","spectral","spectral","scoreMatching","scoreMatching")
count = 1
methods = c("Nelder-Mead","Nelder-Mead","L-BFGS-B","L-BFGS-B","Nelder-Mead","Nelder-Mead")
for(i.case in 1:6){
    file.origin = paste0("data/application2/application_florida_results_",i.case,"_",methods[i.case],".RData")
    load(file.origin,e <- new.env())
    data.florida[count,1:8] = e$fit.result$par;data.florida$time[count] = e$fit.result$time[[3]]
    data.florida$case[count] = cases[i.case]
    data.florida$method[count] = method[i.case]
    data.florida$type[count] = "full"
    data.florida$lik[count] = e$fit.result$value
    count = count + 1 
    for(j in 1:70){
        file = paste0("data/application2/application_florida_results_",i.case,"_",methods[i.case],"_",j,".RData")
        if(file.exists(file)){
            load(file,e <- new.env())
            data.florida[count,1:8] = e$fit.result$par;data.florida$time[count] = e$fit.result$time[[3]]
            data.florida$case[count] = cases[i.case]
            data.florida$method[count] = method[i.case]
            data.florida$type[count] = "partial"
            data.florida$lik[count] = e$fit.result$value
            count = count + 1
        }
    }
}
data.florida$lambda = data.florida$lambda*2/0.03128403
data.florida$theta = data.florida$theta
data.RedSea$theta = data.RedSea$theta
save(data.florida,data.RedSea,file="data/application_results.RData")

## analysis ## 
load("data/application_results.RData")

# Red Sea #
data= subset(data.RedSea,type=="full")
median(subset(data.RedSea,type=="partial" & method=="spectral" & case=="BR")$time)
median(subset(data.RedSea,type=="partial" & method=="spectral" & case=="sBR")$time)
median(subset(data.RedSea,type=="partial" & method=="composite" & case=="BR")$time)
median(subset(data.RedSea,type=="partial" & method=="vecchia" & case=="BR")$time)

# CI for Red Sea #
subset(data.RedSea,type=="full")
round(apply(subset(data.RedSea,type=="partial" & method=="spectral" & case=="BR")[,1:9],2,sd),2)
round(apply(subset(data.RedSea,type=="partial" & method=="spectral" & case=="sBR")[,1:9],2,sd),2)
round(apply(subset(data.RedSea,type=="partial" & method=="composite" & case=="BR")[,1:9],2,sd),2)
round(apply(subset(data.RedSea,type=="partial" & method=="vecchia" & case=="BR")[,1:9],2,sd),2)

RedSea_latex = cbind(paste(data$case,data$method),round(data[,1:9],2),data$lik*2+c(4,9,4,4)*2,round(data$time,2))
xtable::xtable(RedSea_latex)

# Florida Tampa Bay #
data= subset(data.florida,type=="full")
median(subset(data.florida,type=="partial" & method=="scoreMatching" & case == "BR-SUM")$time)
median(subset(data.florida,type=="partial" & method=="scoreMatching" & case == "BR-MAX")$time)
median(subset(data.florida,type=="partial" & method=="spectral" & case == "BR-SUM")$time)
median(subset(data.florida,type=="partial" & method=="spectral" & case == "BR-MAX")$time)
median(subset(data.florida,type=="partial" & method=="spectral" & case == "sBR-SUM")$time)
median(subset(data.florida,type=="partial" & method=="spectral" & case == "sBR-MAX")$time)
Florida_latex = cbind(paste(data$case,data$method),round(data[,1:8],2),data$lik*2+2*c(4,4,8,8,4,4),round(data$time,2))
xtable::xtable(Florida_latex)
# CI for Florida Tampa Bay # score matching
subset(data.florida,type=="full")
jack.knife(as.matrix(subset(data.florida,type=="partial" & method=="scoreMatching" & case == "BR-SUM")[,1:8]),unlist(subset(data.florida,type=="full" & method=="scoreMatching" & case == "BR-SUM")[1:8]))
jack.knife(as.matrix(subset(data.florida,type=="partial" & method=="scoreMatching" & case == "BR-MAX")[,1:8]),unlist(subset(data.florida,type=="full" & method=="scoreMatching" & case == "BR-MAX")[1:8]))
jack.knife(as.matrix(subset(data.florida,type=="partial" & method=="spectral" & case == "BR-SUM")[,1:8]),unlist(subset(data.florida,type=="full" & method=="spectral" & case == "BR-SUM")[1:8]))
jack.knife(as.matrix(subset(data.florida,type=="partial" & method=="spectral" & case == "BR-MAX")[,1:8]),unlist(subset(data.florida,type=="full" & method=="spectral" & case == "BR-MAX")[1:8]))
jack.knife(as.matrix(subset(data.florida,type=="partial" & method=="spectral" & case == "sBR-SUM")[,1:8]),unlist(subset(data.florida,type=="full" & method=="spectral" & case == "sBR-SUM")[1:8]))
jack.knife(as.matrix(subset(data.florida,type=="partial" & method=="spectral" & case == "sBR-MAX")[,1:8]),unlist(subset(data.florida,type=="full" & method=="spectral" & case == "sBR-MAX")[1:8]))



jack.knife <- function(jack,knife){
    n = nrow(jack)
    jack = n*jack - (n-1) * matrix(knife,nrow=n,ncol=ncol(jack),byrow=TRUE)
    est.sd = apply(jack,2,sd)/sqrt(n)
    return(est.sd)
}

