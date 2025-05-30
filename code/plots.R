rm(list=ls())
library(ggplot2)
library(parallel)
library(gridExtra)
library(directlabels)
library(Rfast)
library(RColorBrewer)
library(mev)
library(ggpubr)
library(tidyr)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
set.seed(1)
para.alpha = rbind(c(1,0,0),c(1,-1,-2),c(1,-1,1)) ## slant parameter for skewed norm model ##
d = 32
coord = as.matrix(expand.grid(1:d,1:d))
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))

basis.list <- list()
centers <- rbind(c(0.25,0.25),c(0.25,0.25),c(0.75,0.75))*d
idx.centers <- apply(centers,1,function(x){which.min(apply(coord,1,function(y){sum((x-y)^2)}))})
basis <- sapply(idx.centers,function(x){y=dnorm(diff.mat[x,],mean=0,sd=d*2);y=y-mean(y);y/sqrt(sum(y^2))})
basis[,1] = rep(0,d^2);basis[1:floor(d^2/2),1] = 0.1; basis[(d^2-floor(d^2/2)+1):d^2,1] = -0.1
basis.list[[1]] <- basis    

idx = floor(matrix(seq(1,nrow(coord),length.out=6),ncol=2,3))
basis <- sapply(1:(ncol(para.alpha)),function(x){y <- rep(0,nrow(coord));y[idx[x,]] <- c(-2,2);y})
basis[,1] = rep(0,d^2);basis[1:floor(d^2/2),1] = 0.1; basis[(d^2-floor(d^2/2)+1):d^2,1] = -0.1
basis.list[[2]] <- basis

all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))

#sigma = cov.func(diff.mat,c(8,1.5))
sigma = vario.func(coord,c(4,1))
par.list.1 = apply(para.alpha,1,function(x){alpha2delta(list(sigma,alpha.func(par=x,b.mat=basis.list[[1]])))})
par.list.2 = apply(para.alpha,1,function(x){alpha2delta(list(sigma,alpha.func(par=x,b.mat=basis.list[[2]])))})

### plot the deltas ###
p.list1 <- p.list2 <- list()
for(i in 1:length(par.list.1)){
    data = data.frame(x=coord[,1],y=coord[,2],z=par.list.1[[i]][[2]])
    print(range(data$z))
    p.list1[[i]] <- ggplot(data, aes(x = x, y = y, fill = z)) +
    geom_tile() + ggtitle(paste(c("alpha:",para.alpha[i,]),collapse = " ")) + 
    scale_fill_gradient2(low = "blue",mid="white" ,high = "red",limits=c(min(-0.5,min(data$z)),max(0.5,max(data$z)))) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + coord_fixed()
    data = data.frame(x=coord[,1],y=coord[,2],z=par.list.2[[i]][[2]])
    print(range(data$z))
    p.list2[[i]] <- ggplot(data, aes(x = x, y = y, fill = z)) +
    geom_tile() + ggtitle(paste(c("alpha:",para.alpha[i,]),collapse = " ")) + 
    scale_fill_gradient2(low = "blue",mid="white" ,high = "red") +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + coord_fixed()
}

pdf("figures/deltas_final.pdf",width=4.5*3,height=4.5,onefile=TRUE)
grid.arrange(grobs=c(p.list1),ncol=length(par.list.1))
grid.arrange(grobs=c(p.list2),ncol=length(par.list.2))
dev.off()

## plot the extremal coefficients map##
data = NULL
coord.center = matrix(NA,2,2)
for(idx in 1:2){
    idx.center = rbind(c(11,11),c(25,16))[idx,]
    idx.center = which.min(abs(coord[,1] - idx.center[1]) + abs(coord[,2] - idx.center[2]))
    ind.idx.center = all.pairs[1,] == idx.center |  all.pairs[2,] == idx.center
    idx2 = apply(all.pairs[,ind.idx.center],2,function(x){x[x!=idx.center]})
    coord.center[idx,] = coord[idx.center,]
    for(idx.case in 1:length(par.list.1)){
        true.ext.coef <- unlist(mclapply(all.pairs.list[ind.idx.center],true_extcoef,par=par.list.1[[idx.case]],model="logskew1",mc.cores=5,mc.set.seed = FALSE))
        data <- rbind(data,data.frame( x = coord[idx2,1],
                        y = coord[idx2,2],
                        z = true.ext.coef,idx.case=paste0("b[1]==",para.alpha[idx.case,2],"~b[2]==",para.alpha[idx.case,3]),idx.center=idx))
        true.ext.coef <- unlist(mclapply(all.pairs.list[ind.idx.center],true_extcoef,par=par.list.2[[idx.case]],model="logskew1",mc.cores=5,mc.set.seed = FALSE))

        if(idx.case == 1 ){
            true.ext.coef <- unlist(mclapply(all.pairs.list[ind.idx.center],true_extcoef,par=list(par.list.1[[1]][[1]],rep(0,nrow(coord))),model="logskew1",mc.cores=5,mc.set.seed = FALSE))
            data <- rbind(data,data.frame( x = coord[idx2,1],
                            y = coord[idx2,2],
                            z = true.ext.coef,idx.case="Brown-Resnick",idx.center=idx))
        }
    }
}

loc.df = data.frame(x = coord[idx.centers,][-1,1], y = coord[idx.centers,][-1,2],z=NA)
brks = round(quantile(data$z,probs=seq(0.001,0.999,length.out=5)),4)
data1 = subset(data,idx.center==1)
labels = unique(data$idx.case)[c(2,1,3,4)]
data1$idx.case = factor(data1$idx.case,levels=labels,labels=labels)
levels(data1$idx.case)
p1 <- ggplot(data1, aes(x = x, y = y, z = z))  + 
        geom_tile(aes(fill = z)) + 
        facet_wrap(~ idx.case, labeller = label_parsed) +
        scale_fill_distiller(palette = "RdBu") +
        geom_contour(colour = "black", breaks = brks) + 
        geom_dl(aes(label = sprintf("%.2f",..level..)), method = "bottom.pieces", breaks = brks, stat = "contour") + 
        geom_point(data=loc.df,aes(x=x,y=y),size=2,fill="black") +
        theme(axis.text = element_text(size = 12), 
              strip.text = element_text(size = 16),
              axis.title.x = element_text(size = 16), 
              axis.title.y = element_text(size = 16), 
              plot.title = element_text(hjust = 0.5, size = 16),
              legend.title = element_text(size = 16)) + 
        coord_fixed() + 
        labs(x = expression(s[1]), y = expression(s[2]), fill = expression(theta[2])) +
        # Add a star to the center of the grid
        geom_text(aes(x = coord.center[1,1], y = coord.center[1,2], label = "*"), size = 10, color = "red", vjust = 0.8, hjust = 0.5)
p1

data2 = subset(data,idx.center==2)
data2$idx.case = factor(data2$idx.case,levels=labels,labels=labels)
p2 <- ggplot(data2, aes(x = x, y = y,z=z))  + 
        geom_tile(aes(fill=z)) + facet_wrap(~ idx.case,labeller=label_parsed) +
        scale_fill_distiller(palette="RdBu") +
        geom_contour(colour="black",breaks=brks) + 
        geom_dl(aes(label=sprintf("%.2f",..level..)),method="bottom.pieces",breaks=brks, 
                stat="contour") + 
        geom_point(data=loc.df,aes(x=x,y=y),size=2,fill="black") +
        theme(axis.text = element_text(size=12), 
                            strip.text = element_text(size = 16),
                            axis.title.x = element_text(size=16), 
                            axis.title.y = element_text(size=16), 
                            plot.title = element_text(hjust = 0.5,size=16),legend.title = element_text(size=16)) + coord_fixed() + 
        labs(x = expression(s[1]), y = expression(s[2]),fill=expression(theta[2])) + 
        geom_text(aes(x = coord.center[2,1], y = coord.center[2,2], label = "*"), size = 10, color = "red", vjust = 0.8, hjust = 0.5)
p2

pdf("figures/extcoef_final_logskew.pdf",width=4*4-1,height=4*2-1,onefile=TRUE)
grid.arrange(grobs=list(p1,p2),ncol=2)
dev.off()

p.list = list()
sigma.22 = 10
rho = 1:9/10*sigma.22
BR.values = unlist(lapply(rho,function(x){V_bi_logskew(c(1,1),list(matrix(c(sigma.22,x,x,sigma.22),2,2),c(0,0)))}))
for(i in 1:length(rho)){
    r = sqrt(min(eigen(matrix(c(1,rho[i]/sigma.22,rho[i]/sigma.22,1),2))$values))
    delta = seq(-r,r,length.out=100)
    delta.grid = as.matrix(expand.grid(delta,delta))
    delta.grid.list <- split(delta.grid,row(delta.grid))

    idx.valid = apply(delta.grid,1,function(x){sum(x^2) < r^2}) # & abs(diff(x))<sqrt(2-2*rho)
    
    values <- unlist(lapply(delta.grid.list[idx.valid],function(x){V_bi_logskew(c(1,1),list(sigma=matrix(c(sigma.22,rho[i],rho[i],sigma.22),2,2),delta=x))}))

    data = data.frame(x=delta.grid[idx.valid,1],y=delta.grid[idx.valid,2],z=values)
    
    data2 = data.frame(x=0,y=0,z=BR.values[i])
    p.list[[i]] <- ggplot(data) + geom_point(aes(x=x,y=y,color=z)) + scale_color_distiller(palette="RdBu") + ggtitle(paste("rho",rho[i])) + coord_fixed() + theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + geom_point(data=data2,aes(x=x, y = y), color = "black") + ggtitle(paste("Max:",round(max(data$z),4),"BR:",round(BR.values[i],4)))
}

pdf("figures/bivariate_extcoef_rho.pdf",width=5*3,height = 5*3,onefile = TRUE)
grid.arrange(grobs=p.list,ncol=3,nrow=3)
dev.off()


p.list = list()
sigma.22.vec = c(1,5,10)
for(i in 1:length(sigma.22.vec)){
    sigma.22 = sigma.22.vec[i]
    rho = seq(0,0.999,length.out=100)*sigma.22
    BR.values = unlist(lapply(rho,function(x){V_bi_logskew(c(1,1),list(matrix(c(sigma.22,x,x,sigma.22),2,2),c(0,0)))}))
    alpha = seq(0,5,length.out=100)
    para.grid = as.matrix(expand.grid(alpha,rho))
    values <- lapply(split(para.grid,row(para.grid)),function(x){par.list = alpha2delta(list(matrix(c(sigma.22,x[2],x[2],sigma.22),2,2),alpha=c(x[1],-x[1])));V_bi_logskew(c(1,1),par.list)})
    data = data.frame(x=para.grid[,1],y=para.grid[,2],z=unlist(values))
    data2 = data.frame(x=0,y=para.grid[,2],z=BR.values)

    brks = round(quantile(data$z,probs=seq(0.01,0.99,length.out=4)),2)

    p.list[[i]] <- ggplot(data, aes(x = x, y = y, z=z))  + 
        geom_tile(aes(fill=z)) +
        scale_fill_distiller(palette="RdBu",limits=c(1,2)) +
        geom_contour(colour="black",breaks=brks) + 
        geom_dl(aes(label=sprintf("%.2f",..level..)), breaks=brks,stat="contour",method=list(dl.trans(x = max(0.01,x),y=max(0.01,y)), "smart.grid")) + 
        labs(title=bquote(~sigma[22]==.(sigma.22)),x=expression(eta[1]),y=expression(sigma[12]),fill=expression(theta[2])) +
        theme(axis.text = element_text(size=10), 
                        strip.text = element_text(size = 14),
                        axis.title.x = element_text(size=14), 
                        axis.title.y = element_text(size=14), 
                        plot.title = element_text(hjust = 0.5,size=14),legend.title = element_text(size=14),legend.position = "none")            
}

p.list[[length(sigma.22.vec)+1]] <- as_ggplot(get_legend(p.list[[1]]+theme(legend.position = "right")))

pdf("figures/bivariate_extcoef_rho_alpha.pdf",width=10,height = 3,onefile = TRUE)
grid.arrange(grobs=p.list,nrow=1,widths=c(rep(5,length(sigma.22.vec)),1))
dev.off()

#val.mat = matrix(unlist(values),ncol=length(alpha),byrow=TRUE)

## plot the true extremal coef for the truncated extremal t model ##
sigma.22 = 1
df = c(2,5,10)
rho = seq(0,sigma.22-0.001,length.out=200)
data = NULL
for(i in 1:length(df)){
    par.truncT.list = lapply(rho,function(x){list(matrix(c(sigma.22,x,x,sigma.22),2,2),df[i])})
    true.ext.truncT <- unlist(lapply(par.truncT.list,V_truncT,x=c(1,1)))
    true.ext.t <- unlist(lapply(1:length(par.truncT.list),function(id) mev::expme(z=rep(1,2),par=list(Sigma=par.truncT.list[[id]][[1]],df=par.truncT.list[[id]][[2]]),model="xstud") ))
    data = rbind(data,data.frame(x=rho,truncT=true.ext.truncT,extT = true.ext.t,df=df[i]))
}
new_labels <- paste0("nu==",df)
data$df <- factor(data$df,labels=new_labels)
p <- ggplot(data) + geom_line(aes(x=x,y=truncT,color="Truncated extremal-t"),linewidth=1) + 
geom_line(aes(x=x,y=extT,color="Extremal-t"),linewidth=1) + coord_cartesian(ylim = c(1, 2)) + geom_hline(yintercept=c(1,2),linetype="dashed") +
facet_wrap(~df,labeller =label_parsed) + theme(legend.position = "none") + labs(x=expression(rho),y=expression(theta[2]),color=" ") +
theme(axis.text = element_text(size=10), 
            strip.text = element_text(size = 14), 
            axis.title.x = element_text(size=14), 
            axis.title.y = element_text(size=14), 
            plot.title = element_text(hjust = 0.5,size=14),legend.title = element_text(size=14),legend.text = element_text(size=14))

pdf("figures/bivariate_extcoef_truncT.pdf",width=10,height = 3,onefile = TRUE)
p
dev.off()


## explore the likelihood ## 
## setting 1: covariance function
d = 10
ncores=10
coord = as.matrix(expand.grid(1:d,1:d))
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
nu.list = c(1,2,5,10)
para.range.list = c(4,8)
idx.para=1:3
basis = matrix(0,nrow=nrow(coord),ncol=3)
result2 <- result1 <- create_lists(c(length(nu.list),length(para.range.list)))
for(i in 1:length(nu.list)){
    for(j in 1:length(para.range.list)){
        nu = nu.list[i]
        para.range = para.range.list[j]
        init = cbind(seq(1,20,length.out=100),nu,1,0,0,0)
        data <- simu_logskew(m=10000,par=alpha2delta(list(cov.func(diff.mat,c(para.range,nu,1)),alpha.func(c(0,0,0),basis))),ncores=ncores)
        result2[[i]][[j]] <- apply(init[,idx.para], 1, function(x) fit.model(data=data,init=x,fixed=c(F,F,F),loc=diff.mat,FUN=cov.func,thres=100,model="BR",lb=rep(-Inf,3),ub=rep(Inf,3),ncores=ncores,maxit=1000,trace=FALSE,method="Nelder-Mead",opt=FALSE,hessian=FALSE,idx.para=idx.para)$value)
        result1[[i]][[j]] <- apply(init, 1, function(x) fit.model(data=data,init=x,fixed=c(F,F,F,F,F,F),loc=diff.mat,FUN=cov.func,alpha.func=alpha.func,basis=basis,thres=100,model="logskew",lb=rep(-Inf,6),ub=rep(Inf,6),ncores=ncores,maxit=1000,trace=FALSE,method="Nelder-Mead",opt=FALSE,hessian=FALSE,idx.para=idx.para)$value)
    }
}

save(result2,result1,nu.list,para.range.list,idx.para,diff.mat,file="data/log_likelihood_BR_cov.RData")

load("data/log_likelihood_BR_cov.RData")
pdf("figures/log_likelihood_BR_cov.pdf",width=5*4,height=4*2)
par(mfrow=c(length(para.range.list),length(nu.list)),mgp=c(2.5,1,0),cex.main=2,cex.lab = 2,cex.axis=2,mar=c(5,5,4,1))
range.seq = seq(1,20,length.out=100)
for(j in 1:length(para.range.list)){
    for(i in 1:length(nu.list)){
        plot(range.seq,result2[[i]][[j]],xlab=bquote(lambda),ylab="Log likelihood",main=bquote(sigma==.(nu.list[i])~~lambda==.(para.range.list[j])~~nu==1),type="l",lwd=2)
        lines(range.seq,result1[[i]][[j]],col="red",lwd=2)
        idx.min = which.min(result2[[i]][[j]])
        abline(v=range.seq[idx.min],col="grey",lwd=2)
        abline(v=para.range.list[j],col="blue",lwd=2)
    }
}
dev.off()

## setting 2: variogram function
para.range.list = c(2,4,8,10)
para.smooth.list = c(0.5,1,1.5)
idx.para = 1:2
result2 <- result1 <- create_lists(c(length(para.range.list),length(para.smooth.list)))
for(i in 1:length(para.range.list)){
    for(j in 1:length(para.smooth.list)){
        para.range = para.range.list[i]
        para.smooth = para.smooth.list[j]
        data <- simu_logskew(m=10000,par=alpha2delta(list(vario.func(coord,c(para.range,para.smooth)),alpha.func(c(0,0,0),basis))),ncores=ncores)
        init = cbind(seq(1,20,length.out=100),para.smooth,0,0,0)
        result2[[i]][[j]] <- apply(init[,idx.para], 1, function(x) fit.model(data=data,init=x,fixed=c(F,F),loc=coord,FUN=vario.func,thres=100,model="BR",lb=rep(-Inf,2),ub=rep(Inf,2),ncores=ncores,maxit=1000,trace=FALSE,method="Nelder-Mead",opt=FALSE,hessian=FALSE,idx.para=idx.para)$value)
        result1[[i]][[j]] <- apply(init, 1, function(x) fit.model(data=data,init=x,fixed=c(F,F,F,F,F),loc=coord,FUN=vario.func,alpha.func=alpha.func,basis=basis,thres=100,model="logskew",lb=rep(-Inf,5),ub=rep(Inf,5),ncores=ncores,maxit=1000,trace=FALSE,method="Nelder-Mead",opt=FALSE,hessian=FALSE,idx.para=idx.para)$value)
    }
}

save(result2,result1,para.range.list,para.smooth.list,idx.para,coord,file="data/log_likelihood_BR_vario.RData")

load("data/log_likelihood_BR_vario.RData")
pdf("figures/log_likelihood_BR_vario.pdf",width=5*4,height=4*3)
par(mfrow=c(length(para.smooth.list),length(para.range.list)),mgp=c(2.5,1,0),cex.main=2,cex.lab = 2,cex.axis=2,mar=c(5,5,4,1))
range.seq = seq(1,20,length.out=100)
for(j in 1:length(para.smooth.list)){
    for(i in 1:length(para.range.list)){
        plot(range.seq,result2[[i]][[j]],xlab=bquote(lambda),ylab="Log likelihood",main=bquote(lambda==.(para.range.list[i])~~nu==.(para.smooth.list[j])),type="l",lwd=2)
        lines(range.seq,result1[[i]][[j]],col="red",lwd=2)
        idx.min = which.min(result2[[i]][[j]])
        abline(v=range.seq[idx.min],col="grey",lwd=2)
        abline(v=para.range.list[i],col="blue",lwd=2)
    }
}
dev.off()

## plot the boxplot for the simulation study ## 
load("data/simulation_study_logskew_results_vario_100_new.RData",e1<-new.env())
data = cbind(1,rep(1:nrow(e1$par.skew.normal),times=sapply(e1$est.mat.list[[1]],nrow)),1,do.call(rbind,e1$est.mat.list[[1]]))

data = as.data.frame(data);names(data) = c("thres","id","type","hat(lambda)","hat(vartheta)","hat(b)[1]","hat(b)[2]")
str(data)

# Reshape the data to a long format
data_long <- pivot_longer(data, -c(thres, id, type), names_to = "variable", values_to = "value")
#data_long$facet = paste0("Spline~Type==", data_long$type,"~~",data_long$variable)
data_long$facet = paste0(data_long$variable)
data_long$facet = factor(data_long$facet)
par.skew.normal = as.data.frame(e1$par.skew.normal)
colnames(par.skew.normal) = c("hat(lambda)","hat(vartheta)","hat(b)[1]","hat(b)[2]")
par.skew.normal$id = 1:nrow(par.skew.normal)

par.skew.normal$type = 1
par.skew.normal$thres = 1



par.skew.normal_long <- pivot_longer(par.skew.normal, -c(id,type,thres), names_to = "variable", values_to = "value")
levels=levels(data_long$facet)[c(3,4,1,2,7,8,5,6)]

par.skew.normal_long$facet = factor(paste0(par.skew.normal_long$variable),levels=levels)

p <- ggplot(data_long, aes(x = factor(id), y = value)) +
  geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
  geom_point(data=par.skew.normal_long,aes(x=factor(id),y=value),color="black",size=1,position=position_dodge(width = 1)) +
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
pdf("figures/simulation_est_boxplots_final.pdf",width=10,height=6,onefile=TRUE)
p
dev.off()


load("data/simulation_study_logskew_results_vario_30_BR.RData",e1<-new.env())
load("data/simulation_study_logskew_results_vario_30_BR_comp.RData",e2<-new.env())
data1 = cbind(rep(1:nrow(e1$par.skew.normal),times=sapply(e1$est.mat.list,nrow)),1,do.call(rbind,e1$est.mat.list))
data1 = as.matrix(data1)
data2 = cbind(rep(1:nrow(e2$par.skew.normal),times=sapply(e2$est.mat.list[[1]],nrow)),2,do.call(rbind,e2$est.mat.list[[1]]))
data = rbind(data1,data2)
data = as.data.frame(data);names(data) = c("id","type","hat(lambda)","hat(vartheta)","hat(b)[1]","hat(b)[2]")
data = data[,1:4]
str(data)

# Reshape the data to a long format
data_long <- pivot_longer(data, -c(id, type), names_to = "variable", values_to = "value")
data_long$facet = factor(data_long$variable)
par.skew.normal = as.data.frame(e1$par.skew.normal[,1:2])
colnames(par.skew.normal) = c("hat(lambda)","hat(vartheta)")
par.skew.normal$id = 1:nrow(par.skew.normal)
par.skew.normal = rbind(par.skew.normal,par.skew.normal)
par.skew.normal$type = rep(1:2,each=nrow(par.skew.normal)/2)
par.skew.normal_long <- pivot_longer(par.skew.normal, -c(id,type), names_to = "variable", values_to = "value")
par.skew.normal_long$facet = factor(par.skew.normal_long$variable)

rel_bias_1 <- rel_bias_2 <- matrix(NA,nrow=6,ncol=2)
for(idx in 1:6){
    data_sub_1 = subset(data,id==idx & type==1)
    data_sub_2 = subset(data,id==idx & type==2)
    rel_bias_1[idx,] <- c(mean(data_sub_1[,3]-par.skew.normal[idx,1])/par.skew.normal[idx,1],mean(data_sub_1[,4]-par.skew.normal[idx,2])/par.skew.normal[idx,2])
    rel_bias_2[idx,] <- c(mean(data_sub_2[,3]-par.skew.normal[idx,1])/par.skew.normal[idx,1],mean(data_sub_2[,4]-par.skew.normal[idx,2])/par.skew.normal[idx,2])
}

p <- ggplot(data_long, aes(x = factor(id), y = value,fill=factor(type,labels=c("Angular","Composite")))) +
  geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
  #geom_boxplot(position = position_dodge(width=0.75),width=0.4, color="grey", alpha=1) +
  geom_point(data=par.skew.normal_long,aes(x=factor(id),y=value),color="black",size=1,
  position=position_dodge(width = 1)) +
  facet_wrap(~ facet, scales = "free",ncol=2,nrow=1,labeller = label_parsed) +
  labs(x = "Cases",
       y = "Value")  + 
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16),legend.position = "none")

pdf("figures/simulation_est_boxplots_final_BR_30.pdf",width=5*2,height=4,onefile=TRUE)
p
dev.off()

## plot the boxplot for the simulation study:composite for skewedBR ##
load("data/simulation_comp_results.RData",e<-new.env())

colnames(e$data.est)[1:5] <- c("hat(lambda)","hat(vartheta)","hat(b)[0]","hat(b)[1]","hat(b)[2]")
data_long <- pivot_longer(e$data.est[,-3], -c(method,step,case), names_to = "variable", values_to = "value")
data_long$variable = factor(data_long$variable,levels=names(e$data.est)[c(1,2,4,5)])
names(e$par.skew.normal) <- names(e$data.est)[c(1,2,4,5)]

e$par.skew.normal = do.call(rbind,replicate(4,e$par.skew.normal,simplify = FALSE))
tmp = expand.grid(method=unique(e$data.est$method),step=unique(e$data.est$step),case=unique(e$data.est$case))
e$par.skew.normal = cbind(e$par.skew.normal[tmp$case,],tmp)
colnames(e$par.skew.normal) <- colnames(e$data.est[,-3])

par.skew.normal_long <- pivot_longer(e$par.skew.normal, -c(method,step,case), names_to = "variable", values_to = "value")
par.skew.normal_long$variable <- factor(par.skew.normal_long$variable,levels=names(e$data.est)[c(1,2,4,5)])
p <- ggplot(subset(data_long,case==1 & step==2), aes(x = factor(step), y = value,fill=factor(method,labels=c("Angular","Composite")))) +
  geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
  geom_point(data=subset(par.skew.normal_long,case==1 & step==2),aes(x=factor(step),y=value),color="black",size=2, position=position_dodge(width = 1)) +
  facet_wrap(~variable, scales = "free",ncol=4,nrow=1,labeller = label_parsed) +
  labs(x = "Step",
       y = "Value")  + 
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16),legend.position = "none")
p

pdf("figures/simulation_est_boxplots_final_skewBR_comp.pdf",width=3*4,height=4,onefile=TRUE)
show(p)
dev.off()

## plot the boxplot for the simulation study:composite for truncated-T ##
load("data/simulation_study_truncT_results_final_2000.RData",e1<-new.env())
data = cbind(rep(1:nrow(e1$par.truncT),times=sapply(e1$est.mat.list[[1]],nrow)),1,do.call(rbind,e1$est.mat.list[[1]]))
data = rbind(data,cbind(rep(1:nrow(e1$par.truncT),times=sapply(e1$est.mat.list[[2]],nrow)),2,do.call(rbind,e1$est.mat.list[[2]])))
data = rbind(data,cbind(rep(1:nrow(e1$par.truncT),times=sapply(e1$est.mat.list[[3]],nrow)),3,do.call(rbind,e1$est.mat.list[[3]])))
data = as.data.frame(data);names(data) = c("id","thres","hat(lambda)","hat(vartheta)","hat(nu)")
data = subset(data[,1:4],id %in% 1:2)
str(data)

# Reshape the data to a long format
data_long <- pivot_longer(data, -c(id, thres), names_to = "variable", values_to = "value")
data_long$facet = factor(data_long$variable)
par.truncT = as.data.frame(e1$par.truncT[1:2,1:2])
colnames(par.truncT) = c("hat(lambda)","hat(vartheta)")
par.truncT$id = 1:nrow(par.truncT)
par.truncT = rbind(par.truncT,par.truncT,par.truncT)
par.truncT$thres = rep(1:3,each=nrow(par.truncT)/3)
par.truncT_long <- pivot_longer(par.truncT, -c(id,thres), names_to = "variable", values_to = "value")
par.truncT_long$facet = factor(par.truncT_long$variable)

p <- ggplot(data_long, aes(x = factor(id), y = value,fill=factor(thres,labels=c("98%","95%","90%")))) +
  geom_violin(position = position_dodge(width=1),draw_quantiles = c(0.975,0.5,0.025),width=1.5) + 
  #geom_boxplot(position = position_dodge(width=0.75),width=0.4, color="grey", alpha=1) +
  geom_point(data=par.truncT_long,aes(x=factor(id),y=value),color="black",size=1,
  position=position_dodge(width = 1)) +
  facet_wrap(~ facet, scales = "free",ncol=2,nrow=1,labeller = label_parsed) +
  labs(x = "Cases",
       y = "Value") + 
  theme(axis.text = element_text(size = 14),
        strip.text = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_text(size = 16),legend.position = "none")

pdf("figures/simulation_est_boxplots_final_truncT.pdf",width=5*2,height=4,onefile=TRUE)
p
dev.off()


######################################################################
## plot the extremal coef for the application ########################
######################################################################
library(Matrix)
load("data/data_application.RData")
load("data/application_results_new_2.RData",e<-new.env())
e$results2$par;e$results2$value;e$results2$time
e$results4$par;e$results4$value;e$results4$time

## estimate the CI using jackknife ##
est.mat = matrix(unlist(lapply(e$results22, function(x){x$par})),ncol=length(e$results2$par),byrow=TRUE)
est.mat.BR = matrix(unlist(lapply(e$results42, function(x){x$par})),ncol=length(e$results4$par),byrow=TRUE)
est.jack = nrow(est.mat) * est.mat - (nrow(est.mat) - 1) * matrix(e$results2$par,nrow=nrow(est.mat),ncol=ncol(est.mat),byrow=TRUE)
est.jack.BR = nrow(est.mat.BR) * est.mat.BR - (nrow(est.mat.BR) - 1) * matrix(e$results4$par,nrow=nrow(est.mat.BR),ncol=ncol(est.mat.BR),byrow=TRUE)

est.mean.jack = colMeans(est.jack)
est.mean.jack.BR = colMeans(est.jack.BR)

est.sd.jack = apply(est.jack,2,sd)/sqrt(nrow(est.jack))
est.sd.jack.BR = apply(est.jack.BR,2,sd)/sqrt(nrow(est.jack.BR))

## compute the extremal coefficients ##
par.list.BR = alpha2delta(list(vario.func(e$loc.sub.trans,e$results4$par[1:2]),rep(0,ncol(distmat))))
par.list = list(vario.func(e$loc.sub.trans,e$results2$par[1:2]))
par.list[[2]] = alpha.func(par=c(e$results2$par[-c(1:2)]),b.mat=e$basis)
par.list = alpha2delta(par.list)

pairs = comb_n(1:ncol(distmat),2)

empirical.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=apply(pairs,2,empirical_extcoef,data=maxima.frechet),symmetric = TRUE,dimnames=NULL)

fitted.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x];V_bi_logskew(c(1,1),list(par.list[[1]][x,x],par.list[[2]][x]),alpha.para=FALSE)},mc.cores=5,mc.set.seed = FALSE)),symmetric = TRUE,dimnames=NULL) 

fitted.extcoef.BR.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x];V_bi_logskew(c(1,1),list(par.list.BR[[1]][x,x],c(0,0),alpha.para=FALSE))},mc.cores=5,mc.set.seed = FALSE)),symmetric = TRUE,dimnames=NULL)

diag(fitted.extcoef.mat) <- diag(fitted.extcoef.BR.mat) <- diag(empirical.extcoef.mat) <- 1



sqrt(mean((empirical.extcoef.mat - fitted.extcoef.mat)^2))
sqrt(mean((empirical.extcoef.mat - fitted.extcoef.BR.mat)^2))
diff.mat = abs(empirical.extcoef.mat - fitted.extcoef.mat) - abs(empirical.extcoef.mat - fitted.extcoef.BR.mat) 

diff.col.sums = unlist(lapply(1:ncol(distmat),function(i){mean(diff.mat[i,]<0)}))

sum(diff.col.sums > 0.5)

idx.centers = c(apply(maxima.frechet[which(rowmeans(maxima.frechet)>10),],1,function(x) which.min(x)),apply(maxima.frechet[which(rowmeans(maxima.frechet)>10),],1,function(x) which.max(x)))



save(empirical.extcoef.mat,fitted.extcoef.mat,fitted.extcoef.BR.mat,idx.centers,file="data/extcoef_application_RedSea_MaxStable.RData")
load("data/maps_RedSea.RData")

p1 <- list()
brks.emp <- c(1.2,1.5,1.7,1.75,1.8,1.85,1.9,1.95)#round(quantile(2-emp.extcoef.mat1@x,probs=c(0.005,0.05,0.1,0.3,0.5,0.8,0.9),na.rm=TRUE),4)
ewbreaks <- seq(34.4,41.6,2.4)
nsbreaks <- seq(15.6,26.4,3.6)
ewlabels <- unlist(lapply(ewbreaks, function(x) paste(" ",abs(x), "ºE")))
nslabels <- unlist(lapply(nsbreaks, function(x) paste(" ",x, "ºN")))
basis.centers.geo <- c(idx.centers,sample(1:1043,100))
basis.centers.geo <- basis.centers.geo[order(loc.sub[basis.centers.geo,2])] 
idx.centers = e$idx.centers
loc.df = data.frame(x = loc.sub[idx.centers,1],y = loc.sub[idx.centers,2],label=rep(".",length(idx.centers)))
for(i in 1:length(basis.centers.geo)){
    idx.center <- basis.centers.geo[i]
    data.df <- data.frame(lon=round(loc.sub[,1],5),lat=round(loc.sub[,2],5),
                emp=empirical.extcoef.mat[,idx.center],br=fitted.extcoef.BR.mat[,idx.center],
                sbr=fitted.extcoef.mat[,idx.center])
    data.df[idx.center,-c(1:2)] = 1
    p1[[i]]<-ggmap(map) +
    geom_tile(data=data.df,aes(x=lon,y=lat,fill=emp),alpha=0.8) + 
    colorspace::scale_fill_continuous_divergingx("RdYlBu",limits=c(1,2),mid=exp(1.6),alpha=0.8,name=expression(hat(theta)[2]),trans="exp") +
    coord_fixed() + stat_contour(data=data.df,aes(x=lon,y=lat,z=br),breaks = brks.emp,colour = "black",linetype="dashed") +
    stat_contour(data=data.df,aes(x=lon,y=lat,z=sbr),breaks = brks.emp,colour = "black") + labs(x="Longitude", y="Latitude") + 
    # stat_contour(data=data.df,aes(x=lon,y=lat,z=emp),breaks = brks.emp,colour = "black",linetype="dotted") +
    scale_x_continuous(breaks = ewbreaks, labels = ewlabels,expand=c(0,0),limits = xlim) + 
    scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = ylim) +
    geom_point(data=loc.df,aes(x = x, y = y), size = 2, color = "black") +
    theme(axis.text = element_text(size=10), 
                            strip.text = element_text(size = 14),
                            axis.title.x = element_text(size=14), 
                            axis.title.y = element_text(size=14),
                            axis.text.y = element_text(angle=90,vjust=0.5,hjust=1),
                            legend.title = element_text(size=14),legend.position = "right",
                            plot.margin=unit(c(0.2,0.2,0,0.5),"cm")) 
}

for(i in 1:length(p1)){
    ggsave(paste0("figures/application/RedSea2/RedSea_extcoef_MS_",sprintf(i,fmt="%.3d"),".png"),p1[[i]],width=5,height=6,dpi=300)
}

p1 <- p2 <- p3 <- p5 <- list()
loc.sub = as.data.frame(loc.sub)
#diff.extcoef = c()
for(i in 1:length(idx.centers)){
    idx.center = idx.centers[i]
    pairs.idx = rbind(idx.center,1:ncol(distmat))[,-idx.center]
    pairs.idx.list = split(pairs.idx,col(pairs.idx))
    true.ext.coef <- fitted.extcoef.mat[idx.center,-idx.center]
    true.ext.coef.BR <- fitted.extcoef.BR.mat[idx.center,-idx.center]
    empirical.extcoef <- empirical.extcoef.mat[idx.center,-idx.center]

    brks = round(quantile(c(true.ext.coef.BR,true.ext.coef),probs=seq(0.05,0.95,length.out=8)),3)
    global_min <- min(min(empirical.extcoef),min(true.ext.coef), min(true.ext.coef.BR))
    global_max <- max(max(empirical.extcoef),max(true.ext.coef), max(true.ext.coef.BR))

    data <- data.frame( x = loc.sub[-idx.center,1],
                            y = loc.sub[-idx.center,2],
                            z = true.ext.coef)

    loc.df = data.frame(x = loc.sub[c(idx.centers[i],e$idx.centers),1],y = loc.sub[c(idx.centers[i],e$idx.centers),2],label=c("*",rep(".",5)))

    p1[[i]] <- ggplot() + geom_tile(data=data,aes(x=x,y=y,fill=z)) +
            scale_fill_distiller(palette="RdBu",limits=c(global_min,global_max)) +
            geom_contour(data=data,aes(x=x,y=y,z=z),colour="black",breaks=brks) + 
            labs( x = "Longitude", y = "Latitude") + geom_point(data=loc.df[-1,],aes(x = x, y = y), size = 2, color = "black") + 
            geom_point(data=loc.df[1,],aes(x = x, y = y),shape=8, size = 2, color = "black") +
            theme(axis.text = element_text(size=10), 
                            strip.text = element_text(size = 14),
                            axis.title.x = element_text(size=14), 
                            axis.title.y = element_text(size=14), 
                            plot.title = element_text(hjust = 0.5,size=14),legend.title = element_text(size=14),legend.position = "none") + coord_fixed() 


    data$z = true.ext.coef.BR
 
    p2[[i]] <- ggplot() + geom_tile(data=data,aes(x=x,y=y,fill=z)) +
            scale_fill_distiller(palette="RdBu",limits=c(global_min,global_max)) +
            geom_contour(data=data,aes(x=x,y=y,z=z),colour="black",breaks=brks) + 
            labs( x = "Longitude", y = "Latitude") + geom_point(data=loc.df[-1,],aes(x = x, y = y), size = 2, color = "black") + 
            geom_point(data=loc.df[1,],aes(x = x, y = y),shape=8, size = 2, color = "black") + 
            theme(axis.text = element_text(size=12), 
                            strip.text = element_text(size = 14),
                            axis.title.x = element_text(size=14), 
                            axis.title.y = element_text(size=14), 
                            plot.title = element_text(hjust = 0.5,size=14),legend.title = element_text(size=14),legend.position = "none") + coord_fixed() 

   
    data$z = empirical.extcoef

    p3[[i]] <- ggplot() + geom_tile(data=data,aes(x=x,y=y,fill=z)) +
            scale_fill_distiller(palette="RdBu",limits=c(global_min,global_max)) +
            # geom_contour(data=data,aes(x=x,y=y,z=z),colour="black",breaks=brks) + 
            geom_point(data=loc.df[-1,],aes(x = x, y = y), size = 2, color = "black") + 
            geom_point(data=loc.df[1,],aes(x = x, y = y),shape=8, size = 2, color = "black") +
            labs(x = "Longitude", y = "Latitude",fill=expression(hat(theta)[2])) + 
            theme(axis.text = element_text(size=12), 
                            strip.text = element_text(size = 14),
                            axis.title.x = element_text(size=14), 
                            axis.title.y = element_text(size=14), 
                            plot.title = element_text(hjust = 0.5,size=14),legend.title = element_text(size=14),legend.position = "right") + coord_fixed() 

    data$z = abs(true.ext.coef - empirical.extcoef) - abs(true.ext.coef.BR - empirical.extcoef)
    p5[[i]] <- ggplot(data, aes(x = x, y = y, fill=as.factor(z<0)))  + 
            geom_tile() +
            scale_fill_manual(values=c("red","blue")) +
            theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot",legend.position = "right") + 
            coord_fixed() + 
            labs(title = paste("Skewed BR is closer to the Empirical?",round(mean(data$z<0)*100,1),"%"), x = "Longitude", y = "Latitude",fill="Values") +
            geom_point(aes(x = loc.sub[idx.centers[i],1], y = loc.sub[idx.centers[i],2]),shape=1, size = 2, color = "black") +
            geom_point(aes(x = loc.sub[e$idx.centers,1], y = loc.sub[e$idx.centers,2]),shape=8, size = 2, color = "black")
}

pdf("figures/extcoef_application_2_2.pdf",width=3.5*3+1,height=5,onefile=TRUE)
for(i in 1:length(idx.centers)){
    grid.arrange(grobs=list(p1[[i]],p2[[i]],p3[[i]]),nrow=1,widths=c(3.5,3.5,4.35))
}
dev.off()

p4 <- list()
data = data.frame( x = loc.sub[,1],
                    y = loc.sub[,2],
                    z = diff.col.sums)
                    
p4[[1]] <- ggplot(data, aes(x = x, y = y, fill=as.factor(z>0.5)))  + 
                geom_tile() +
                scale_fill_manual(values=c("red","blue")) +
                theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + 
                coord_fixed() + 
                labs(title = paste("Skewed BR is closer to the Empirical?",round(mean(data$z>0.5)*100,1),"%"), x = "X", y = "Y",fill="Values") 

pdf("figures/extcoef_application_diff.pdf",width=5,height=5,onefile=TRUE)
p4[[1]]
dev.off()

## compare the fitted bivariate extremal coef vs the empirical extremal coef ##
data = data.frame(x=c(row(empirical.extcoef.mat)),y=c(col(empirical.extcoef.mat)),z=as.vector(empirical.extcoef.mat))
png("figures/extcoef_application_empirical.png",width=800,height=800)
p4[[2]] <- ggplot(data,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_distiller(palette="RdBu") + coord_fixed() + labs(title="Empirical Extremal Coef",x="X",y="Y")
p4[[2]]
dev.off()

data$z = as.vector(fitted.extcoef.mat) 
png("figures/extcoef_application_fitted.png",width=800,height=800)
p4[[3]] <- ggplot(data,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_distiller(palette="RdBu") + coord_fixed() + labs(title="Fitted Skewed Brown Resnick",x="X",y="Y")
p4[[3]]
dev.off()

data$z = as.vector(fitted.extcoef.BR.mat)
png("figures/extcoef_application_fitted_BR.png",width=800,height=800)
p4[[4]] <- ggplot(data,aes(x=x,y=y,z=z)) + geom_tile(aes(fill=z)) + scale_fill_distiller(palette="RdBu") + coord_fixed() + labs(title="Fitted Brown Resnick",x="X",y="Y")
p4[[4]]
dev.off()

data$z = as.vector(fitted.extcoef.BR.mat - fitted.extcoef.mat) 
png("figures/extcoef_application_diff_models.png",width=800,height=800)
p4[[6]] <- ggplot(data,aes(x=x,y=y,fill=z)) + geom_tile() + scale_fill_distiller(palette="RdBu") + coord_fixed() + labs(title="Differences",x="X",y="Y",fill="Values")
p4[[6]]
dev.off()

save(idx.centers,est.sd.jack,est.sd.jack.BR,est.jack,est.jack.BR,p1,p2,p3,p4,par.list,par.list.BR,empirical.extcoef.mat,fitted.extcoef.BR.mat,fitted.extcoef.mat,file="data/plot_application_2.RData")




















































