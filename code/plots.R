rm(list=ls())
library(ggplot2)
library(parallel)
library(gridExtra)
library(directlabels)
library(Rfast)
library(RColorBrewer)
library(mev)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
set.seed(1)
para.alpha = rbind(c(0,0),c(-1,-2),c(-1,1)) ## slant parameter for skewed norm model ##
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
basis <- sapply(1:(ncol(para.alpha)+1),function(x){y <- rep(0,nrow(coord));y[idx[x,]] <- c(-2,2);y})
basis[,1] = rep(0,d^2);basis[1:floor(d^2/2),1] = 0.1; basis[(d^2-floor(d^2/2)+1):d^2,1] = -0.1
basis.list[[2]] <- basis

all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))

sigma = cov.func(coord,c(8,1.5))
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



idx.center = c(16,16)
idx.center = which.min(abs(coord[,1] - idx.center[1]) + abs(coord[,2] - idx.center[2]))
ind.idx.center = all.pairs[1,] == idx.center |  all.pairs[2,] == idx.center
idx = apply(all.pairs[,ind.idx.center],2,function(x){x[x!=idx.center]})
idx.case = 3;p1.list <- p2.list <- list()
for(idx.case in 1:length(par.list.1)){
    true.ext.coef <- unlist(mclapply(all.pairs.list[ind.idx.center],true_extcoef,par=par.list.1[[idx.case]],model="logskew1",mc.cores=5,mc.set.seed = FALSE))
    data <- data.frame( x = coord[idx,1],
                        y = coord[idx,2],
                        z = true.ext.coef )

    brks = round(quantile(data$z,probs=seq(0.001,0.999,length.out=5)),4)

    p1 <- ggplot(data, aes(x = x, y = y, z=z))  + 
        geom_tile(aes(fill=z)) +
        scale_fill_distiller(palette="RdBu") +
        geom_contour(colour="black",breaks=brks) + 
        geom_dl(aes(label=..level..),method="bottom.pieces",breaks=brks, 
                stat="contour") + 
        theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
        labs(title = paste("Bivariate Extremal Coef"), x = "X", y = "Y")
    
    p1.list[[idx.case]] <- p1
    
    true.ext.coef <- unlist(mclapply(all.pairs.list[ind.idx.center],true_extcoef,par=par.list.2[[idx.case]],model="logskew1",mc.cores=5,mc.set.seed = FALSE))
    data <- data.frame( x = coord[idx,1],
                        y = coord[idx,2],
                        z = true.ext.coef )

    brks = round(quantile(data$z,probs=seq(0.001,0.999,length.out=5)),4)
    p2 <- ggplot(data, aes(x = x, y = y, z=z))  + 
        geom_tile(aes(fill=z)) +
        scale_fill_distiller(palette="RdBu") +
        geom_contour(colour="black",breaks=brks) + 
        geom_dl(aes(label=..level..),method="bottom.pieces",breaks=brks, 
                stat="contour") + 
        theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
        labs(title = paste("Bivariate Extremal Coef"), x = "X", y = "Y")
    p2.list[[idx.case]] <- p2
}

true.ext.coef <- unlist(mclapply(all.pairs.list[ind.idx.center],true_extcoef,par=list(par.list.1[[1]][[1]],rep(0,nrow(coord))),model="logskew1",mc.cores=5,mc.set.seed = FALSE))
data <- data.frame( x = coord[idx,1],
                        y = coord[idx,2],
                        z = true.ext.coef )

brks = round(quantile(data$z,probs=seq(0.001,0.999,length.out=5)),4)
p2 <- ggplot(data, aes(x = x, y = y, z=z))  + 
        geom_tile(aes(fill=z)) +
        scale_fill_distiller(palette="RdBu") +
        geom_contour(colour="black",breaks=brks) + 
        geom_dl(aes(label=..level..),method="bottom.pieces",breaks=brks, 
                stat="contour") + 
        theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
        labs(title = paste("Bivariate Extremal Coef"), x = "X", y = "Y")
p2.list[[1]] <- p2

pdf("figures/extcoef_final_logskew.pdf",width=5*3,height=5*2,onefile=TRUE)
# grid.arrange(grobs=p1.list,ncol=3)
# grid.arrange(grobs=p2.list,ncol=3)
grid.arrange(grobs=c(p1.list,p2.list),ncol=3,nrow=2)
dev.off()

p.list = list()
sigma.22 = 1
rho = 1:9/10*sigma.22
BR.values = unlist(lapply(rho,function(x){V_bi_logskew(c(1,1),delta=c(0,0),sigma=matrix(c(sigma.22,x,x,sigma.22),2,2))}))
for(i in 1:length(rho)){
    r = sqrt(min(eigen(matrix(c(1,rho[i]/sigma.22,rho[i]/sigma.22,1),2))$values))
    delta = seq(-r,r,length.out=100)
    delta.grid = as.matrix(expand.grid(delta,delta))
    delta.grid.list <- split(delta.grid,row(delta.grid))

    idx.valid = apply(delta.grid,1,function(x){sum(x^2) < r^2}) # & abs(diff(x))<sqrt(2-2*rho)
    
    values <- unlist(lapply(delta.grid.list[idx.valid],function(x){V_bi_logskew(c(1,1),delta=x,sigma=matrix(c(sigma.22,rho[i],rho[i],sigma.22),2,2))}))

    data = data.frame(x=delta.grid[idx.valid,1],y=delta.grid[idx.valid,2],z=values)
    
    data2 = data.frame(x=0,y=0,z=BR.values[i])
    p.list[[i]] <- ggplot(data) + geom_point(aes(x=x,y=y,color=z)) + scale_color_distiller(palette="RdBu") + ggtitle(paste("rho",rho[i])) + coord_fixed() + theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + geom_point(data=data2,aes(x=x, y = y), color = "black") + ggtitle(paste("Max:",round(max(data$z),4),"BR:",round(BR.values[i],4)))
}

pdf("figures/bivariate_extcoef_rho.pdf",width=5*3,height = 5*3,onefile = TRUE)
grid.arrange(grobs=p.list,ncol=3,nrow=3)
dev.off()

## plot the true extremal coef for the truncated extremal t model ##
d=1000
coord.trunc = as.matrix(expand.grid(0,(1:d)/d))
all.pairs.trunc = rbind(1,2:d)
all.pairs.list.trunc = split(all.pairs.trunc,col(all.pairs.trunc))
par.truncT.list = list(cov.func(coord.trunc,c(0.5,1.5)),nu=2)

true.ext.truncT <- unlist(lapply(all.pairs.list.trunc,true_extcoef,par=par.truncT.list,model="truncT2"))

true.ext.t <- unlist(lapply(all.pairs.list.trunc,function(id) mev::expme(z=rep(1,2),par=list(Sigma=par.truncT.list[[1]][id,id],df=2),model="xstud") ))

plot(x=coord.trunc[-1,2],y=true.ext.truncT,type="l",col="black",ylim=c(1,2),xlab="coordinate",ylab="Bivariate extremal coeffient")
lines(x=coord.trunc[-1,2],y=true.ext.t,col="red")


## plot the extremal coef for the application ##
load("data/data_application.RData")
load("data/application_results2.RData",e<-new.env())
par.list.BR = alpha2delta(list(vario.func(loc.sub.trans,e$results4$par[1:2]),rep(0,ncol(distmat))))
loc.sub.trans = apply(loc.sub.trans,2,function(x) x-mean(x))
par.list = list(vario.func(loc.sub.trans,e$results2$par[1:2]))
par.list[[2]] = alpha.func(par=e$results2$par[-c(1:2)],b.mat=e$basis/sqrt(diag(par.list[[1]])))
par.list = alpha2delta(par.list)
idx.centers = c(195,536,944)
#idx.centers = c(212, 508, 843)
#idx.centers = which(rank(colSums(distmat)) %in% c(1,100,200))
#idx.centers = 1:ncol(distmat)
p1 <- p2 <- p3 <- p4 <- list()
#diff.extcoef = c()
for(i in 1:length(idx.centers)){
    idx.center = idx.centers[i]
    pairs.idx = rbind(idx.center,1:ncol(distmat))[,-idx.center]
    pairs.idx.list = split(pairs.idx,col(pairs.idx))
    true.ext.coef <- unlist(lapply(pairs.idx.list,function(x){V_bi_logskew(c(1,1),delta = par.list[[2]][x],sigma=par.list[[1]][x,x])}))
    true.ext.coef.BR <- unlist(lapply(pairs.idx.list,function(x){V_bi_logskew(c(1,1),delta = par.list.BR[[2]][x],sigma=par.list[[1]][x,x])}))
    empirical.extcoef <- apply(pairs.idx,2,empirical_extcoef,data=maxima.frechet)
#     diff.extcoef[i] = mean(abs(true.ext.coef.BR - empirical.extcoef) - abs(true.ext.coef - empirical.extcoef))
#     print(i)
# }
    data <- data.frame( x = loc.sub[-idx.center,1],
                            y = loc.sub[-idx.center,2],
                            z = true.ext.coef)

    brks = round(quantile(data$z,probs=seq(0.001,0.999,length.out=5)),3)
    
    p1[[i]] <- ggplot(data, aes(x = x, y = y, z=z))  + 
            geom_tile(aes(fill=z)) +
            scale_fill_distiller(palette="RdBu") +
            geom_contour(colour="black",breaks=brks) + 
            # geom_dl(aes(label=after_stat(level)),method="bottom.pieces",breaks=brks, stat="contour") + 
            theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
            labs(title = paste("Skewed Brown-Resnick"), x = "X", y = "Y") 


    data$z = true.ext.coef.BR

    brks = round(quantile(data$z,probs=seq(0.001,0.999,length.out=5)),3)

    p2[[i]] <- ggplot(data, aes(x = x, y = y, z=z))  + 
            geom_tile(aes(fill=z)) +
            scale_fill_distiller(palette="RdBu") +
            geom_contour(colour="black",breaks=brks) + 
            # geom_dl(aes(label=after_stat(level)),method="bottom.pieces",breaks=brks, stat="contour") + 
            theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
            labs(title = paste("Brown-Resnick"), x = "X", y = "Y")

   
    data$z = empirical.extcoef

    brks = round(quantile(data$z,probs=seq(0.001,0.999,length.out=5)),3)

    p3[[i]] <- ggplot(data, aes(x = x, y = y, z=z))  + 
            geom_tile(aes(fill=z)) +
            scale_fill_distiller(palette="RdBu") +
            #geom_contour(colour="black",breaks=brks) + 
            #geom_dl(aes(label=after_stat(level)),method="bottom.pieces",breaks=brks, stat="contour") + 
            theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
            labs(title = paste("Empirical"), x = "X", y = "Y") 

    
    data$z = abs(true.ext.coef.BR - empirical.extcoef) - abs(true.ext.coef - empirical.extcoef) 

    #brks = round(quantile(data$z,probs=seq(0.001,0.999,length.out=5)),4)

    # p4[[i]] <- ggplot(data, aes(x = x, y = y, z=z))  + 
    #         geom_tile(aes(fill=z)) +
    #         scale_fill_distiller(palette="RdBu") +
    #         theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
    #         labs(title = paste("Differences"), x = "X", y = "Y") 
    p4[[i]] <- ggplot(data, aes(x = x, y = y, fill=factor(z>0)))  + 
                geom_tile() +
                scale_fill_manual(values=c("red","blue")) +
                theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + 
                coord_fixed() + 
                labs(title = paste("Differences",round(mean(data$z),3)), x = "X", y = "Y",fill="Values") 

}

layout = matrix(1:(4*length(idx.centers)),ncol=length(idx.centers),byrow=TRUE)
grid.arrange(grobs=c(p1,p2,p3,p4),layout_matrix=layout)


## compare the fitted bivariate extremal coef vs the empirical extremal coef ##
pairs = comb_n(1:ncol(distmat),2)
empirical.extcoef.all <- apply(pairs,2,empirical_extcoef,data=maxima.frechet)
fitted.extcoef <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x];V_bi_logskew(c(1,1),delta = par.list[[2]][x],sigma=par.list[[1]][x,x])},mc.cores=5,mc.set.seed = FALSE))
fitted.extcoef.BR <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x];V_bi_logskew(c(1,1),delta = c(0,0),sigma=par.list[[1]][x,x])},mc.cores=5,mc.set.seed = FALSE))

sqrt(mean((empirical.extcoef.all - fitted.extcoef)^2))
sqrt(mean((empirical.extcoef.all - fitted.extcoef.BR)^2))
save(idx.centers,p1,p2,p3,p4,par.list,par.list.BR,empirical.extcoef.all,file="data/plot_application.RData")
