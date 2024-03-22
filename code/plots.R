rm(list=ls())
library(ggplot2)
library(parallel)
library(gridExtra)
library(directlabels)
library(RColorBrewer)
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
    scale_fill_gradient2(low = "blue",mid="white" ,high = "red",limits=c(-0.5,0.5)) +
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1),plot.title = element_text(hjust = 0.5)) + coord_fixed()
    data = data.frame(x=coord[,1],y=coord[,2],z=par.list.2[[i]][[2]])
    print(range(data$z))
    p.list2[[i]] <- ggplot(data, aes(x = x, y = y, fill = z)) +
    geom_tile() + ggtitle(paste(c("alpha:",para.alpha[i,]),collapse = " ")) + 
    scale_fill_gradient2(low = "blue",mid="white" ,high = "red",limits=c(-0.5,0.5)) +
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
        scale_color_gradient(low="blue",high = "red") +
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
        scale_color_gradient(low="blue",high = "red") +
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
        scale_color_gradient(low="blue",high = "red") +
        theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
        labs(title = paste("Bivariate Extremal Coef"), x = "X", y = "Y")
p2.list[[1]] <- p2

grid.arrange(grobs=p1.list,ncol=3)
grid.arrange(grobs=p2.list,ncol=3)
grid.arrange(grobs=c(p1.list,p2.list),ncol=3,nrow=2)

## plot the true extremal coef for the truncated extremal t model ##
d=5
coord.trunc = as.matrix(expand.grid(1:d,1:d))
all.pairs.trunc = combn(1:nrow(coord.trunc),2)
all.pairs.list.trunc = split(all.pairs.trunc,col(all.pairs.trunc))
idx.center = c(d/2,d/2)
idx.center = which.min(abs(coord.trunc[,1] - idx.center[1]) + abs(coord.trunc[,2] - idx.center[2]))
ind.idx.center = all.pairs.trunc[1,] == idx.center |  all.pairs.trunc[2,] == idx.center
idx = apply(all.pairs.trunc[,ind.idx.center],2,function(x){x[x!=idx.center]})
par.truncT.list = list(cov.func(coord.trunc,c(5,1.5)),nu=2)
T_j.list = a_fun(par.truncT.list,ncores=5)

## plot the true extremal coef for the truncated extremal t model ##
true.ext.coef1 <- true_extcoef(all.pairs.list.trunc[ind.idx.center],par.truncT.list,model="truncT1",T_j=T_j.list)
data <- data.frame( x = coord.trunc[idx,1],
                        y = coord.trunc[idx,2],
                        z = true.ext.coef)
brks = round(quantile(data$z,probs=seq(0.001,0.999,length.out=5)),4)
p3 <- ggplot(data, aes(x = x, y = y, z=z))  + 
        geom_tile(aes(fill=z)) +
        scale_fill_distiller(palette="RdBu") +
        geom_contour(colour="black",breaks=brks) + 
        geom_dl(aes(label=..level..),method="bottom.pieces",breaks=brks, 
                stat="contour") + 
        scale_color_gradient(low="blue",high = "red") +
        theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
        labs(title = paste("Bivariate Extremal Coef"), x = "X", y = "Y")
p3

true.ext.coef2 <- unlist(mclapply(all.pairs.list.trunc[ind.idx.center],true_extcoef,par=par.truncT.list,model="truncT2",mc.cores=5,mc.set.seed = FALSE))
data <- data.frame( x = coord.trunc[idx,1],
                        y = coord.trunc[idx,2],
                        z = true.ext.coef)

brks = round(quantile(data$z,probs=seq(0.001,0.999,length.out=5)),4)
p4 <- ggplot(data, aes(x = x, y = y, z=z))  + 
        geom_tile(aes(fill=z)) +
        scale_fill_distiller(palette="RdBu") +
        geom_contour(colour="black",breaks=brks) + 
        geom_dl(aes(label=..level..),method="bottom.pieces",breaks=brks, 
                stat="contour") + 
        scale_color_gradient(low="blue",high = "red") +
        theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
        labs(title = paste("Bivariate Extremal Coef"), x = "X", y = "Y")
p4

