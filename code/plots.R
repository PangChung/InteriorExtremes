rm(list=ls())
library(ggplot2)
library(parallel)
library(gridExtra)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
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

sigma = cov.func(coord,c(4,1))
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


library(directlabels)
idx.center = c(16,16)
idx.case = 3
idx.center = which.min(abs(coord[,1] - idx.center[1]) + abs(coord[,2] - idx.center[2]))
ind.idx.center = all.pairs[1,] == idx.center |  all.pairs[2,] == idx.center
true.ext.coef <- unlist(mclapply(all.pairs.list[ind.idx.center],true_extcoef,par=par.list.1[[idx.case]],model="logskew1",mc.cores=5,mc.set.seed = FALSE))
idx = apply(all.pairs[,ind.idx.center],2,function(x){x[x!=idx.center]})
data <- data.frame( x = coord[idx,1],
                    y = coord[idx,2],
                    z = true.ext.coef )
p1 <- ggplot(data, aes(x = x, y = y, z=z)) +
    geom_contour(bins=5) +
    geom_dl(aes(label=after_stat(nlevel)),method="bottom.pieces", 
             stat="contour") + 
    scale_color_gradient(low="blue",high = "red") +
    theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
    labs(title = paste("Bivariate Extremal Coef"), x = "X", y = "Y")
p1
p2 <- ggplot(data, aes(x = x, y = y, fill=z)) +
    geom_tile() + scale_fill_distiller(palette="YlOrRd") +
    theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
    labs(title = paste("Bivariate Extremal Coef."), x = "X", y = "Y")
p2
