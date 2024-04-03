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

#sigma = cov.func(diff.mat,c(8,1.5))
sigma = vario.func(coord,c(8,1))
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

idx.center = c(8,8)
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
        labs(title = paste("Skewed Brown-Resnick"), x = expression(s[1]), y = expression(s[2]),fill=expression(theta[2]))
    
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
        labs(title = paste("Skewed Brown-Resnick"), x = expression(s[1]), y = expression(s[2]),fill=expression(theta[2]))
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
        labs(title = paste("Brown-Resnick"), x = expression(s[1]), y = expression(s[2]),fill=expression(theta[2]))
p2.list[[1]] <- p2

pdf("figures/extcoef_final_logskew2.pdf",width=5*3,height=5*2,onefile=TRUE)
# grid.arrange(grobs=p1.list,ncol=3)
# grid.arrange(grobs=p2.list,ncol=3)
grid.arrange(grobs=c(p1.list,p2.list),ncol=3,nrow=2)
dev.off()

p.list = list()
sigma.22 = 10
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


p.list = list()
sigma.22 = 10
rho = seq(0.1,sigma.22-0.1,length.out=200)
BR.values = unlist(lapply(rho,function(x){V_bi_logskew(c(1,1),delta=c(0,0),sigma=matrix(c(sigma.22,x,x,sigma.22),2,2))}))
alpha = seq(0,5,length.out=100)
para.grid = as.matrix(expand.grid(alpha,rho))
values <- lapply(split(para.grid,row(para.grid)),function(x){par.list = alpha2delta(list(matrix(c(sigma.22,x[2],x[2],sigma.22),2,2),alpha=c(x[1],-x[1])));V_bi_logskew(c(1,1),delta=par.list[[2]],sigma=par.list[[1]])})

data = data.frame(x=para.grid[,1],y=para.grid[,2],z=unlist(values))
data2 = data.frame(x=0,y=para.grid[,2],z=BR.values)

brks = round(quantile(data$z,probs=seq(0.01,0.99,length.out=4)),2)
p <- ggplot(data, aes(x = x, y = y, z=z))  + 
        geom_tile(aes(fill=z)) +
        scale_fill_distiller(palette="RdBu",limits=c(1,2)) +
        geom_contour(colour="black",breaks=brks) + 
        geom_dl(aes(label=..level..),method="bottom.pieces",breaks=brks,stat="contour") + 
        labs(title="Bivariate extremal coefficients",x=expression(alpha[1]),y=expression(sigma[12]),fill=expression(theta[2])) +
        theme(axis.text = element_text(size=10), 
                            axis.title.x = element_text(size=14), 
                            axis.title.y = element_text(size=14), 
                            plot.title = element_text(hjust = 0.5,size=14),legend.title = element_text(size=14))

pdf("figures/bivariate_extcoef_rho_alpha.pdf",width=5,height = 4,onefile = TRUE)
p
dev.off()

#val.mat = matrix(unlist(values),ncol=length(alpha),byrow=TRUE)

## plot the true extremal coef for the truncated extremal t model ##
sigma.22 = 1
df = c(2,5,10)
rho = seq(0.001,sigma.22-0.001,length.out=200)
data = NULL
for(i in 1:length(df)){
    par.truncT.list = lapply(rho,function(x){list(matrix(c(sigma.22,x,x,sigma.22),2,2),df[i])})
    true.ext.truncT <- unlist(lapply(par.truncT.list,V_truncT,x=c(1,1)))
    true.ext.t <- unlist(lapply(1:length(par.truncT.list),function(id) mev::expme(z=rep(1,2),par=list(Sigma=par.truncT.list[[id]][[1]],df=par.truncT.list[[id]][[2]]),model="xstud") ))
    data = rbind(data,data.frame(x=rho,truncT=true.ext.truncT,extT = true.ext.t,df=df[i]))
}
new_labels <- paste0("nu==",df)
data$df <- factor(data$df,labels=new_labels)
p <- ggplot(data) + geom_line(aes(x=x,y=truncT,color="Truncated extremal-t"),,linewidth=1) + 
geom_line(aes(x=x,y=extT,color="Extremal-t"),linewidth=1) + coord_cartesian(ylim = c(1, 2)) + geom_hline(yintercept=c(1,2),linetype="dashed") +
facet_wrap(~df,labeller =label_parsed) + theme(legend.position = "bottom") + labs(title="Bivariate extremal coefficient",x=expression(rho),y=expression(theta[2]),color="Model") +
theme(axis.text = element_text(size=10), 
            strip.text = element_text(size = 14), 
            axis.title.x = element_text(size=14), 
            axis.title.y = element_text(size=14), 
            plot.title = element_text(hjust = 0.5,size=14),legend.title = element_text(size=14),legend.text = element_text(size=14))

pdf("figures/bivariate_extcoef_truncT.pdf",width=8,height = 3,onefile = TRUE)
p
dev.off()

## plot the extremal coef for the application ##
load("data/data_application.RData")
load("data/application_results2_6.RData",e<-new.env())
e$results2$par
e$results4$par
par.list.BR = alpha2delta(list(vario.func(e$loc.sub.trans,e$results4$par[1:2]),rep(0,ncol(distmat))))
par.list = list(vario.func(e$loc.sub.trans,e$results2$par[1:2]))
par.list[[2]] = alpha.func(par=e$results2$par[-c(1:2)],b.mat=e$basis/sqrt(diag(par.list[[1]])))
par.list = alpha2delta(par.list)
pairs = comb_n(1:ncol(distmat),2)
library(Matrix)
empirical.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=apply(pairs,2,empirical_extcoef,data=maxima.frechet),symmetric = TRUE,dimnames=NULL)
fitted.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x];V_bi_logskew(c(1,1),delta = par.list[[2]][x],sigma=par.list[[1]][x,x])},mc.cores=5,mc.set.seed = FALSE)),symmetric = TRUE,dimnames=NULL) 
fitted.extcoef.BR.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x];V_bi_logskew(c(1,1),delta = c(0,0),sigma=par.list.BR[[1]][x,x])},mc.cores=5,mc.set.seed = FALSE)),symmetric = TRUE,dimnames=NULL)
diag(fitted.extcoef.mat) <- diag(fitted.extcoef.BR.mat) <- diag(empirical.extcoef.mat) <- 1
sqrt(mean((empirical.extcoef.mat - fitted.extcoef.mat)^2))
sqrt(mean((empirical.extcoef.mat - fitted.extcoef.BR.mat)^2))
diff.mat = abs(empirical.extcoef.mat - fitted.extcoef.mat) - abs(empirical.extcoef.mat - fitted.extcoef.BR.mat) 
#diff.col.sums = colMeans(diff.mat)
diff.col.sums = unlist(lapply(1:ncol(distmat),function(i){mean(diff.mat[i,]<0)}))
sum(diff.col.sums < 0)
#idx.centers = e$idx.centers
idx.centers = c(200,538,800)
diff.col.sums[idx.centers]
#idx.centers = which(rank(colSums(distmat)) %in% c(1,100,200))
#idx.centers = 1:ncol(distmat)
p1 <- p2 <- p3 <- p5 <- list()
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

    p1[[i]] <- ggplot(data, aes(x = x, y = y, z=z))  + 
            geom_tile(aes(fill=z)) +
            scale_fill_distiller(palette="RdBu",limits=c(global_min,global_max)) +
            geom_contour(colour="black",breaks=brks) + 
            # geom_dl(aes(label=after_stat(level)),method="bottom.pieces",breaks=brks, stat="contour") + 
            theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
            labs(title = paste("Skewed Brown-Resnick"), x = "X", y = "Y") 


    data$z = true.ext.coef.BR

    p2[[i]] <- ggplot(data, aes(x = x, y = y, z=z))  + 
            geom_tile(aes(fill=z)) +
            scale_fill_distiller(palette="RdBu",limits=c(global_min,global_max)) +
            geom_contour(colour="black",breaks=brks) + 
            # geom_dl(aes(label=after_stat(level)),method="bottom.pieces",breaks=brks, stat="contour") + 
            theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
            labs(title = paste("Brown-Resnick"), x = "X", y = "Y")

   
    data$z = empirical.extcoef

    p3[[i]] <- ggplot(data, aes(x = x, y = y, z=z))  + 
            geom_tile(aes(fill=z)) +
            scale_fill_distiller(palette="RdBu",limits=c(global_min,global_max)) +
            #geom_contour(colour="black",breaks=brks) + 
            #geom_dl(aes(label=after_stat(level)),method="bottom.pieces",breaks=brks, stat="contour") + 
            theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + coord_fixed() + 
            labs(title = paste("Empirical"), x = "X", y = "Y") 

    data$z = abs(true.ext.coef - empirical.extcoef) - abs(true.ext.coef.BR - empirical.extcoef)
    p5[[i]] <- ggplot(data, aes(x = x, y = y, fill=as.factor(z<0)))  + 
            geom_tile() +
            scale_fill_manual(values=c("red","blue")) +
            theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + 
            coord_fixed() + 
            labs(title = paste("Skewed BR is closer to the Empirical?",round(mean(data$z<0)*100,1),"%"), x = "X", y = "Y",fill="Values")
}


pdf("figures/extcoef_application.pdf",width=5*4,height=5*length(idx.centers),onefile=TRUE)
layout = matrix(1:(4*length(idx.centers)),ncol=length(idx.centers),byrow=TRUE)
grid.arrange(grobs=c(p1,p2,p3,p5),layout_matrix=layout)
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

data$z = par.list[[2]]
data.new = data.frame(x1=loc.sub[idx.centers,1],y1=loc.sub[idx.centers,2])
png("figures/extcoef_application_delta.png",width=800,height=800)
p4[[9]] <- ggplot(data, aes(x = x, y = y))  + 
                geom_tile(aes(fill=z)) +
                scale_fill_distiller(palette="RdBu") +
                theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + 
                coord_fixed() + 
                labs(title = paste("Delta"), x = "X", y = "Y",fill="Values") #
p4[[9]] <- p4[[9]] + geom_point(data=data.new,aes(x=x1,y=y1),color="black",size=2,shape=20)
p4[[9]]
dev.off()

data$z = delta2alpha(par.list)[[2]]
data.new = data.frame(x1=loc.sub[idx.centers,1],y1=loc.sub[idx.centers,2])
png("figures/extcoef_application_alpha.png",width=800,height=800)
p4[[9]] <- ggplot(data, aes(x = x, y = y))  + 
                geom_tile(aes(fill=z)) +
                scale_fill_distiller(palette="RdBu") +
                theme(plot.title = element_text(hjust = 0.5), plot.title.position = "plot") + 
                coord_fixed() + 
                labs(title = paste("Alpha"), x = "X", y = "Y",fill="Values") #
p4[[9]] <- p4[[9]] + geom_point(data=data.new,aes(x=x1,y=y1),color="black",size=2,shape=20)
p4[[9]]
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

data$z = as.vector(diff.mat) 
png("figures/extcoef_application_diff.png",width=800,height=800)
p4[[5]] <- ggplot(data,aes(x=x,y=y,fill=z)) + geom_tile() + scale_fill_distiller(palette="RdBu")  + coord_fixed() + labs(title="Differences",x="X",y="Y",fill="Values")
p4[[5]]
dev.off()

data$z = as.vector(fitted.extcoef.BR.mat - fitted.extcoef.mat) 
png("figures/extcoef_application_diff_models.png",width=800,height=800)
p4[[6]] <- ggplot(data,aes(x=x,y=y,fill=z)) + geom_tile() + scale_fill_distiller(palette="RdBu") + coord_fixed() + labs(title="Differences",x="X",y="Y",fill="Values")
p4[[6]]
dev.off()

data$z = as.vector(fitted.extcoef.mat - empirical.extcoef.mat) 
png("figures/extcoef_application_diff_skew_BR.png",width=800,height=800)
p4[[7]] <- ggplot(data,aes(x=x,y=y,fill=z)) + geom_tile() + scale_fill_distiller(palette="RdBu") + coord_fixed() + labs(title="Differences",x="X",y="Y",fill="Values")
p4[[7]]
dev.off()

data$z = as.vector(fitted.extcoef.BR.mat - empirical.extcoef.mat) 
png("figures/extcoef_application_diff_BR.png",width=800,height=800)
p4[[8]] <- ggplot(data,aes(x=x,y=y,fill=z)) + geom_tile() + scale_fill_distiller(palette="RdBu") + coord_fixed() + labs(title="Differences",x="X",y="Y",fill="Values")
p4[[8]]
dev.off()


save(idx.centers,p1,p2,p3,p4,par.list,par.list.BR,empirical.extcoef.mat,fitted.extcoef.BR.mat,fitted.extcoef.mat,file="data/plot_application.RData")















