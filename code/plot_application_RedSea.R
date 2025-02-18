rm(list=ls())
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(partitions)
library(Rfast)
library(evd)
library(matrixStats)
library(Matrix)
library(sf)
library(ggplot2)
library(ggmap)
library(RColorBrewer)
library(Matrix)
# load the data##
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
load("data/data_application.RData")
load("data/Trends_fits.RData")

load("data/maps_RedSea.RData")
pairs <- comb_n(1:nrow(loc.sub),2)

# register_google(key=system("echo $g_key",intern = TRUE))
# ylim = c(12,30)
# xlim = c(32,44)
# xy.center = c(mean(xlim),mean(ylim))
# map <- get_googlemap(xy.center,zoom=5,maptype = "terrain",style = "feature:all|element:all|saturation:-100|lightness:50")
# ggmap(map) + theme(plot.title = element_text(hjust=0.5)) + ggtitle("Red Sea") + xlim(xlim) + ylim(ylim) 

# save(map,loc.sub,loc.sub.trans,xy.center,xlim,ylim,file="data/maps_RedSea.RData")

idx.centers = unlist(lapply(quantile(loc.sub.trans[,1],seq(0.1,0.9,length.out=5)),function(x){ idx = abs(loc.sub.trans[,1] - x) < 5; which(idx)[which.min(abs(loc.sub.trans[idx,2] - median(loc.sub.trans[idx,2])))]}))
basis <- sapply(idx.centers,function(x){y=dnorm(distmat[x,],mean=0,sd=ncol(distmat)*2);y=y-mean(y);y/sqrt(sum(y^2))})

pairs <- comb_n(1:nrow(loc.sub),2)

load("data/application_results.RData")
fit.results <- subset(data.RedSea,type=="full")
fitted.extcoef.mat.list <- list()
for(i in 1:4){
    alpha <- alpha.func(as.numeric(fit.results[i,5:9]),basis)
    cov.mat <- vario.func2(loc.sub.trans,as.numeric(fit.results[i,1:4]),ncores=5)
    par.list <- alpha2delta(list(cov.mat,alpha))
    fitted.extcoef <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x];V_bi_logskew(c(1,1),list(par.list[[1]][x,x],par.list[[2]][x]),alpha.para=FALSE)},mc.cores=5,mc.set.seed = FALSE))
    range(fitted.extcoef)
    fitted.extcoef.mat.list[[i]] <- sparseMatrix(i=pairs[2,],j=pairs[1,],x=fitted.extcoef,symmetric = TRUE,dimnames=NULL) 
}
save(fitted.extcoef.mat.list,file="data/fitted_extcoef_RedSea.RData")


data <- residuals2
data[data<0] = 0
data <- apply(data,2,function(x) {ind = x>0;x[ind] = qgpd(rank(x[ind])/(sum(ind)+1),1,1,1); x })
data = apply(data,1,function(x){list(which(x>0),x[x>0])})
idx.data = which(sapply(data,function(x) length(x[[2]]))>0)
data = data[idx.data]
data.sum = sapply(data,function(x) mean(x[[2]]))
data.max = sapply(data,function(x) max(x[[2]]))

data.fit = data[data.sum>quantile(data.sum,0.75)]
len.row <- unlist(lapply(1:length(data.fit),function(i){length(data.fit[[i]][[1]])}))
data.pareto.mat <- sparseMatrix(i=rep(1:length(data.fit),times=len.row),j=unlist(lapply(1:length(data.fit),function(i){data.fit[[i]][[1]]})),x=unlist(lapply(1:length(data.fit),function(i){data.fit[[i]][[2]]})),dimnames=NULL,symmetric = FALSE)
data.pareto.mat <- as.matrix(data.pareto.mat)

empirical.extcoef <- function(data){
    x = data[,1]
    y= data[,2]
    u = 50
    return( sum(x>u & y>u)/(sum(x>u)+sum(y>u))*2)
}
emp.extcoef1 <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x]; empirical.extcoef(data.pareto.mat[,x])},mc.cores=5,mc.set.seed = FALSE))
emp.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=emp.extcoef1,symmetric = TRUE,dimnames=NULL)

data.fit <- data[data.max>quantile(data.max,0.75)]
len.row <- unlist(lapply(1:length(data.fit),function(i){length(data.fit[[i]][[1]])}))
data.pareto.mat <- sparseMatrix(i=rep(1:length(data.fit),times=len.row),j=unlist(lapply(1:length(data.fit),function(i){data.fit[[i]][[1]]})),x=unlist(lapply(1:length(data.fit),function(i){data.fit[[i]][[2]]})),dimnames=NULL,symmetric = FALSE)
data.pareto.mat <- as.matrix(data.pareto.mat)
emp.extcoef2 <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x]; empirical.extcoef(data.pareto.mat[,x])},mc.cores=5,mc.set.seed = FALSE))
emp.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=emp.extcoef2,symmetric = TRUE,dimnames=NULL)

save(emp.extcoef1,emp.extcoef2,file="data/application_RedSea/application_RedSea_results_emp.RData")

load("data/extcoef_application_RedSea_MaxStable.RData",e<-new.env())
# load("data/application_RedSea2/application_RedSea_results_1_L-BFGS-B.RData",e1<-new.env())
# load("data/application_RedSea2/application_RedSea_results_2_L-BFGS-B.RData",e2<-new.env())
# load("data/application_RedSea2/application_RedSea_results_3_L-BFGS-B.RData",e3<-new.env())
# load("data/application_RedSea2/application_RedSea_results_4_L-BFGS-B.RData",e4<-new.env())
load("data/fitted_extcoef_RedSea.RData")
# emp.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=apply(pairs,2,empirical_extcoef,data=maxima.frechet),symmetric = TRUE,dimnames=NULL)
emp.extcoef.mat <- e$empirical.extcoef.mat #sparseMatrix(i=pairs[1,],j=pairs[2,],x=e$emp.extcoef1,symmetric = TRUE,dimnames=NULL)
# emp.extcoef.mat2 <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=e$emp.extcoef2,symmetric = TRUE,dimnames=NULL)


p1 <- p2 <- list()
brks.emp1 <- c(1.35,1.5,1.65,1.8,1.9,1.95)#round(quantile(2-emp.extcoef.mat1@x,probs=c(0.005,0.05,0.1,0.3,0.5,0.8,0.9),na.rm=TRUE),4)
brks.emp2 <- c(1.35,1.5,1.65,1.8,1.9,1.95)#round(quantile(2-emp.extcoef.mat2@x,probs=c(0.005,0.05,0.1,0.3,0.5,0.8,0.9),na.rm=TRUE),4)
brks.emp <- c(1.2,1.5,1.7,1.75,1.8,1.85,1.9,1.95)
ewbreaks <- seq(34.4,41.6,2.4)
nsbreaks <- seq(15.6,26.4,3.6)
ewlabels <- unlist(lapply(ewbreaks, function(x) paste(" ",abs(x), "ºE")))
nslabels <- unlist(lapply(nsbreaks, function(x) paste(" ",x, "ºN")))
basis.centers.geo <- sample(1:1043,100)
basis.centers.geo <- basis.centers.geo[order(loc.sub[basis.centers.geo,2])] 
loc.df = data.frame(lon = loc.sub[idx.centers,1],lat = loc.sub[idx.centers,2])
for(i in 1:length(basis.centers.geo)){
    idx.center <- basis.centers.geo[i]
    data.df <- data.frame(lon=round(loc.sub[,1],5),lat=round(loc.sub[,2],5),
                emp=emp.extcoef.mat[,idx.center],br=fitted.extcoef.mat.list[[1]][,idx.center],
                sbr=fitted.extcoef.mat.list[[2]][,idx.center],br2=fitted.extcoef.mat.list[[3]][,idx.center],br3=fitted.extcoef.mat.list[[4]][,idx.center])
    data.df[idx.center,-c(1:2)] = NA
    p1[[i]]<-ggmap(map) +
    geom_tile(data=data.df,aes(x=lon,y=lat,fill=emp),alpha=0.8) + 
    colorspace::scale_fill_continuous_divergingx("RdYlBu",limits=c(1,2),mid=exp(1.6),alpha=0.8,name=expression(hat(theta)[2]),trans="exp") +
    coord_fixed() + stat_contour(data=data.df,aes(x=lon,y=lat,z=sbr),breaks = brks.emp,colour = "black",linetype="dashed") +
    stat_contour(data=data.df,aes(x=lon,y=lat,z=br),breaks = brks.emp,colour = "black",linetype="solid") + labs(x="Longitude", y="Latitude") + 
    #stat_contour(data=data.df,aes(x=lon,y=lat,z=br2),breaks = brks.emp1,colour = "black",linetype="dotted") +
    scale_x_continuous(breaks = ewbreaks, labels = ewlabels,expand=c(0,0),limits = xlim) + 
    scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = ylim) +
    geom_point(data=loc.df,aes(x=lon,y=lat),size=2,fill="black") +
    theme(axis.text = element_text(size=10), 
                            strip.text = element_text(size = 14),
                            axis.title.x = element_text(size=14), 
                            axis.title.y = element_text(size=14),
                            axis.text.y = element_text(angle=90,vjust=0.5,hjust=1),
                            legend.title = element_text(size=14),legend.position = "right",
                            plot.margin=unit(c(0.2,0.2,0,0.5),"cm")) 
    
    p2[[i]]<-ggmap(map) + 
    geom_tile(data=data.df,aes(x=lon,y=lat,fill=emp),alpha=0.8) + 
    colorspace::scale_fill_continuous_divergingx("RdYlBu",limits=c(1,2),alpha=0.8,mid=exp(1.6),name=expression(hat(theta)[2]),trans="exp") + 
    coord_fixed() + stat_contour(data=data.df,aes(x=lon,y=lat,z=br),breaks = brks.emp,colour = "black",linetype="solid") + 
    stat_contour(data=data.df,aes(x=lon,y=lat,z=br2),breaks = brks.emp,colour = "black",linetype="dotted") + labs(x="Longitude",y="Latitude") + 
    stat_contour(data=data.df,aes(x=lon,y=lat,z=br3),breaks = brks.emp,colour = "black",linetype="dashed") +
    scale_x_continuous(breaks = ewbreaks, labels = ewlabels,expand=c(0,0),limits=xlim) + 
    scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = ylim) +
    theme(axis.text = element_text(size=10), 
                            strip.text = element_text(size = 14),
                            axis.title.x = element_text(size=14), 
                            axis.text.y = element_text(angle=90,vjust=0.5,hjust=1),
                            axis.title.y = element_text(size=14), 
                            legend.title = element_text(size=14),legend.position = "right",plot.margin=unit(c(0.2,0.2,0,0.5),"cm"))
}

for(i in 1:length(p1)){
    ggsave(paste0("figures/application/RedSea2/RedSea_extcoef2_1_",sprintf(i,fmt="%.3d"),".png"),p1[[i]],width=5,height=6,dpi=300)
}

for(i in 1:length(p2)){
    ggsave(paste0("figures/application/RedSea2/RedSea_extcoef2_2_",sprintf(i,fmt="%.3d"),".png"),p2[[i]],width=5,height=6,dpi=300)
}

system("magick -delay 20 -loop 0 figures/application/RedSea2/RedSea_extcoef2_1_*.png figures/application/RedSea2/RedSea_combined1_1.gif")

system("magick -delay 20 -loop 0 figures/application/RedSea2/RedSea_extcoef2_2_*.png figures/application/RedSea2/RedSea_combined2_1.gif")

x = distmat[t(pairs)]
png("figures/application/RedSea2/RedSea_extcoef_scatter.png", width=4*3, height=4*2, units="in",res=300)
par(mfrow=c(2,3), mar=c(4,4,2,1),mgp=c(2,1,0))
y = e$empirical.extcoef.mat[t(pairs)]
plot(x,y , pch=20, cex=0.01, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency

plot(x, y, pch=20, cex=0.01, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency
points(x, fitted.extcoef.mat.list[[1]][t(pairs)], pch=20, cex=0.01, col=rgb(0, 0, 1, 0.3))  # Blue with transparency

plot(x, y, pch=20, cex=0.01, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency
points(x, fitted.extcoef.mat.list[[2]][t(pairs)], pch=20, cex=0.01, col=rgb(1, 0, 0, 0.3))  # Red with transparency


plot(x, y, pch=20, cex=0.01, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency

plot(x, y, pch=20, cex=0.01, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency
points(x, fitted.extcoef.mat.list[[3]][t(pairs)], pch=20, cex=0.01, col=rgb(0, 0, 1, 0.3))  # Blue with transparency

plot(x, y, pch=20, cex=0.01, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency
points(x, fitted.extcoef.mat.list[[4]][t(pairs)], pch=20, cex=0.01, col=rgb(1, 0, 0, 0.3))  # Red with transparency
dev.off()

angles.pairs <- (apply(pairs,2,function(x){x1=loc.sub[x[1],1]-loc.sub[x[2],1];y1=loc.sub[x[1],2]-loc.sub[x[2],2];a=atan2(x1,y1)*180/pi}) + 180) %% 180

png("figures/application/RedSea2/RedSea_extcoef_angle/%02d.png", width=4*3, height=4*2,units="in",res=300)
for( angle in seq(0,180,10)){
    par(mfrow=c(2,3), mar=c(4,4,2,1),mgp=c(2,1,0))
    idx = angles.pairs < angle+2 & angles.pairs > angle-2
    plot(x[idx], y[idx], pch=20, cex=0.1, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency

    plot(x[idx], y[idx], pch=20, cex=0.1, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency
    points(x[idx], fitted.extcoef.mat.list[[1]][t(pairs)][idx], pch=20, cex=1, col=rgb(0, 0, 1, 0.5))  # Blue with transparency

    plot(x[idx], y[idx], pch=20, cex=0.1, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency
    points(x[idx], fitted.extcoef.mat.list[[2]][t(pairs)][idx], pch=20, cex=1, col=rgb(1, 0, 0, 0.5))  # Red with transparency


    plot(x[idx], y[idx], pch=20, cex=0.1, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency

    plot(x[idx], y[idx], pch=20, cex=0.1, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency
    points(x[idx], fitted.extcoef.mat.list[[3]][t(pairs)][idx], pch=20, cex=1, col=rgb(0, 0, 1, 0.5))  # Blue with transparency

    plot(x[idx], y[idx], pch=20, cex=0.1, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 2306)) # Black with transparency
    points(x[idx], fitted.extcoef.mat.list[[4]][t(pairs)][idx], pch=20, cex=1, col=rgb(1, 0, 0, 0.5))  # Red with transparency
}
dev.off()



system("cp figures/application/RedSea2/RedSea_extcoef2_[1-2]_0{35,82}.png figures/application/png")




