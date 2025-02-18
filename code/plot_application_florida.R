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
load("data/application_florida_list.RData")
load("data/maps_florida.RData")

grid.sf <- read_sf("data/Florida/DOPGrid/DOPGrid.shp")
intb.sf <- read_sf("data/Florida/INTB_Basins/INTB_Basins.shp")
fill_basin <- as.factor(rep(1:8,length.out=172))
all_IDs <- names(read.csv("data/Florida/PixelRain15min_1995.csv", header = TRUE, nrows = 1))[-1]
all_IDs_num <- as.numeric(stringi::stri_extract_first(all_IDs, regex = "[0-9]+"))
fill_grid = grid.sf$PIXEL %in% all_IDs_num


# register_google(key=system("echo $g_key",intern = TRUE))
# map <- get_googlemap(center=c(-82.273528,28.209394),zoom=9,maptype = "terrain",style = "feature:all|element:all|saturation:-100|lightness:50")
# intb.transform <- st_transform(intb.sf, crs = 4326)
# grid.transform <- st_transform(grid.sf, crs = 4326)

# save(intb.transform,grid.transform,map,grid.sf,intb.sf,file="data/maps_florida.RData")

# ggmap(map) + theme_void() + 
# ggtitle("Tampa Bay") + theme(plot.title = element_text(hjust = 0.5))  + 
# geom_sf(data=grid.transform, aes(colour=as.factor(fill_grid)),alpha=0.5,inherit.aes = FALSE) + scale_color_brewer(palette = "Set1",name="Grid")

coord.geo <- as.data.frame(st_coordinates(grid.transform))
idx.pixel <- unlist(lapply(all_IDs_num,function(x){which(grid.transform$PIXEL==x)}))
coord.geo <- matrix(unlist(lapply(idx.pixel, function(i){apply(coord.geo[coord.geo$L2==i,1:2][1:4,],2,mean)})),ncol=2,byrow=TRUE)
coord.grid = coord.grid/60000
coord.grid.new <- apply(coord.grid,2,function(x){x-median(x)})
pairs <- comb_n(1:nrow(coord.grid),2)
coord.ratio <- 0.01796407/0.01912047
basis.centers <- as.matrix(expand.grid(quantile(coord.grid[,1],c(0.2,0.8)),quantile(coord.grid[,2],c(0.2,0.8))))
# basis.centers <- rbind(basis.centers,apply(basis.centers,2,function(x){median(x)}))
basis <- lapply(1:nrow(basis.centers),function(i){
    y=dnorm(sqrt((coord.grid[,1]-basis.centers[i,1])^2 + (coord.grid[,2]-basis.centers[i,2])^2),mean=0,sd=ncol(coord.grid)*2)
    y=y-mean(y)
    y/sqrt(sum(y^2))
})
# basis.centers <- expand.grid(quantile(coord.grid[,1],c(0.2,0.8)),quantile(coord.grid[,2],c(0.2,0.8)))
# basis <- lapply(1:nrow(basis.centers),function(i){
#     y=dnorm(sqrt((coord.grid[,1]-basis.centers[i,1])^2 + (coord.grid[,2]-basis.centers[i,2])^2),mean=0,sd=ncol(coord.grid)*2)
#     y=y-mean(y)
#     y/sqrt(sum(y^2))
# })
basis <- matrix(unlist(basis),nrow=nrow(coord.grid),byrow=FALSE)

## fitted extcoef ##
load("data/application_results.RData")
unit = 2/0.03128403
fit.results <- subset(data.florida,type=="full")
fit.results$lambda = fit.results$lambda/unit
fitted.extcoef.mat.list <- list()

for(i in 1:6){
    alpha <- alpha.func(as.numeric(fit.results[i,5:8]),basis)
    cov.mat <- vario.func2(coord.grid.new,as.numeric(fit.results[i,1:4]),ncores=5)
    par.list <- alpha2delta(list(cov.mat,alpha))
    fitted.extcoef <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x];V_bi_logskew(c(1,1),list(par.list[[1]][x,x],par.list[[2]][x]),alpha.para=FALSE)},mc.cores=5,mc.set.seed = FALSE))
    range(fitted.extcoef)
    fitted.extcoef.mat.list[[i]] <- sparseMatrix(i=pairs[2,],j=pairs[1,],x=fitted.extcoef,symmetric = TRUE,dimnames=NULL) 
}
save(fitted.extcoef.mat.list,file="data/fitted_extcoef_florida.RData")

# empirical extcoef
data.pareto <- mclapply(data,function(x){list(x[[1]],qgpd(x[[2]],1,1,1))},mc.cores=4)

data.sum <- unlist(mclapply(data.pareto,function(x){mean(x[[2]])},mc.cores=4))

data.max <- unlist(mclapply(data.pareto,function(x){max(x[[2]])},mc.cores=4))

data.fit <- data.pareto[data.sum>quantile(data.sum,0.95)]

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


data.fit <- data.pareto[data.max>quantile(data.max,0.95)]
len.row <- unlist(lapply(1:length(data.fit),function(i){length(data.fit[[i]][[1]])}))
data.pareto.mat <- sparseMatrix(i=rep(1:length(data.fit),times=len.row),j=unlist(lapply(1:length(data.fit),function(i){data.fit[[i]][[1]]})),x=unlist(lapply(1:length(data.fit),function(i){data.fit[[i]][[2]]})),dimnames=NULL,symmetric = FALSE)
data.pareto.mat <- as.matrix(data.pareto.mat)

emp.extcoef2 <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x]; empirical.extcoef(data.pareto.mat[,x])},mc.cores=5,mc.set.seed = FALSE))

emp.extcoef.mat <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=emp.extcoef2,symmetric = TRUE,dimnames=NULL)

save(emp.extcoef1,emp.extcoef2,file="data/application_florida/application_florida_results_emp.RData")


load("data/application_florida/application_florida_results_emp.RData",e<-new.env())
# load("data/application_florida/application_florida_results_1_L-BFGS-B.RData",e1<-new.env())
# load("data/application_florida/application_florida_results_2_L-BFGS-B.RData",e2<-new.env())
# load("data/application_florida/application_florida_results_3_L-BFGS-B.RData",e3<-new.env())
# load("data/application_florida/application_florida_results_4_L-BFGS-B.RData",e4<-new.env())
load("data/fitted_extcoef_florida.RData")

emp.extcoef.mat1 <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=e$emp.extcoef1,symmetric = TRUE,dimnames=NULL)
emp.extcoef.mat2 <- sparseMatrix(i=pairs[1,],j=pairs[2,],x=e$emp.extcoef2,symmetric = TRUE,dimnames=NULL)
basis.centers <- as.matrix(expand.grid(quantile(coord.geo[,1],c(0.2,0.8)),quantile(coord.geo[,2],c(0.2,0.8))))
loc.df <- data.frame(lon=basis.centers[,1],lat=basis.centers[,2])
p1 <- p2 <- p3 <- p4 <- list()
brks.emp <- c(1.2,1.5,1.7,1.75,1.8,1.85,1.9,1.95)#round(quantile(2-emp.extcoef.mat1@x,probs=c(0.005,0.05,0.1,0.3,0.5,0.8,0.9),na.rm=TRUE),4)
ewbreaks <- c(-82.9,-82.5,-82.1,-81.6)
nsbreaks <- c(27.7, 28, 28.4, 28.8)
ewlabels <- unlist(lapply(-ewbreaks, function(x) paste(" ",abs(x), "ºW")))
nslabels <- unlist(lapply(nsbreaks, function(x) paste(" ",x, "ºN")))
basis.centers.geo <- sample(1:4449,50)
basis.centers.geo <- basis.centers.geo[order(coord.geo[basis.centers.geo,2])] 
for(i in 1:length(basis.centers.geo)){
    idx.center = basis.centers.geo[i]
    data.df <- data.frame(lon=round(coord.geo[,1],5),lat=round(coord.geo[,2],5),
                emp1=2-emp.extcoef.mat1[,idx.center],emp2=2-emp.extcoef.mat2[,idx.center],br=fitted.extcoef.mat.list[[1]][,idx.center],
                sbr=fitted.extcoef.mat.list[[3]][,idx.center],br2=fitted.extcoef.mat.list[[2]][,idx.center],sbr2=fitted.extcoef.mat.list[[4]][,idx.center],br2.1 = fitted.extcoef.mat.list[[5]][,idx.center],br2.2=fitted.extcoef.mat.list[[6]][,idx.center])
    data.df[idx.center,-c(1:2)] = NA
    p1[[i]]<-ggmap(map) +
    geom_tile(data=data.df,aes(x=lon,y=lat,fill=emp1),alpha=0.8) + 
    colorspace::scale_fill_continuous_divergingx("RdYlBu",limits=c(1.2,2),mid=exp(1.6),alpha=0.8,name=expression(hat(theta)[2]),trans="exp") +
    coord_fixed(ratio=1/coord.ratio) + stat_contour(data=data.df,aes(x=lon,y=lat,z=br),breaks = brks.emp,colour = "black",linetype="dashed") +
    stat_contour(data=data.df,aes(x=lon,y=lat,z=sbr),breaks = brks.emp,colour = "black") + labs(x="Longitude", y="Latitude") + 
    scale_x_continuous(breaks = ewbreaks, labels = ewlabels,expand=c(0,0),limits=c(-82.9,-81.6)) + 
    scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(27.5,28.9)) + 
    geom_point(data=loc.df,aes(x=lon,y=lat),size=2,fill="black") +
    theme(axis.text = element_text(size=10), 
                            strip.text = element_text(size = 14),
                            axis.title.x = element_text(size=14), 
                            axis.title.y = element_text(size=14),
                            axis.text.y = element_text(angle=90,vjust=0.5,hjust=1),
                            legend.title = element_text(size=14),legend.position = "right",
                            plot.margin=unit(c(0.2,0.2,0,0.5),"cm")) 
    
    p2[[i]]<-ggmap(map) + 
    geom_tile(data=data.df,aes(x=lon,y=lat,fill=emp2),alpha=0.8) + 
    colorspace::scale_fill_continuous_divergingx("RdYlBu",limits=c(1.2,2),alpha=0.8,mid=exp(1.6),name=expression(hat(theta)[2]),trans="exp") + 
    coord_fixed(ratio=1/coord.ratio) + stat_contour(data=data.df,aes(x=lon,y=lat,z=br2),breaks = brks.emp,colour = "black",linetype="solid") + 
    stat_contour(data=data.df,aes(x=lon,y=lat,z=sbr2),breaks = brks.emp,colour = "black",linetype="dashed") + labs(x="Longitude",y="Latitude") + 
    scale_x_continuous(breaks = ewbreaks, labels = ewlabels,expand=c(0,0),limits=c(-82.9,-81.6)) + 
    scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(27.5,28.9)) +
    geom_point(data=loc.df,aes(x=lon,y=lat),size=2,fill="black") +
    theme(axis.text = element_text(size=10), 
                            strip.text = element_text(size = 14),
                            axis.title.x = element_text(size=14), 
                            axis.text.y = element_text(angle=90,vjust=0.5,hjust=1),
                            axis.title.y = element_text(size=14), 
                            legend.title = element_text(size=14),legend.position = "right",plot.margin=unit(c(0.2,0.2,0,0.5),"cm"))

    p3[[i]]<-ggmap(map) +
    geom_tile(data=data.df,aes(x=lon,y=lat,fill=emp1),alpha=0.8) + 
    colorspace::scale_fill_continuous_divergingx("RdYlBu",limits=c(1.2,2),mid=exp(1.6),alpha=0.8,name=expression(hat(theta)[2]),trans="exp") +
    coord_fixed(ratio=1/coord.ratio) + stat_contour(data=data.df,aes(x=lon,y=lat,z=br),breaks = brks.emp,colour = "black",linetype="solid") +
    stat_contour(data=data.df,aes(x=lon,y=lat,z=br2.1),breaks = brks.emp,colour = "black",linetype="dashed") + labs(x="Longitude", y="Latitude") + 
    scale_x_continuous(breaks = ewbreaks, labels = ewlabels,expand=c(0,0),limits=c(-82.9,-81.6)) + 
    scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(27.5,28.9)) + 
    # geom_point(data=loc.df,aes(x=lon,y=lat),size=2,fill="black") +
    theme(axis.text = element_text(size=10), 
                            strip.text = element_text(size = 14),
                            axis.title.x = element_text(size=14), 
                            axis.title.y = element_text(size=14),
                            axis.text.y = element_text(angle=90,vjust=0.5,hjust=1),
                            legend.title = element_text(size=14),legend.position = "right",
                            plot.margin=unit(c(0.2,0.2,0,0.5),"cm")) 
                            
    p4[[i]]<-ggmap(map) +
    geom_tile(data=data.df,aes(x=lon,y=lat,fill=emp2),alpha=0.8) + 
    colorspace::scale_fill_continuous_divergingx("RdYlBu",limits=c(1.2,2),mid=exp(1.6),alpha=0.8,name=expression(hat(theta)[2]),trans="exp") +
    coord_fixed(ratio=1/coord.ratio) + stat_contour(data=data.df,aes(x=lon,y=lat,z=br2),breaks = brks.emp,colour = "black",linetype="solid") +
    stat_contour(data=data.df,aes(x=lon,y=lat,z=br2.2),breaks = brks.emp,colour = "black",linetype="dashed") + labs(x="Longitude", y="Latitude") + 
    scale_x_continuous(breaks = ewbreaks, labels = ewlabels,expand=c(0,0),limits=c(-82.9,-81.6)) + 
    scale_y_continuous(breaks = nsbreaks, labels = nslabels, expand = c(0, 0), limits = c(27.5,28.9)) + 
    # geom_point(data=loc.df,aes(x=lon,y=lat),size=2,fill="black") +
    theme(axis.text = element_text(size=10), 
                            strip.text = element_text(size = 14),
                            axis.title.x = element_text(size=14), 
                            axis.title.y = element_text(size=14),
                            axis.text.y = element_text(angle=90,vjust=0.5,hjust=1),
                            legend.title = element_text(size=14),legend.position = "right",
                            plot.margin=unit(c(0.2,0.2,0,0.5),"cm")) 
}

for(i in 1:length(p1)){
    ggsave(paste0("figures/application/florida/florida_extcoef2_1_",sprintf(i,fmt="%.3d"),".png"),p1[[i]],width=6.4,height=6,dpi=300)
    ggsave(paste0("figures/application/florida/florida_extcoef2_2_",sprintf(i,fmt="%.3d"),".png"),p2[[i]],width=6.4,height=6,dpi=300)
    ggsave(paste0("figures/application/florida/florida_extcoef2_3_",sprintf(i,fmt="%.3d"),".png"),p3[[i]],width=6.4,height=6,dpi=300)
    ggsave(paste0("figures/application/florida/florida_extcoef2_4_",sprintf(i,fmt="%.3d"),".png"),p4[[i]],width=6.4,height=6,dpi=300)
}

system("magick -delay 20 -loop 0 figures/application/florida/florida_extcoef2_1_*.png figures/application/florida/combined1_1.gif")

system("magick -delay 20 -loop 0 figures/application/florida/florida_extcoef2_2_*.png figures/application/florida/combined2_1.gif")


unit = 2/0.03128403
x = as.matrix(dist(coord.grid))[t(pairs)]*unit

png("figures/application/florida/florida_extcoef.png", width=4, height=4.5*2,units="in",res=300)
par(mfrow=c(2,1), mar=c(3,3.5,0,0),mgp=c(2,1,0))

plot(x, 2-e$emp.extcoef1, pch=20, cex=0.01, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 196)) # Black with transparency
points(x, fitted.extcoef.mat.list[[3]][t(pairs)], pch=20, cex=0.02, col=rgb(0, 0, 1, 0.3))  # Blue with transparency
points(x, fitted.extcoef.mat.list[[1]][t(pairs)], pch=20, cex=0.02, col=rgb(1, 0, 0, 0.3))  # Red with transparency

plot(x, 2-e$emp.extcoef2, pch=20, cex=0.01, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 196)) # Black with transparency
points(x, fitted.extcoef.mat.list[[4]][t(pairs)], pch=20, cex=0.02, col=rgb(0, 0, 1, 0.3))  # Blue with transparency
points(x, fitted.extcoef.mat.list[[2]][t(pairs)], pch=20, cex=0.02, col=rgb(1, 0, 0, 0.3))  # Red with transparency
dev.off()

angles.pairs <- (apply(pairs,2,function(x){x1=coord.geo[x[1],1]-coord.geo[x[2],1];y1=coord.geo[x[1],2]-coord.geo[x[2],2];a=atan2(x1,y1)*180/pi}) + 180) %% 180


png("figures/application/florida/florida_extcoef_angle/%02d.png", width=4, height=4.5*2,units="in",res=300)
for( angle in seq(0,180,10)){
    par(mfrow=c(2,1), mar=c(3,3.5,0,0),mgp=c(2,1,0))
    idx = angles.pairs < angle+2 & angles.pairs > angle-2

    plot(x[idx], 2-e$emp.extcoef1[idx], pch=20, cex=0.1, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 196)) # Black with transparency
    points(x[idx], fitted.extcoef.mat.list[[3]][t(pairs)][idx], pch=20, cex=0.2, col=rgb(0, 0, 1, 0.5))  # Blue with transparency
    points(x[idx], fitted.extcoef.mat.list[[1]][t(pairs)][idx], pch=20, cex=0.2, col=rgb(1, 0, 0, 0.5))  # Red with transparency
    

    plot(x[idx], 2-e$emp.extcoef2[idx], pch=20, cex=0.1, xlab="Distance (km)", ylab=expression(hat(theta)[2]),col=rgb(0, 0, 0, 0.5),ylim=c(1,2), xlim=c(0, 196)) # Black with transparency
    points(x[idx], fitted.extcoef.mat.list[[4]][t(pairs)][idx], pch=20, cex=0.2, col=rgb(0, 0, 1, 0.5))  # Blue with transparency
    points(x[idx], fitted.extcoef.mat.list[[2]][t(pairs)][idx], pch=20, cex=0.2, col=rgb(1, 0, 0, 0.5))  # Red with transparency
}
dev.off()

summary(abs(2-e$emp.extcoef1-fitted.extcoef.mat.list[[3]][t(pairs)])) - summary(abs(2-e$emp.extcoef1-fitted.extcoef.mat.list[[3]][t(pairs)]))
summary(abs(2-e$emp.extcoef2-fitted.extcoef.mat.list[[4]][t(pairs)])) - summary(abs(2-e$emp.extcoef1-fitted.extcoef.mat.list[[4]][t(pairs)]))

## collect jackknife results ## 
# est.list <- list()
# est.sd.list <- list()
# for(i in 1:4){
#     est.mat <- matrix(NA,nrow=70,ncol=8)
#     file = paste0("data/application_florida/application_florida_results_",i,"_L-BFGS-B.RData")
#     load(file, e.tmp<-new.env())
#     est <- e.tmp$fit.result$par
#     est[3] = est[3] %% pi/4
#     for(j in 1:70){
#         file = paste0("data/application_florida/application_florida_results_",i,"_L-BFGS-B_",j,".RData")
#         load(file, e.tmp<-new.env())
#         est.mat[j,] <- e.tmp$fit.result$par
#     }
#     est.jack = nrow(est.mat) * est.mat - (nrow(est.mat) - 1) * matrix(est,nrow=nrow(est.mat),ncol=ncol(est.mat),byrow=TRUE)
#     est.jack[,3] <- est.jack[,3] %% pi/4
#     est.sd = apply(est.jack,2,sd)/sqrt(nrow(est.jack))
#     est.list[[i]] <- est
#     est.sd.list[[i]] <- est.sd
# }

# est <- do.call(rbind,est.list)
# est.sd <- do.call(rbind,est.sd.list)
# est[,1] <- est[,1]*unit
# est.sd[,1] <- est.sd[,1]*unit
# round(est*100,2)
# round(est.sd*100,2)

# a = matrix(NA,nrow=4,ncol=8)
# for(i in 1:4){
#     for(j in 1:8){
#         a[i,j] = paste0(round(est[i,j]*100,2)," (",round(est.sd[i,j]*100,2),")")
#     }
# }
# xtable(a) 

system("cp figures/application/florida/florida_extcoef2_[1-4]_0{34,05}.png figures/application/png")
