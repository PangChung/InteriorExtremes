rm(list=ls())
library(sf)
library(ggplot2)
library(ggmap)
library(isoband)
library(Matrix)
## maps of the bivariate extremal coefficients for the Florida rainfall data at the 3x3 block (~5.7km) hourly resolution

load("data/maps_florida.RData") # background map
load("data/application_florida_5km.RData")
load("data/application_florida_5km_diag.RData")
load("data/application_florida_5km_spectral.RData")

## recover the block indices in the rotated radar grid (same construction as in application_florida_5km.R)
load("data/application_florida.RData")
ids <- as.numeric(sub("P","",names(tmp)[-1])); rm(tmp); invisible(gc())
grid.sf <- read_sf("data/Florida/DOPGrid/DOPGrid.shp")
xy <- suppressWarnings(st_coordinates(st_centroid(grid.sf)))[match(ids,grid.sf$PIXEL),1:2]
phi <- atan2(14,1858); R <- matrix(c(cos(phi),sin(phi),-sin(phi),cos(phi)),2)
xr <- xy %*% t(R); xr.min <- apply(xr,2,min)
col <- round((xr[,1]-xr.min[1])/1888); row <- round((xr[,2]-xr.min[2])/1991)
b <- 3
key <- paste(col%/%b,row%/%b); tb <- table(key); keep.key <- names(tb)[tb >= b*b/2]
stopifnot(length(keep.key)==ncol(Z))
bc <- as.integer(sapply(strsplit(keep.key," "),`[`,1)); br <- as.integer(sapply(strsplit(keep.key," "),`[`,2))

# block-index space (i,j) -> longitude/latitude
idx2lonlat <- function(i,j){
    rot = cbind(xr.min[1]+1888*(b*i+(b-1)/2), xr.min[2]+1991*(b*j+(b-1)/2))
    utm = rot %*% R
    st_coordinates(st_transform(st_as_sf(data.frame(x=utm[,1],y=utm[,2]),coords=c("x","y"),crs=st_crs(grid.sf)),4326))
}
# polygons of the blocks
corners <- rbind(c(-0.5,-0.5),c(0.5,-0.5),c(0.5,0.5),c(-0.5,0.5))
poly.df <- do.call(rbind,lapply(1:length(bc),function(k){
    ll = idx2lonlat(bc[k]+corners[,1],br[k]+corners[,2])
    data.frame(id=k,lon=ll[,1],lat=ll[,2])
}))
# contour lines of a field defined on the blocks
contour.df <- function(z,levels){
    M = matrix(NA,max(bc)+1,max(br)+1); M[cbind(bc+1,br+1)] = z
    iso = isolines(x=0:max(bc),y=0:max(br),z=t(M),levels=levels)
    do.call(rbind,lapply(seq_along(iso),function(l){
        if(length(iso[[l]]$x)==0) return(NULL)
        ll = idx2lonlat(iso[[l]]$x,iso[[l]]$y)
        data.frame(lon=ll[,1],lat=ll[,2],group=paste(l,iso[[l]]$id))
    }))
}

brks <- c(1.2,1.5,1.7,1.75,1.8,1.85,1.9,1.95)
ref.ll <- idx2lonlat(bc[idx.ref],br[idx.ref])
# kernel centres of the sBR model (loc = coord.utm/60000 centred at the median)
centers.utm <- sweep(basis.centers*60000,2,apply(coord.utm,2,median),"+")
centers.ll <- st_coordinates(st_transform(st_as_sf(data.frame(x=centers.utm[,1],y=centers.utm[,2]),coords=c("x","y"),crs=st_crs(grid.sf)),4326))
coord.ratio <- 0.01796407/0.01912047
ewbreaks <- c(-82.9,-82.5,-82.1,-81.6); nsbreaks <- c(27.7,28,28.4,28.8)
ewlabels <- unlist(lapply(-ewbreaks,function(x) paste(" ",abs(x),"ºW")))
nslabels <- unlist(lapply(nsbreaks,function(x) paste(" ",x,"ºN")))

plot.map <- function(r,z.solid,z.dashed,centers=FALSE){
    emp = emp.extcoef[,r]; emp[idx.ref[r]] = NA
    fill.df = poly.df; fill.df$emp = emp[fill.df$id]
    ggmap(map) +
    geom_polygon(data=fill.df,aes(x=lon,y=lat,group=id,fill=emp),alpha=0.8,colour=NA) +
    colorspace::scale_fill_continuous_divergingx("RdYlBu",limits=c(1.2,2),mid=exp(1.6),alpha=0.8,name=expression(hat(theta)[2]),trans="exp",na.value="grey40") +
    geom_path(data=contour.df(z.solid,brks),aes(x=lon,y=lat,group=group),colour="black",linetype="solid") +
    geom_path(data=contour.df(z.dashed,brks),aes(x=lon,y=lat,group=group),colour="black",linetype="dashed") +
    geom_point(data=data.frame(lon=ref.ll[r,1],lat=ref.ll[r,2]),aes(x=lon,y=lat),shape=4,size=3,stroke=1.2) +
    {if(centers) geom_point(data=data.frame(lon=centers.ll[,1],lat=centers.ll[,2]),aes(x=lon,y=lat),size=2,colour="black")} +
    # zoom with the coordinate system so that polygons crossing the frame are clipped rather than dropped
    coord_fixed(ratio=1/coord.ratio,xlim=c(-82.9,-81.6),ylim=c(27.5,28.9),expand=FALSE) + labs(x="Longitude",y="Latitude") +
    scale_x_continuous(breaks=ewbreaks,labels=ewlabels) +
    scale_y_continuous(breaks=nsbreaks,labels=nslabels) +
    theme(axis.text = element_text(size=10),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        axis.text.y = element_text(angle=90,vjust=0.5,hjust=1),
        legend.title = element_text(size=14),legend.position = "right",
        plot.margin=unit(c(0.2,0.2,0,0.5),"cm"))
}

dir.out <- "InteriorExtremes_Manuscript/JCGS/figures/florida/"
for(nm in names(extcoef)){
    tag = gsub("[^0-9a-zA-Z]+","_",nm)
    for(r in 1:2){
        # BR (solid) vs sBR (dashed), both fitted with the spectral likelihood
        ggsave(paste0(dir.out,"florida5km_extcoef_sbr_",tag,"_ref",r,".png"),plot.map(r,extcoef[[nm]]$BR[,r],extcoef[[nm]]$sBR[,r],centers=TRUE),width=6.4,height=6,dpi=300)
        # BR fitted with the spectral likelihood (solid) vs score matching (dashed)
        ggsave(paste0(dir.out,"florida5km_extcoef_sm_",tag,"_ref",r,".png"),plot.map(r,extcoef[[nm]]$BR[,r],extcoef[[nm]]$BR.sm[,r]),width=6.4,height=6,dpi=300)
    }
}
