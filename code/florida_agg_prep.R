# block-aggregated hourly Florida radar totals for several block sizes
suppressMessages({library(sf); library(Matrix)})
load("data/application_florida.RData")
ids <- as.numeric(sub("P","",names(tmp)[-1]))
grid.sf <- read_sf("data/Florida/DOPGrid/DOPGrid.shp")
xy <- suppressWarnings(st_coordinates(st_centroid(grid.sf)))[match(ids,grid.sf$PIXEL),1:2]
X <- as.matrix(tmp[,-1]); rm(tmp); invisible(gc())
phi <- atan2(14,1858)
xr <- xy %*% t(matrix(c(cos(phi),sin(phi),-sin(phi),cos(phi)),2))
col <- round((xr[,1]-min(xr[,1]))/1888); row <- round((xr[,2]-min(xr[,2]))/1991)
for(b in c(2,3,4,6,8)){
    key <- paste(col%/%b,row%/%b); tb <- table(key); keep.key <- names(tb)[tb >= b*b/2]
    keep <- key %in% keep.key; f <- factor(key[keep],levels=keep.key); cnt <- as.numeric(table(f))
    M <- sparseMatrix(i=which(keep),j=as.integer(f),x=1/cnt[as.integer(f)],dims=c(ncol(X),length(keep.key)))
    Y <- as.matrix(X %*% M)
    Y <- matrix(colSums(matrix(Y,nrow=4)),ncol=ncol(Y)) # hourly totals
    coord.utm <- as.matrix(t(M) %*% xy)
    save(Y,coord.utm,cnt,file=sprintf("data/florida_agg/hourly_b%d.RData",b))
    cat("b =",b," D =",ncol(Y)," Nt =",nrow(Y),"\n")
}
