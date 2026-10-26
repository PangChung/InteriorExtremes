# summary of the missing-data BR spectral fits across aggregations
f <- "data/florida_agg/"
all <- list()
if(file.exists(paste0(f,"fits_missing_BR.RData"))){load(paste0(f,"fits_missing_BR.RData")); all <- c(all,res)}
if(file.exists(paste0(f,"fits_missing_BR_p3p10.RData"))){load(paste0(f,"fits_missing_BR_p3p10.RData")); all <- c(all,res2)}
if(file.exists(paste0(f,"fits_missing_BR_heavy.RData"))){load(paste0(f,"fits_missing_BR_heavy.RData")); all <- c(all,list(res.heavy))}
all <- all[sapply(all,function(r) is.null(r$err))]
tab <- do.call(rbind,lapply(all,function(r) data.frame(
    agg=sprintf("%dx%d,%dh",r$job$b,r$job$b,r$job$h),D=r$D,p=r$job$p,n.p=r$n.p,n=r$n,ncl=r$ncl,miss=round(r$miss,3),
    lambda.km=round(r$par[1]*60,2),se.l=round(r$se$sand[1]*60,2),vartheta=round(r$par[2],2),se.v=round(r$se$sand[2],2),
    angle=round(r$par[3],2),se.a=round(r$se$sand[3],2),aniso=round(r$par[4],2),se.an=round(r$se$sand[4],2),
    sand.hess=round(mean(r$se$sand/r$se$hess),2),
    ret.emp=round(r$ret.emp,3),ret.mod=round(r$ret.model,3),rmse=round(r$rmse,3),rmse50=round(r$rmse50,3),
    min=round(r$time/60,1),conv=r$conv)))
tab <- tab[order(tab$D,-tab$n.p,tab$agg,tab$p),]
options(width=250); print(tab,row.names=FALSE)
cat("\nmean (empirical - fitted) extremal coefficient by distance (km):\n")
bins <- do.call(rbind,lapply(all,function(r) round(r$theta.bins,3)))
rownames(bins) <- sapply(all,function(r) sprintf("%dx%d,%dh,p=%s",r$job$b,r$job$b,r$job$h,r$job$p))
print(bins[order(rownames(bins)),])
