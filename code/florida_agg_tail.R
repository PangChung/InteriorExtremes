# tail checks for r_0 under unit Pareto margins (dry = 0): t P(r_0 > t) -> 1 and E[log(r_0/t) | r_0 > t] -> 1
suppressMessages(library(matrixStats))
source(textConnection(paste(readLines("code/florida_agg_fit.R")[2:28],collapse="\n"))) # pwm.gpd, to.pareto, lp.norm
lp.norm <- function(Z,p){ m = rowMaxs(Z); if(is.infinite(p)) m else m*rowSums((Z/m)^p)^(1/p) }
out <- list()
for(cf in list(c(3,1),c(4,3),c(6,3),c(4,6),c(6,6),c(8,6))){
    b <- cf[1]; h <- cf[2]
    load(sprintf("data/florida_agg/hourly_b%d.RData",b))
    Yh <- if(h==1) Y else rowsum(Y,rep(seq_len(nrow(Y)/h),each=h),reorder=FALSE)
    wet <- Yh > 0; Z <- apply(Yh,2,to.pareto); Z[!wet] <- 0; D <- ncol(Z); Nt <- nrow(Z)
    r0 <- rowSums(Z)/D
    # thresholds u_p used in the fits (86% quantile of r_p), and fixed exceedance proportions of r_0
    up <- sapply(c(3,10,Inf),function(p){rp <- lp.norm(Z,p)/(if(is.infinite(p)) 1 else D^(1/p)); rp[!is.finite(rp)] <- 0; quantile(rp,0.86,names=FALSE,type=1)})
    tt <- c(setNames(up,c("u3","u10","uInf")), setNames(quantile(r0,1-c(0.05,0.02,0.01,0.005),names=FALSE),c("q95","q98","q99","q995")))
    for(k in seq_along(tt)){ t <- tt[k]; ex <- r0 > t
        out[[length(out)+1]] <- data.frame(agg=sprintf("%dx%d,%dh",b,b,h),thr=names(tt)[k],t=round(t,1),n=sum(ex),pct=round(100*mean(ex),2),
            tP=round(t*mean(ex),2),meanlog=round(mean(log(r0[ex]/t)),2),se=round(sd(log(r0[ex]/t))/sqrt(sum(ex)),2)) }
}
options(width=200); print(do.call(rbind,out),row.names=FALSE)
