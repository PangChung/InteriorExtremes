# two-step selection: r_p exceedances of a high quantile u_p, then r_0 > u_p (Theorem 3)
suppressMessages(library(matrixStats))
source(textConnection(paste(readLines("code/florida_agg_fit.R")[2:28],collapse="\n")))
lp.norm <- function(Z,p){ m = rowMaxs(Z); if(is.infinite(p)) m else m*rowSums((Z/m)^p)^(1/p) }
out <- list()
for(cf in list(c(3,1),c(4,3),c(6,3),c(4,6),c(6,6),c(8,6),c(8,24))){
    b <- cf[1]; h <- cf[2]
    load(sprintf("data/florida_agg/hourly_b%d.RData",b))
    Yh <- if(h==1) Y else rowsum(Y,rep(seq_len(nrow(Y)/h),each=h),reorder=FALSE)
    wet <- Yh > 0; Z <- apply(Yh,2,to.pareto); Z[!wet] <- 0; D <- ncol(Z); Nt <- nrow(Z)
    r0 <- rowSums(Z)/D
    for(p in c(3,10,Inf)){
        rp <- lp.norm(Z,p)/(if(is.infinite(p)) 1 else D^(1/p)); rp[!is.finite(rp)] <- 0
        for(pe in c(0.05,0.02,0.01,0.005)){
            u <- quantile(rp,1-pe,names=FALSE,type=1); s1 <- rp > u; s2 <- r0 > u
            idx <- which(s2)
            out[[length(out)+1]] <- data.frame(agg=sprintf("%dx%d,%dh",b,b,h),p=p,top.pct=100*pe,u=round(u,1),n.p=sum(s1),n=sum(s2),
                ret=round(mean(s2)/mean(s1),3),ncl=if(length(idx)) max(cumsum(c(1,diff(idx)>24/h))) else 0,
                miss=if(length(idx)) round(1-mean(wet[idx,]),3) else NA,
                mlog.rp=round(mean(log(rp[s1]/u)),2),mlog.r0=if(length(idx)) round(mean(log(r0[s2]/u)),2) else NA)
        }
    }
}
res <- do.call(rbind,out); save(res,file="data/florida_agg/twostep.RData")
options(width=200); print(res,row.names=FALSE)
