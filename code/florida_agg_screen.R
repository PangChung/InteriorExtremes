# screening of spatial (b) and temporal (h hours) aggregation with dry sites treated as missing
suppressMessages({library(matrixStats)})
pwm.gpd <- function(y){
    y <- sort(y); n <- length(y); p <- (seq_len(n)-0.35)/n
    a0 <- mean(y); a1 <- mean((1-p)*y)
    c(2-a0/(a0-2*a1), 2*a0*a1/(a0-2*a1))
}
to.pareto <- function(y){
    n <- length(y); Fe <- rank(y,ties.method="average")/(n+1)
    thr <- quantile(y[y>0],0.95,names=FALSE); ex <- y>thr
    if(sum(ex)>=30){
        pr <- pwm.gpd(y[ex]-thr)
        z <- pmax(1+pr[1]*(y[ex]-thr)/pr[2],1e-300)
        Fe[ex] <- 1-pmax(mean(ex)*z^(-1/pr[1]),0.1/(n+1))
    }
    1/(1-Fe)
}
lp.norm <- function(Z,p){ m = rowMaxs(Z); if(is.infinite(p)) m else m*rowSums((Z/m)^p)^(1/p) }
res <- list()
for(b in c(2,3,4,6,8)){
    load(sprintf("data/florida_agg/hourly_b%d.RData",b))
    for(h in c(1,3,6,24)){
        Yh <- if(h==1) Y else rowsum(Y,rep(seq_len(nrow(Y)/h),each=h),reorder=FALSE)
        D <- ncol(Yh); Nt <- nrow(Yh)
        wet <- Yh > 0
        Z <- apply(Yh,2,to.pareto); Z[!wet] <- 0
        r0 <- rowSums(Z)/D
        for(q in c(0.86,0.912)) for(p in c(3,10,Inf)){
            rp <- lp.norm(Z,p)/(if(is.infinite(p)) 1 else D^(1/p)); rp[!is.finite(rp)] <- 0
            u <- quantile(rp,q,names=FALSE,type=1)
            idx.p <- which(rp>u); idx <- which(r0>u)
            cl <- cumsum(c(1,diff(idx) > 24/h))
            res[[length(res)+1]] <- data.frame(b=b,h=h,D=D,Nt=Nt,dry.all=round(1-mean(wet),3),alldry.t=round(mean(rowSums(wet)==0),3),
                q=q,p=p,u=signif(u,3),n.p=length(idx.p),n=length(idx),ret=round(length(idx)/length(idx.p),3),
                miss.sel=round(1-mean(wet[idx,]),3),miss.p=round(1-mean(wet[idx.p,]),3),
                U.med=median(rowSums(wet[idx,,drop=FALSE])),ncl=max(cl))
        }
        cat(b,h,"\n")
    }
}
res <- do.call(rbind,res)
save(res,file="data/florida_agg/screen.RData")
options(width=250); print(res,row.names=FALSE)
