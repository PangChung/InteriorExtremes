# LaTeX rows of Table 2 (Florida, dry hours missing, u_p = 95% quantile of r_p) from the result files
load("data/application_florida_5km_missing_spectral.RData")
has.sm <- file.exists("data/application_florida_5km_missing_score.RData")
if(has.sm) load("data/application_florida_5km_missing_score.RData")
f <- function(x,d) formatC(x,format="f",digits=d)
row <- function(model,p,n,par,se,se.hess,aic,k){
    est = c(f(par[1]*60,2),f(par[2],3),f(par[3],2),f(par[4],3),if(k==8) f(par[5:8],2) else rep("-",4))
    ses = c(f(se[1]*60,2),f(se[2],3),f(se[3],2),f(se[4],3),if(k==8) f(se[5:8],2) else rep("-",4))
    # the sandwich needs the inverse Hessian: not reported if -H is not positive definite (zero Hessian SE)
    if(!is.null(se.hess) && any(se.hess[1:k] <= 0)) ses[ses!="-"] = "--"
    ses = ifelse(ses=="-","-",paste0("(",ses,")"))
    c(paste(c(model,p,n,est,if(is.na(aic)) "-" else f(aic,1)),collapse=" & "),
      paste(c(" & & ",ses,""),collapse=" & "))
}
for(nm in names(res.spectral)){
    r = res.spectral[[nm]]; p = sets[[nm]]$p; pl = if(is.infinite(p)) "$\\infty$" else p
    out = c()
    # score-matching fits whose optimisation did not converge are not reported
    if(has.sm && res.score[[nm]]$fit$convergence==0){ s = res.score[[nm]]; out = c(out,row("$\\textrm{BR}^{sm}$",pl,s$n,s$fit$par,s$fit$se,NULL,NA,4)) }
    out = c(out,row("BR",pl,r$n,r$br$par,r$br$se,r$br$se.hess,r$br$aic,4),row("sBR",pl,r$n,r$sbr$par,r$sbr$se,r$sbr$se.hess,r$sbr$aic,8))
    cat(paste0(out," \\\\"),sep="\n")
}
cat("\n% total CPU minutes (BR^sm / BR / sBR) and convergence codes\n")
for(nm in names(res.spectral)){
    r = res.spectral[[nm]]
    cat("%",nm,":",if(has.sm) f(res.score[[nm]]$fit$time,1),"/",f(r$br$time,1),"/",f(r$sbr$time,1),
        "  conv BR/sBR",r$br$convergence,r$sbr$convergence,"  BR start values",paste(signif(r$br.values,8),collapse=", "),"\n")
}
