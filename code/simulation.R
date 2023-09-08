###################################################################################
######## This file contains functions to simulate the max-stable processes ########  
###################################################################################

## Simulate a truncated extremal-t max-stable process
simu <- function(loc,par,model,parallel=TRUE,ncores){    
    nu = par[[1]];sigma=par[[2]]
    n = nrow(sigma)
    x = tmvtnorm::rtmvnorm(1,sigma=sigma,lower=rep(0,n),algorithm = "rejection")
    

}


## Simulate a log-skew normal based max-stable process
simu <- function(loc,par,model,parallel=TRUE,ncores){  
    alpha = par[[1]]; sigma = par[[2]]
    n = nrow(sigma)
    
}
