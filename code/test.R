rm(list=ls())
.libPaths("../src")
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
library(partitions)
library(Rfast)
library(evd)
library(matrixStats)
library(Matrix)
# load the data##
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
load("data/application_florida_list.RData")
pairs <- comb_n(1:nrow(coord.grid),2)
data.pareto <- mclapply(data,function(x){list(x[[1]],qgpd(x[[2]],1,1,1))},mc.cores=4)
len.row <- unlist(lapply(1:length(data.pareto),function(i){length(data.pareto[[i]][[1]])}))
data.pareto.mat <- sparseMatrix(i=rep(1:length(data.pareto),times=len.row),j=unlist(lapply(1:length(data.pareto),function(i){data.pareto[[i]][[1]]})),x=unlist(lapply(1:length(data.pareto),function(i){data.pareto[[i]][[2]]})),dimnames=NULL,symmetric = FALSE)


empirical.extcoef <- function(data){
    x = data[,1]
    y = data[,2]
    u = quantile(c(x,y),0.9)
    return( sum(x>u & y>u)/(sum(x>u)+sum(y>u))*2)
}

data.pareto.mat.nonsparse <- as.matrix(data.pareto.mat)

system.time({emp.extcoef <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x]; empirical.extcoef(data.pareto.mat.nonsparse[,x])},mc.cores=detectCores()/2,mc.set.seed = FALSE))})

save(data.pareto.mat,emp.extcoef,file="data/application_florida_results_ext_1.RData")
