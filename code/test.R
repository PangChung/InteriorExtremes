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
data.pareto.mat <- sparseMatrix(i=rep(1:length(data.pareto),times=len.row),j=unlist(lapply(1:length(data.pareto),function(i){data.pareto[[i]][[1]]})),x=unlist(lapply(1:length(data.pareto),function(i){data.pareto[[i]][[2]]})),dimnames=NULL,symmetric = FALSE)


empirical.extcoef <- function(data){
    u = 10
    x = data[,1]
    y = data[,2]
    return( sum(x>u & y>u)/(sum(x>u)+sum(y>u))*2)
}

data.pareto.mat.nonsparse <- as.matrix(data.pareto.mat)

system.time({emp.extcoef <- unlist(mclapply(1:ncol(pairs),function(x){x=pairs[,x]; empirical.extcoef(data.pareto.mat.nonsparse[,x])},mc.cores=detectCores(),mc.set.seed = FALSE))})

save(data.pareto.mat,emp.extcoef,file="data/application_florida/application_florida_results_ext_1.RData")
