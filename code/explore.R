args <- commandArgs(TRUE)
computer = "ws"
id = 1
d <- 15 ## 10 * 10 grid on [0,1]^2
m <- 1000 ## number of samples
model = "logskew"; # "logskew" or "truncT"
#model = "truncT"; # "logskew" or "truncT"
for (arg in args) eval(parse(text = arg))
switch(computer,
    "ws" = {DataPath<-"~/Desktop/InteriorExtremes/"},
    "hpc" = {DataPath<-"/srv/scratch/z3536974/";.libPaths("../src")},
    "local" = {DataPath<-"~/Documents/Github/InteriorExtremes/"}
)
# settings 
coord = as.matrix(expand.grid(0:(d-1),0:(d-1))/d)
diff.vector <- cbind(as.vector(outer(coord[,1],coord[,1],'-')),as.vector(outer(coord[,2],coord[,2],'-'))) 
diff.mat <- matrix(apply(diff.vector, 1, function(x) sqrt(sum(x^2))), ncol=nrow(coord))
para.range = c(1,2) #c(0.5,1,2) ## range for the correlation function ##      
para.nu = 1 #c(0.5,1,1.5) ## smoothness parameter for the correlation function ##
para.alpha = rbind(c(0,0,0),c(-1,-2,-3),c(-2,-1,4),c(2,1,4)) ## slant parameter for skewed norm model ##
para.deg = c(2,3) ## degree of the freedom for the truncated t model ##
all.pairs = combn(1:nrow(coord),2)
all.pairs.list = split(all.pairs,col(all.pairs))
thres = c(0.95,0.9)
# loading library and setting path
library(parallel)
library(mvtnorm)
library(TruncatedNormal)
library(evd)
#library(gridExtra)
library(partitions)
#library(ggplot2)
library(Rfast)
library(matrixStats)
library(splines)
source("code/simulation.R")
source("code/exponent_functions.R")
source("code/likelihood_inference.R")
ncores=detectCores()
file2save = paste0(DataPath,"data/simulation_study_",model,"_",id,"_",m,".RData")
init.seed = as.integer((as.integer(Sys.time())/id + sample.int(10^5,1))%%10^5)
set.seed(init.seed)

