thres.list = c(80,85,90,95,99)
files.list <- list()
files.list1000 <- list()
for(i in 1:length(thres.list)){
    files.list[[i]] <- list.files(path = "data", pattern = paste0("simulation_study_\\d+_",thres.list[[i]],".RData"), full.names = TRUE, recursive = FALSE)
    files.list1000[[i]] <- list.files(path = "data", pattern = paste0("simulation_study_\\d+_",thres.list[[i]],"_1000.RData"), full.names = TRUE, recursive = FALSE)
}

extract_results <- function(files){
    fit.results <- list()
    for(i in 1:length(files)){
        load(files[[i]],e<-new.env())
        fit.results[[i]] <- e$fit.logskew.angular
    }
    est.mat.list <- list()
    for(i in 1:27){
        est.mat.list[[i]] <- matrix(unlist(lapply(fit.results,function(x){x[[i]]$par})),ncol=5,byrow=TRUE)
    }
    return(est.mat.list)
}



est.mat.list <- lapply(files.list,extract_results)
est.mat.list1000 <- lapply(files.list1000,extract_results)

save(est.mat.list,est.mat.list1000,files.list,files.list1000,file="data/simulation_study_logskew_results.RData")
