files <- list.files(path="/srv/scratch/z3536974/data/",pattern=".RData",full.name=TRUE)
count = length(files)
for(f in files){
    load(f,envir=e<-new.env())
    save(fit.logskew.angular,envir=e,file=f)
    print(count <- count - 1)
}
