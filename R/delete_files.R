realm = "Mar"

# sp list
SDMsSpList <- list.files(sdm_dir(realm))

# delete files
for(i in 1:length(SDMsSpList)){ cat("\rsps",i,"from",length(SDMsSpList))
    myfiles <- list.files(here::here(sdm_dir(realm),SDMsSpList[i],gsub("_",".",SDMsSpList[i])),full.names = T)
    test <- strsplit(myfiles,"/")
    test <- sapply(test, function(x){
        tmp = x[8]
        any(duplicated(strsplit(tmp," ")[[1]]))
    })
    if(any(test)){
        unlink(myfiles[which(test)],recursive = TRUE)
    }
}

##############################
# sp list
existing <- list.files(here::here("/media/seagate/boliveira/SDMs/MaxNet",realm))
missing_ones <- existing[-which(existing %in% SDMsSpList)]

# delete files
for(i in 1:length(missing_ones)){
    mydir <- here::here("/media/seagate/boliveira/SDMs/MaxNet",realm,missing_ones[i])
    unlink(mydir,recursive = TRUE)
}

###############################
# sp list
SDMsSpList <- list.files(here::here("/media/seagate/boliveira/SDMs/MaxNet",realm))

# delete files
for(i in 1:length(SDMsSpList)){
    myfiles <- list.files(here::here("/media/seagate/boliveira/SDMs/MaxNet",realm,SDMsSpList[i]),
                          full.names = T, pattern = ".tif")
    myfiles <- myfiles[grep(" bios",myfiles)]
    rem <- sapply(myfiles, function(x){
        tmp <- rast(x)
        any(names(tmp)=="mean")
    })
    
    if(any(rem)){
        myfiles <- myfiles[which(rem)]
        unlink(myfiles)
    }
}

###############################
# move files

files2move <- list.files(here::here("Exposure-SDM/Data"), pattern = "_P")

for(i in 1:length(files2move)){ cat("\r",i, "from", length(files2move))
    system(paste("mv", here::here("Exposure-SDM/Data",files2move[i]), here::here("/scratch/boliveira/tmp",files2move[i])))
}
