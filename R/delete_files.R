
# sp list
SDMsSpList <- list.files(here::here("/media/seagate/boliveira/SDMs/MaxNet",realm))

# delete files
for(i in 1:length(SDMsSpList)){
    myfiles <- list.files(here::here("/media/seagate/boliveira/SDMs/MaxNet",realm,SDMsSpList[i]),full.names = T)
    if(any(grepl("sensitivity_shift",myfiles))){
        myfiles <- myfiles[grep("sensitivity_shift",myfiles)]
        unlink(myfiles)
    }
}


# sp list
existing <- list.files(here::here("/media/seagate/boliveira/SDMs/MaxNet",realm))
missing_ones <- existing[-which(existing %in% SDMsSpList)]

# delete files
for(i in 1:length(missing_ones)){
    mydir <- here::here("/media/seagate/boliveira/SDMs/MaxNet",realm,missing_ones[i])
    unlink(mydir,recursive = TRUE)
}
