rm(list=ls())
gc()


library(tictoc)
library(here)

########################
# set computer
computer = "matrics"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}

work_dir <- getwd()

source("R/settings.R")

########################

## 1) get list of marine species I have SDMs
sdms_mar <- list.dirs(sdm_dir("Mar"), recursive = FALSE, full.names = FALSE)

## 1.1) Load in validation datasets
CV_mar <- pbapply::pblapply(sdms_mar, function(x){
    file2go <- here::here(sdm_dir("Mar"),x,paste0(x,"_CV_ens_all.csv"))
    if(file.exists(file2go)){
        tmp <- read.csv(file2go)
    } else{
        file2go <- here::here(sdm_dir("Mar"),x,paste0(x,"_CV_ens_RUN.csv"))
        if(file.exists(file2go)){
            tmp <- read.csv(file2go)
        } 
    }
})
rem <- sapply(CV_mar, class)
rem <- which(rem=="try-error")
if(any(rem)){
    CV_mar <- CV_mar[-rem]
}
CV_mar <- data.table::rbindlist(CV_mar)
CV_mar$ECO = "Mar"

## 1) get list of terrestrial species I have SDMs
sdms_ter <- list.dirs(sdm_dir("Ter"), recursive = FALSE, full.names = FALSE)

## 1.1) Load in validation datasets
CV_ter <- pbapply::pblapply(sdms_ter, function(x){
    file2go <- here::here(sdm_dir("Ter"),x,paste0(x,"_CV_ens_all.csv"))
    if(file.exists(file2go)){
        tmp <- read.csv(file2go)
    } else{
        file2go <- here::here(sdm_dir("Ter"),x,paste0(x,"_CV_ens_RUN.csv"))
        if(file.exists(file2go)){
            tmp <- read.csv(file2go)
        } 
    }
})
rem <- sapply(CV_ter, class)
rem <- which(rem=="try-error")
if(any(rem)){
    CV_ter <- CV_ter[-rem]
}
CV_ter <- data.table::rbindlist(CV_ter)
CV_ter$ECO = "Ter"

## 2) group all
sdms_CV <- rbind(CV_mar, CV_ter)

## 3) save
write.csv(sdms_CV, here::here("Data","sdms_CV.csv"), row.names = FALSE)
