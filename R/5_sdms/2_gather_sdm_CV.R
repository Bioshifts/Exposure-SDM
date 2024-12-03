rm(list=ls())
gc()


library(tictoc)
library(here)
library(biomod2)

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
source("R/my_functions.R")

########################
bioshifts <- read.csv(here(Bioshifts_dir,Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_sdms_selection(bioshifts)

## 1) get list of marine species I have SDMs
sdms_mar <- list.dirs(sdm_dir("Mar"), recursive = FALSE, full.names = FALSE)
sdms_mar <- sdms_mar[which(sdms_mar %in% bioshifts$sp_name_std)]

## 1.1) Load in validation datasets
CV_mar <- pbapply::pblapply(sdms_mar, function(x){
    # file2go <- here::here(sdm_dir("Mar"),x,paste0(x,"_CV_ens_all.csv"))
    # if(file.exists(file2go)){
    #     tmp <- read.csv(file2go)
    # } else{
    file2go <- here::here(sdm_dir("Mar"),x,paste0(x,"_CV.csv"))
    if(file.exists(file2go)){
        tmp <- read.csv(file2go)
    }
    # }
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
sdms_ter <- sdms_ter[which(sdms_ter %in% bioshifts$sp_name_std)]

## 1.1) Load in validation datasets
CV_ter <- pbapply::pblapply(sdms_ter, function(x){
    # file2go <- here::here(sdm_dir("Ter"),x,paste0(x,"_CV_ens_all.csv"))
    # if(file.exists(file2go)){
    #     tmp <- read.csv(file2go)
    # } else{
    file2go <- here::here(sdm_dir("Ter"),x,paste0(x,"_CV.csv"))
    if(file.exists(file2go)){
        tmp <- read.csv(file2go)
    }
    # }
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


####################### BOYCE #########################

library(parallel)
cl <- makeCluster(detectCores()-3)
clusterExport(cl, c("sdms_mar", "sdms_ter","apply_gsub_to_s4","env_data_dir","computer","work_dir","sdm_dir","scratch_dir"))


# Calculate Boyce index for presence only data
Boyce_mar <- pbapply::pblapply(sdms_mar, function(x){
    
    try({
        
        # Load predictions
        # output dir
        output_dir <- here::here(sdm_dir("Mar"), x)
        if(!dir.exists(output_dir)){
            dir.create(output_dir,recursive = TRUE)
        }
        files_sdms <- list.files(
            here::here(output_dir,gsub("_",".",x)), 
            full.names = TRUE,
            pattern = "models.out")  
        
        pos <- grep("ensemble",files_sdms)
        files_sdms_ensemble <- files_sdms[pos]
        files_sdms <- files_sdms[-pos]
        
        model_sp <- get(load(files_sdms))
        ens_model_sp <- get(load(files_sdms_ensemble))
        
        model_sp <- apply_gsub_to_s4(model_sp, pattern = "/lustre/oliveirab", replacement = "/scratch/boliveira")
        ens_model_sp <- apply_gsub_to_s4(ens_model_sp, pattern = "/lustre/oliveirab", replacement = "/scratch/boliveira")
        
        ens_pred <- ens_model_sp@models.prediction@val
        
        # Load observations
        PresAbs <- qs::qread(here::here(sdm_dir("Mar"),x,paste0(x,"_PresAbs_PC.qs")))
        PresAbs <- PresAbs$pa
        
        PresAbs <- PresAbs[ens_pred$points]
        ens_pred <- ens_pred$pred
        
        # length(ens_pred)
        # length(PresAbs)
        
        data.frame(species = x,
                   biomod2::bm_FindOptimStat(metric.eval = "BOYCE", obs = PresAbs, fit = ens_pred))
    },
    silent = TRUE)
    
}, cl = cl)
rem <- sapply(Boyce_mar, class)
rem <- which(rem=="try-error")
if(any(rem)){
    Boyce_mar <- Boyce_mar[-rem]
}

Boyce_mar <- data.table::rbindlist(Boyce_mar, fill = TRUE)
Boyce_mar$ECO = "Mar"

nrow(Boyce_mar)
nrow(na.omit(Boyce_mar))
length(sdms_mar)


# Calculate Boyce index for presence only data
Boyce_ter <- pbapply::pblapply(sdms_ter, function(x){
    
    try({
        
        # Load predictions
        # output dir
        output_dir <- here::here(sdm_dir("Ter"), x)
        if(!dir.exists(output_dir)){
            dir.create(output_dir,recursive = TRUE)
        }
        files_sdms <- list.files(
            here::here(output_dir,gsub("_",".",x)), 
            full.names = TRUE,
            pattern = "models.out")  
        
        pos <- grep("ensemble",files_sdms)
        files_sdms_ensemble <- files_sdms[pos]
        files_sdms <- files_sdms[-pos]
        
        model_sp <- get(load(files_sdms))
        ens_model_sp <- get(load(files_sdms_ensemble))
        
        model_sp <- apply_gsub_to_s4(model_sp, pattern = "/lustre/oliveirab", replacement = "/scratch/boliveira")
        ens_model_sp <- apply_gsub_to_s4(ens_model_sp, pattern = "/lustre/oliveirab", replacement = "/scratch/boliveira")
        
        ens_pred <- ens_model_sp@models.prediction@val
        
        # Load observations
        PresAbs <- qs::qread(here::here(sdm_dir("Ter"),x,paste0(x,"_PresAbs_PC.qs")))
        PresAbs <- PresAbs$pa
        
        PresAbs <- PresAbs[ens_pred$points]
        ens_pred <- ens_pred$pred
        
        # length(ens_pred)
        # length(PresAbs)
        
        data.frame(species = x,
                   biomod2::bm_FindOptimStat(metric.eval = "BOYCE", obs = PresAbs, fit = ens_pred))
    },
    silent = TRUE)
    
}, cl = cl)
rem <- sapply(Boyce_ter, class)
rem <- which(rem=="try-error")
if(any(rem)){
    Boyce_ter <- Boyce_ter[-rem]
}

Boyce_ter <- data.table::rbindlist(Boyce_ter, fill = TRUE)
Boyce_ter$ECO = "Ter"

nrow(Boyce_ter)
nrow(na.omit(Boyce_ter))
length(sdms_ter)

## 2) group all
sdms_Boyce <- rbind(Boyce_mar, Boyce_ter)

## 3) save
write.csv(sdms_Boyce, here::here(getwd(),"Data","sdms_boyce.csv"), row.names = FALSE)

stopCluster(cl)
