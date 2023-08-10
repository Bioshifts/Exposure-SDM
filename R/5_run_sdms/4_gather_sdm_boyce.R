rm(list=ls())
gc()


library(tictoc)
library(here)

########################
# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    
    work_dir <- getwd()
    scratch_dir <- "/lustre/oliveirab"
    sdm_dir <- here::here(scratch_dir,"SDMs")
}

########################

## 1) get list of marine species I have SDMs
sdms_mar <- list.dirs(here::here(sdm_dir,"Mar"), recursive = FALSE, full.names = FALSE)

## 2) loop trough species

sdms_mar_boyce <- list()
for(i in 1:length(sdms_mar)){
    
    sdms_mar_boyce[[i]] <- try({
        tmp_dir <- here::here(sdm_dir,"Mar",sdms_mar[i],gsub("_",".",sdms_mar[i]))
        tmp <- list.files(tmp_dir, pattern = "models.out")
        tmp_ens <- tmp[grep("ensemble",tmp)]
        tmp_out <- tmp[-grep("ensemble",tmp)]
        # if more then one model use the most updated one
        if(length(tmp_ens)>1){
            tmp2 <- sapply(tmp_ens, function(x) {
                strsplit(x,"[.]")[[1]][3]
            })
            tmp_ens <- tmp_ens[order(tmp2,decreasing = TRUE)[1]]
        }
        if(length(tmp_out)>1){
            tmp2 <- sapply(tmp_out, function(x) {
                strsplit(x,"[.]")[[1]][3]
            })
            tmp_out <- tmp_out[order(tmp2,decreasing = TRUE)[1]]
        }
        
        bm.mod = get(load(here::here(tmp_dir,tmp_out)))
        bm.em = get(load(here::here(tmp_dir,tmp_ens)))
        
        res_i = BIOMOD_PresenceOnly(bm.mod = bm.mod,
                                    bm.em = bm.em)
        
    })
    
}


