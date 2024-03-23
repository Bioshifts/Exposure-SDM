rm(list=ls())
gc()


library(tictoc)
library(here)

########################
# set computer
computer = "matrics"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    
    work_dir <- getwd()
    scratch_dir <- "/lustre/oliveirab"
    sdm_dir <- here::here(scratch_dir,"SDMs")
}

########################

## get list of marine species I have SDMs
sdms_mar <- list.dirs(here::here(sdm_dir,"Mar"), recursive = FALSE, full.names = FALSE)
sdms_ter <- list.dirs(here::here(sdm_dir,"Ter"), recursive = FALSE, full.names = FALSE)

sdms_mar <- sapply(1:length(sdms_mar), function(i) {
    
    # check if there a file with ensemble model outputs
    tmp = list.files(here::here(sdm_dir,
                                "Mar",
                                sdms_mar[i]))
    tmp = any(grepl("ens_SDM",tmp))
    
    tmp = !tmp
    
    return(tmp)
    
})


if(any(sdms_mar)){ # any missing dir?
    cat("There are", length(sdms_mar), "marine species:", 
        "\nSDMs fitted for:",length(which(!sdms_mar)), "species",
        "\nSDMs missing for:",length(which(sdms_mar)), "species\n")
} else {
    cat("I have SDMs for all",length(which(sdms_mar)), "marine species from which environmental data was extracted")
}


## get list of terrestrial species I have SDMs
sdms_ter <- sapply(1:nrow(ter_sps), function(i) {
    
    # check if dir is missing
    tmp = !dir.exists(here::here(sdm_dir,
                                 "Ter",
                                 gsub(" ","_",ter_sps$sps[i]),
                                 gsub(" ",".",ter_sps$sps[i])))
    
    if(tmp){
        tmp
    } else {
        # check if there a file with ensemble model outputs
        tmp = list.files(here::here(sdm_dir,
                                    "Ter",
                                    gsub(" ","_",ter_sps$sps[i])))
        tmp = any(grepl("ens_SDM",tmp))
        
        tmp = !tmp
    }
    
    return(tmp)
    
})

if(any(sdms_ter)){ # any missing dir?
    cat("There are", length(sdms_ter), "terrestrial species:", 
        "\nSDMs fitted for:",length(which(!sdms_ter)), "species",
        "\nSDMs missing for:",length(which(sdms_ter)), "species\n")
} else {
    cat("I have SDMs for all",length(which(sdms_ter)), "terrestrial species from which environmental data was extracted")
}