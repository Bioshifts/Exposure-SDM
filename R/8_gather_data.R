# Setup
rm(list=ls())
gc()

library(pbapply)
library(parallel)
library(tictoc)
library(dplyr)
library(ggrepel)

source(here::here("R/my_functions.R"))

# Run code for terrestrials or marines?
realm <- "Mar"

# Input directory
input.dir <- here::here("/media/seagate/boliveira/SDMs/MaxNet",realm)

# sp list
SDMsSpList <- list.files(input.dir)

# get bioshifts
sp_bioshifts <- list()
for(i in 1:length(SDMsSpList)){ cat("\r", i)
    sptogo <- SDMsSpList[i]
    bioshifts <- read.csv(here::here(input.dir,sptogo,paste(sptogo, "bioshift.csv")))
    sp_bioshifts[[i]] <- bioshifts
}
sp_bioshifts <- rbindlist(sp_bioshifts)

# get SDM shift
sp_SDM_shift <- pblapply(1:length(SDMsSpList), function(i) { 
    sptogo <- SDMsSpList[i]
    try({sh_SA <- read.csv(here::here(input.dir,sptogo,paste(sptogo, "shifts_SA.csv")))
    names(sh_SA)[2:4] <- paste0(names(sh_SA)[2:4],'_SA')
    sh_BG <- read.csv(here::here(input.dir,sptogo,paste(sptogo, "shifts_BG.csv")))
    names(sh_BG)[2:4] <- paste0(names(sh_BG)[2:4],'_BG')
    return(data.frame(sh_SA[,1:4],sh_BG[,-1]))
    })
})
# remove errors
rem <- sapply(sp_SDM_shift,class)
rem <- grep("error",rem)
sp_SDM_shift <- sp_SDM_shift[-rem]
sp_SDM_shift <- rbindlist(sp_SDM_shift)

# get exposure
sp_exposure <- pblapply(1:length(SDMsSpList), function(i) { 
    sptogo <- SDMsSpList[i]
    try({
        exp_bg <- read.csv(here::here(input.dir,sptogo,paste(sptogo, "exposure_BG.csv")))
        exp_sa <- read.csv(here::here(input.dir,sptogo,paste(sptogo, "exp_SA.csv")))
        sp_info <- read.csv(here::here(input.dir,sptogo,paste(sptogo, "shift_info.csv")))
        exp_sp <- cbind(exp_bg,exp_sa,sp_info[,c(1:3,5)])
        names(exp_sp)[1:2] <- c("exp_BG","exp_SA")
        return(exp_sp)
    })
})
# remove errors
rem <- sapply(sp_exposure,class)
rem <- grep("error",rem)
sp_exposure <- sp_exposure[-rem]
sp_exposure <- rbindlist(sp_exposure)

# range position
sp_ran_pos <- pblapply(1:length(SDMsSpList), function(i) { 
    sptogo <- SDMsSpList[i]
    try({
        rp <- read.csv(here::here(input.dir,sptogo,paste(sptogo, "range_position.csv")))
        shift_info <- read.csv(here::here(input.dir,sptogo,paste(sptogo, "shift_info.csv")))
        rp <- cbind(rp,shift_info)
        return(rp)
    })
})
# remove errors
rem <- sapply(sp_ran_pos,class)
rem <- grep("error",rem)
sp_ran_pos <- sp_ran_pos[-rem]
sp_ran_pos <- rbindlist(sp_ran_pos)

# SDM output
sp_SDM_out <- pblapply(1:length(SDMsSpList), function(i) { 
    sptogo <- SDMsSpList[i]
    try({
        rp <- read.csv(here::here(input.dir,sptogo,paste0(sptogo, "_SDM_CV.csv")))
        rp$Species <- sptogo
        return(rp)
    })
})
# remove errors
rem <- sapply(sp_SDM_out,class)
rem <- grep("error",rem)
sp_SDM_out <- sp_SDM_out[-rem]
sp_SDM_out <- rbindlist(sp_SDM_out)

######################
# Save data

# Bioshifts data
write.csv(sp_bioshifts, 
          here::here("Output",paste(realm,"sps_bioshifts.csv")), 
          row.names = FALSE)
# SDM outputs data
write.csv(sp_SDM_out, 
          here::here("Output",paste(realm,"sp_SDM_out.csv")), 
          row.names = FALSE)
# range position
write.csv(sp_ran_pos, 
          here::here("Output",paste(realm,"sp_range_pos.csv")), 
          row.names = FALSE)
# Exposure
write.csv(sp_exposure, 
          here::here("Output",paste(realm,"sp_exposure.csv")), 
          row.names = FALSE)
# SDM shifts
write.csv(sp_SDM_shift, 
          here::here("Output",paste(realm,"sp_SDM_shift.csv")), 
          row.names = FALSE)

######################
                                   