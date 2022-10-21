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

# sp list
SDMsSpList <- list.files(here::here("/media/seagate/boliveira/SDMs/MaxNet",realm))

# submit job
n_cores <- 60

# submit job
cl <- makeCluster(n_cores)
clusterExport(cl, c("systemjob","SDMsSpList","realm"))

tic()
pblapply(1:length(SDMsSpList), function(i){
    
    # sptogo <- "Ampelisca_abdita"
    sptogo <- gsub(" ","_",SDMsSpList[i])
    args = paste(sptogo, realm)
    
    # create dir if doesn't exist
    jobs.dir <- here::here("/media/seagate/boliveira/SDMs/MaxNet",realm,sptogo)
    if(!dir.exists(jobs.dir)){
        dir.create(jobs.dir,recursive = T)
    }
    
    systemjob(args = args,
              code_dir = here::here("R/7_get_shifts.R"),
              Rout_file = here::here(jobs.dir, paste0(sptogo, "_job.Rout")))
    
}, cl =cl)

stopCluster(cl)

toc()
