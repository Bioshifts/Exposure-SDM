# Setup
rm(list=ls())
gc()

library(pbapply)
library(parallel)
library(tictoc)

source(here::here("R/my_functions.R"))

# Run code for terrestrials or marines?
realm <- "Mar"

# sp list
SDMsSpList <- list.files(here::here("/media/seagate/boliveira/SDMs/Env_data",realm))
infos <- lapply(SDMsSpList, function(x) strsplit(x,"_")[[1]])
SDMsSpList <- sapply(infos, function(x) paste(x[1],x[2],sep = "_"))

# submit job
n_jobs = 20 # n cores used at each job (internal parallelization at a sps-job)
# n_cores <- round((detectCores()-10) / n_jobs) # n of species possible to parallelize 
n_cores <- 10
    
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
              code_dir = here::here("R/5_get_sdm_run.R"),
              Rout_file = here::here(jobs.dir, paste0(sptogo, "_job.Rout")))
    
}, cl =cl)

stopCluster(cl)

toc()

################################
# check if job failed for any species

gotdata <- c()
for(i in 1:length(SDMsSpList)){
    gotdata[i] <- file.exists(here::here("/media/seagate/boliveira/SDMs/MaxNet",realm,SDMsSpList[i],paste0(SDMsSpList[i],"_model.RDS")))
}

gotdataSps <- SDMsSpList[which(gotdata)]

# missing species
missing_ones <- SDMsSpList[which(!gotdata)]

# N fitted SDMs
length(gotdataSps) # 344
# N missing SDMs
length(which(!gotdata)) # 68

if(length(missing_ones)>0){
    # submit job
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("systemjob","missing_ones","realm"))
    
    tic()
    pblapply(1:length(missing_ones), function(i){
        
        sptogo <- gsub(" ","_",missing_ones[i])
        args = paste(sptogo, realm)
        
        # create dir if doesn't exist
        jobs.dir <- here::here("/media/seagate/boliveira/SDMs/MaxNet",realm,sptogo)
        if(!dir.exists(jobs.dir)){
            dir.create(jobs.dir,recursive = T)
        }
        
        systemjob(args = args,
                  code_dir = here::here("R/5_get_sdm_run.R"),
                  Rout_file = here::here(jobs.dir, paste0(sptogo, "_job.Rout")))
        
    }, cl =cl)
    
    stopCluster(cl)
    
}

###