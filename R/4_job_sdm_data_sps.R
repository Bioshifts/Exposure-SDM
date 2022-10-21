# Setup
rm(list=ls())
gc()

library(pbapply)
library(parallel)
library(tictoc)
library(dplyr)

# Run code for terrestrials or marines?
realm <- "Mar"

source(here::here("R/my_functions.R"))
source(here::here("R/settings.R"))

# sp list
SDMsSpList <- readRDS(here::here("Data/SDMsSpList.RDS"))

if(realm == "Mar"){
    SDMsSpList <- unique(SDMsSpList$species[which(SDMsSpList$ECO=="M")])
}
if(realm == "Ter"){
    SDMsSpList <- unique(SDMsSpList$species[which(SDMsSpList$ECO=="T")])
}
SDMsSpList <- gsub(" ","_",SDMsSpList)

# Get data for species with observed shifts documented for periods that we have environmental data & for LAT or ALT shifts
## Collect info on first and last periods
biov1 <- read.table(here::here("Data/Shifts2018_checkedtaxo.txt"),header = T)

# Use LAT ELE shifts
biov1$Type[which(!is.na(biov1$Azimuth))] <- "LAT" # All obs type == HOR contain Azimuth value
biov1 <- biov1[which(biov1$Type=="ELE" | biov1$Type=="LAT"),]

# Use temporal period from the environmental data
biov1 <- biov1 %>% filter(START >= temporal_range_env_data[1] + n_yr_bioclimatic)

keep <- unique(biov1$New_name)

SDMsSpList <- SDMsSpList[which(gsub(" ","_",SDMsSpList) %in% keep)]

# create dir if doesn't exist
jobs.dir <- here::here("/media/seagate/boliveira/SDMs/Env_data/jobs", realm)
if(!dir.exists(jobs.dir)){
    dir.create(jobs.dir,recursive = T)
}

# submit job
if(realm == "Mar"){
    n_tiles = 10 # for internal parallelization
    # n_cores <- round((detectCores()-10) / n_tiles) # species running at a time
    n_cores <- 20
}
if(realm == "Ter"){
    n_tiles = 40 # for internal parallelization
    # n_cores <- round((detectCores()-10) / n_tiles) # species running at a time
    n_cores <- 10 # species running at a time
}

# submit job
cl <- makeCluster(n_cores)
clusterExport(cl, c("systemjob","SDMsSpList","jobs.dir","realm"))

tic()
pblapply(1:length(SDMsSpList), function(i){
    
    sptogo <- gsub(" ","_",SDMsSpList[i])
    args = paste(sptogo, realm)
    
    systemjob(args = args,
              code_dir = here::here("R/4_get_sdm_data_sps.R"),
              Rout_file = here::here(jobs.dir, paste0(i,"_",sptogo, ".Rout")))
    
}, cl =cl)

stopCluster(cl)

toc()

################################
# check if job failed for any species

gotdata <- list.files(here::here("/media/seagate/boliveira/SDMs/Env_data",realm))
infos <- lapply(gotdata, function(x) strsplit(x,"_")[[1]])
gotdata <- sapply(infos, function(x) paste(x[1],x[2],sep = "_"))
# missing species
missing_ones <- SDMsSpList[which(!SDMsSpList %in% gotdata)]

any(missing_ones %in% gotdata)
gotdata[which(!gotdata %in% SDMsSpList)]

if(length(missing_ones)>0){
    # submit job
    cl <- makeCluster(n_cores)
    clusterExport(cl, c("systemjob","missing_ones","jobs.dir","realm"))
    
    tic()
    pblapply(1:length(missing_ones), function(i){
        
        sptogo <- gsub(" ","_",missing_ones[i])
        args = paste(sptogo, realm)
        
        systemjob(args = args,
                  code_dir = here::here("R/4_get_sdm_data_sps.R"),
                  Rout_file = here::here(jobs.dir, paste0(i,"_",sptogo, ".Rout")))
        
    }, cl =cl)
    
    stopCluster(cl)
    
    toc()
}


