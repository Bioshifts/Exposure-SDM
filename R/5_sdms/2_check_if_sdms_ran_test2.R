## 1) Run this in singularity 
## 2) save the missing species matrix
## 3) Run regular R (outside singularity)
## 4) load missing species matrix
## 5) refit sdms on missing species


library(terra)

# get list of missing species
# these are those with errors loading a projection raster from sdms or with NAs in raster values

computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    
    work_dir <- getwd()
    
    scratch_dir <- "/lustre/oliveirab"
    sdm_dir <- here::here(scratch_dir,"SDMs")
}


# sp list with environmental data
all_mar <- list.dirs(here::here(sdm_dir,"Mar"), recursive = FALSE, full.names = FALSE)

test <- pbapply::pbsapply(1:length(all_mar), function(i){
    x <- here::here(here::here(sdm_dir,"Mar",all_mar[i],gsub("_",".",all_mar[i])))
    tmp <- list.files(x,full.names = TRUE)
    tmp <- tmp[grep(" SA start ens",tmp)][1]
    if(is.na(tmp)){
        FALSE
    } else {
        tmp <- list.files(tmp, pattern = ".out", full.names = TRUE)
        if(file.exists(tmp)){
            tmp <- get(load(tmp))
            tmp = tmp@proj.out@val
            tmp = rast(tmp)
            any(is.na(terra::minmax(tmp[[1]])))
        } else {
            FALSE
        }
    }
})
missing_mar <- all_mar[test]

all_ter <- list.dirs(here::here(sdm_dir,"Ter"), recursive = FALSE, full.names = FALSE)

test <- pbapply::pbsapply(1:length(all_ter), function(i){
    x <- here::here(here::here(sdm_dir,"Ter",all_ter[i],gsub("_",".",all_ter[i])))
    tmp <- list.files(x,full.names = TRUE)
    tmp <- tmp[grep(" SA start ens",tmp)][1]
    if(is.na(tmp)){
        FALSE
    } else {
        tmp <- list.files(tmp, pattern = ".out", full.names = TRUE)
        if(file.exists(tmp)){
            tmp <- get(load(tmp))
            tmp = tmp@proj.out@val
            tmp = rast(tmp)
            any(is.na(terra::minmax(tmp[[1]])))
        } else {
            FALSE
        }
    }
}, cl = 10)
test <- ifelse(test=="FALSE", FALSE, TRUE)
missing_ter <- all_ter[test]


# missing species
missing_sps <- rbind(
    data.frame(realm = "Mar",
               sps = missing_mar),
    data.frame(realm = "Ter",
               sps = missing_ter)
)

dim(missing_sps)


write.csv(missing_sps,"missing_sps.csv",row.names = FALSE)
