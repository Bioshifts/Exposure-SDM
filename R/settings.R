require("here")

# -----------------------------
## Settings

# set resolution terrestrial environmental data
# must be 1km or 5km
my_res <- "1km"

# set terrestrial environmental data source
# must be cruts or chelsa
ter_data <- "cruts"

# set marine environmental data source
# must be oras or copernicus
mar_data <- "oras"

# Bioshifts database
Bioshifts_DB_v1 <- "biov1_fixednames.csv"
Bioshifts_DB_v2 <- "biov2_fixednames.csv"
Bioshifts_DB_v3 <- "BIOSHIFTS_v3_24032024_methodo_rate.csv"

# basis of records for downloading GBIF data
basisOfRecord = c("HUMAN_OBSERVATION", "OBSERVATION", "OCCURRENCE")

# N years to calculate bioclimatics
n_yr_bioclimatic <- 1

# limit records for fitting SDMs
limit_recs = 20000

# Eckert equal-area projection
Eckt <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" 

# -------------------------------
## Directories 

# Working directory
work_dir <- if(computer == "muse"){
    "/storage/simple/projects/t_cesab/brunno/Exposure-SDM"
}else{
    if(computer=="matrics"){
        "/users/boliveira/Exposure-SDM"
    } else {
        if(computer=="personal"){
            here::here()
        }
    }
}

# Scratch data
scratch_dir <- if(computer == "muse"){
    "/lustre/oliveirab"
}else{
    if(computer=="matrics"){
        "/scratch/boliveira"
    } else {
        if(computer=="personal"){
            here::here()
        }
    }
}

# Temporary data
tmp_dir <- here::here(scratch_dir,"tmp")

# Bioshifts dir
Bioshifts_dir <- here::here(work_dir,"Data/Bioshifts")

# Env data dir
env_data_dir <- function(realm){here::here(work_dir,"Data/Env_data",realm)}

# Env variables dir
vars_dir <- function(realm){
    if(realm == "Ter"){
        return(here::here(scratch_dir,"Data/Land",ter_data))
    }
    if(realm == "Mar"){
        return(here::here(scratch_dir,"Data/Marine",mar_data))
    }
}

# Bioclimatics dir
bios_dir <- function(realm){
    if(realm == "Ter"){
        return(paste0(vars_dir(realm),"/bio_proj_",my_res))
    } 
    if(realm == "Mar"){
        return(paste0(vars_dir(realm),"/bio_proj"))
    }
}

# Bioclimatics at the SA dir
bios_SA_dir <- function(realm){
    return(paste0(bios_dir(realm),"_SA"))
}

# Occurrence data dir
occ_dir <- here::here(work_dir,"Data/GBIF_data")

# SA velocity dir
velocity_SA_dir <- here::here(work_dir,"Data/Velocity_SA")

# SA velocity script dir
velocity_SA_scrit_dir <- here::here(work_dir,"R/7_velocity")

# SA shapefiles dir
SA_shps_dir <- here::here(Bioshifts_dir,"ShapefilesBioShiftsv3")

# Range shift dir
shift_dir <- function(realm){here::here(work_dir,"Data/SHIFT",realm)}

# Range shift plot
shiftplot_dir <- function(realm){here::here(work_dir,"Data/SHIFTplot",realm)}

# Range shift scripts dir
shift_script_dir <- here::here(work_dir,"R/6_shifts")

# SDMs dir
sdm_dir <- function(realm){here::here(scratch_dir,"SDMs",realm)}

# SDMs scripts dir
sdm_script_dir <- here::here(work_dir,"R/5_sdms")

# Singularity image
singularity_image <- if(computer == "muse"){
    "/storage/simple/projects/t_cesab/brunno/brunnospatial.sif"
}else{
    if(computer=="matrics"){
        "/users/boliveira/brunnospatial.sif"
    }
}

# -------------------------------
## Functions
temporal_range_env_data <- function(realm){
    if(realm == "Ter"){c(1901, 2016)} else {c(1958, 2019)}
}
min_year_layers <- function(realm){
    min(temporal_range_env_data(realm))
}
myvars <- function(realm){
    if(realm == "Ter"){c("tasmax", "tasmin", "tas", "pr")} else {"SST"}
}

n_tiles <- function(realm){
    if(realm == "Ter"){40} else {10}
}




