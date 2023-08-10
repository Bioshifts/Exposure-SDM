require("here")

# set resolution terrestrial environmental data
# must be 1km or 5km
my_res <- "1km"

# set terrestrial environmental data source
# must be cruts or chelsa
ter_data <- "cruts"
# set marine environmental data source
# must be oras or copernicus
mar_data <- "oras"

# Working directory
work_dir <- "/storage/simple/projects/t_cesab/brunno/Exposure-SDM"

# Bioshifts folder
Bioshifts_dir <- here::here("Data/Bioshifts")

# Bioshifts database
Bioshifts_DB <- "bioshifts_v3_raw.csv"

# Env data folder
env_data_dir <- function(realm){here::here(work_dir,"Data/Env_data",realm)}

# Occurrence data folder
occ_dir <- here::here(work_dir,"Data/GBIF_data")

# SA velocity folder
velocity_SA_dir <- here::here(work_dir,"Data/Velocity_SA")

# SA shapefiles folder
SA_shps_dir <- here::here(Bioshifts_dir,"ShapefilesBioShiftsv3")

temporal_range_env_data <- function(realm){
    if(realm == "Ter"){c(1901, 2016)} else {c(1958, 2019)}
}
min_year_layers <- function(realm){
    min(temporal_range_env_data(realm))
}
myvars <- function(realm){
    if(realm == "Ter"){c("tasmax", "tasmin", "tas", "pr")} else {"SST"}
}
vars_dir <- function(realm){
    if(realm == "Ter"){
        if(ter_data == "cruts"){
            here::here("/lustre/oliveirab/Data/Land/cruts")
        } else {
            here::here("/lustre/oliveirab/Data/Land/chelsa")
        }
    } else {
        if(mar_data == "oras"){
            here::here("/lustre/oliveirab/Data/Marine/oras")
        } else {
            here::here("/lustre/oliveirab/Data/Marine/copernicus")
        }
    }
}

n_tiles <- function(realm){
    if(realm == "Ter"){40} else {10}
}

# directory where range shifts results are saved
shift_dir <- function(realm){here::here("Data/SHIFT",realm)}

# directory for scratch data
scratch_dir <- "/lustre/oliveirab"

# directory where SDMs are saved
sdm_dir <- function(realm){here::here(scratch_dir,"SDMs",realm)}

# directory where range shift scripts are located
shift_script_dir <- here::here("R/6_shifts")

# directory where SDMs scripts are located
sdm_script_dir <- here::here("R/5_run_sdms")

#################

# basis of records for downloading GBIF data
basisOfRecord = c("HUMAN_OBSERVATION", "OBSERVATION", "OCCURRENCE")

# N years to calculate bioclimatics
n_yr_bioclimatic <- 1

# limit records for fitting SDMs
limit_recs = 20000

