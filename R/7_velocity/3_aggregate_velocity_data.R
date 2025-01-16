
### P.S.: Cannot run this inside singularity container
### On terminal, open new screen, load R and run this!!!

rm(list=ls())
gc()

library(tictoc)
library(here)
library(terra)

########################
# set computer
computer = "matrics"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    work_dir <- getwd()
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
    work_dir <- getwd()
}

# source settings
source("R/settings.R")
source("R/my_functions.R")

# Load study areas v3
v3_polygons <- gsub(".shp","",list.files(SA_shps_dir,pattern = ".shp"))
length(v3_polygons)

Bioshifts_DB <- read.csv(here::here(Bioshifts_dir,Bioshifts_DB_v3))
Bioshifts_DB <- bioshifts_fix_columns(Bioshifts_DB)
Bioshifts_DB <- Bioshifts_DB %>%
    mutate(new_eco = ifelse(
        Eco == "Mar", "Mar", ifelse(Eco != "Mar", "Ter", NA))) %>%
    mutate(Eco = new_eco)

# get useful columns
Bioshifts_DB <- Bioshifts_DB %>%
    dplyr::select(ID, Eco, Areakm2, Type) %>%
    unique()

# load velocities
# get info on resolution
vel_SA_resolution <- list.files(velocity_SA_dir, pattern = ".csv")
vel_SA_resolution <- sapply(vel_SA_resolution, function(x){
    strsplit(gsub(".csv","",x),"_")[[1]][3]
})
# load in all velocities
vel_SA <- list.files(velocity_SA_dir, pattern = ".csv", full.names = TRUE)
vel_SA <- lapply(vel_SA, read.csv)
vel_SA <- data.table::rbindlist(vel_SA, fill = TRUE)
# add info on resolution
vel_SA$res <- vel_SA_resolution
dim(vel_SA)
# add info on Eco and area
vel_SA <- merge(vel_SA, unique(Bioshifts_DB[,c("ID","Eco","Areakm2")]), by = "ID")

write.csv(vel_SA, 
          here(Bioshifts_dir,"vel_SA_all.csv"), 
          row.names = FALSE)


# #################################################
# # recalculate SA summary if any is missing
# 
# mar_cols <- colnames(vel_SA)
# mar_cols <- mar_cols[grep("sst",mar_cols)]
# mar_cols <- mar_cols[grep("v.",mar_cols)]
# 
# ter_cols <- colnames(vel_SA)
# ter_cols <- ter_cols[grep("mat|map",ter_cols)]
# ter_cols <- ter_cols[grep("v.",ter_cols)]
# 
# for(i in 1:nrow(vel_SA)){
#     
#     SA_i <- vel_SA$ID[i]
#     Eco_i <- vel_SA$Eco[i]
#     res_i <- vel_SA$res[i]
#     
#     SA_i <- "A46_P1"
#     Eco_i <- "Mar"
#     res_i <- "110km"
#     
#     if(Eco_i == "Ter"){
#         cols2go <- ter_cols
#     } else {
#         cols2go <- mar_cols
#     }
#     
#     vel_SA_i <- vel_SA %>% filter(ID == SA_i, res == res_i) %>% dplyr::select(cols2go) %>% data.frame
#     
#     if(any(is.na(vel_SA_i))){
#         
#         for(j in 1:ncol(vel_SA_i)){
#             
#             vel_SA_i_var_j <- vel_SA_i[,j]
#             vel_SA_i_var_j_name <- names(vel_SA_i)[j]
#             
#             if(is.na(vel_SA_i_var_j)){
#                 
#                 # make var
#                 vel_SA_i_var_j_info <- strsplit(vel_SA_i_var_j_name,".",fixed = TRUE)[[1]]
#                 
#                 # var velocity
#                 if(any("lat" %in% vel_SA_i_var_j_info)){
#                     var_vel <- "gVelLat"
#                 } else {
#                     var_vel <- "gVel"
#                 }
#                 
#                 # var climate
#                 if(any("sst" %in% vel_SA_i_var_j_info)){
#                     var_clim <- "sst"
#                 } else {
#                     if(any("mat" %in% vel_SA_i_var_j_info)){
#                         var_clim <- "mat"
#                     } else {
#                         var_clim <- "map"
#                     }
#                 }
#                 
#                 # var summary
#                 if(any("mean" %in% vel_SA_i_var_j_info)){
#                     var_summ <- "mean"
#                 } else {
#                     var_summ <- "median"
#                 }
#                 
#                 
#                 vel_SA_i_var_j_get <- paste0(paste(SA_i,var_clim,res_i,var_vel,sep = "_"),".tif")
#                 # var j 
#                 var_j <- rast(here(velocity_SA_dir,vel_SA_i_var_j_get))
#                 var_j
#                 plot(var_j);dev.off()
#             }
#             
#         }
#         
#     }
# }

