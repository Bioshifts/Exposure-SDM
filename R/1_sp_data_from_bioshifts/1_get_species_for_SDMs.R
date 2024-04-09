# title: "Get species data for SDMs"
# author: "Brunno F Oliveira & Bioshifts group"


# Setup

rm(list=ls())
gc()

list.of.packages <- c("stringr", "rgbif", "parallel", "pbapply", "raster", "dplyr", "here", "tidyr",
                      "data.table", "ggplot2", "readxl", "knitr","googledrive","googleAuthR",
                      "sqldf","RSQLite", "scales","wikitaxa")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

########################
# set computer
computer = "personal"

# Load functions
source("R/settings.R")
source("R/bioshifts_fix_columns.R")

# load sp list
splist <- read.csv(here::here(Bioshifts_dir, "splist_v3.csv"))
length(unique(splist$scientificName)) # 12410

# Load raw bioshifts
# Load v3
bioshifts <- read.csv(here::here(Bioshifts_dir, Bioshifts_DB_v3), header = T)
bioshifts <- bioshifts_fix_columns(bioshifts)
length(unique(bioshifts$sp_name_std)) # 12410

# Select latitude shifts 
bioshifts <- bioshifts %>%
    filter(Type == "LAT") 
length(unique(bioshifts$sp_name_std)) # 6257

# select shifts within the temporal extent of environmental data
bioshifts <- bioshifts %>%
    filter(Eco == "Ter" & (Start >= (temporal_range_env_data("Ter")[1] + n_yr_bioclimatic)) | # Shifts Marine or Terrestrial + within time period of the environmental data
               (Eco == "Mar" & (Start >= (temporal_range_env_data("Mar")[1] + n_yr_bioclimatic)))) 
length(unique(bioshifts$sp_name_std)) # 5508

all_sps <- unique(bioshifts$sp_name_std)

splist <- splist %>% filter(scientificName %in% all_sps)

splist <- merge(splist, bioshifts[,c("sp_name_std","Eco")],
                by.x = "scientificName",
                by.y = "sp_name_std",
                all.x = TRUE)

splist <- unique(splist)
length(unique(splist$db_code)) # 5508

######################
# N occurrences GBIF

# Retrieve a summary table with N occurrence per species.
# We only retrieve N occurrences for records based:
# 1) "HUMAN_OBSERVATION" "OBSERVATION" "OCCURRENCE"
# 2) from which there was coordinates recorded
# 3) Registered after 1901 


ncores = parallelly::availableCores()

splist_GBIF <- splist %>% filter(db == "GBIF")

cl <- makeCluster(ncores)
clusterExport(cl, c("splist_GBIF","basisOfRecord","temporal_range_env_data"))

N_OCC <- pbsapply(1:nrow(splist_GBIF), function(i) {
    
    code_i <- gsub("GBIF:","",splist_GBIF$db_code[i])
    tr_i <- temporal_range_env_data(splist_GBIF$Eco[i])
    tr_i <- paste(tr_i[1],tr_i[2],sep = ",")
    
    records <- sapply(basisOfRecord, function(b) {
        rgbif::occ_count(taxonKey = code_i,
                         hasGeospatialIssue = FALSE,
                         basisOfRecord = b,
                         year = tr_i) 
    })
    
    return(sum(records))
}
,
cl = cl)


stopCluster(cl)

# Took 11 min my computer

splist_GBIF$N_OCC <- N_OCC

write.csv(splist_GBIF, here("Data/N_OCC.csv"), row.names = F)

############################

