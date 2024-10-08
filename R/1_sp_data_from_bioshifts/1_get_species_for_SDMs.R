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
source("R/my_functions.R")

# Load bioshifts v3
bioshifts <- read.csv(here::here(Bioshifts_dir, Bioshifts_DB_v3), header = T)
table(bioshifts$class)

bioshifts <- bioshifts_sdms_selection(bioshifts)
table(bioshifts$class)

bioshifts$sp_name_std <- gsub("_"," ",bioshifts$sp_name_std)
all_sps <- unique(bioshifts$sp_name_std)
length(unique(all_sps)) # 4637

# load sp list
splist <- read.csv(here::here(Bioshifts_dir, "splist_v3.csv"))
all(all_sps %in% splist$scientificName)

splist <- splist %>% filter(db == "GBIF")
length(unique(splist$scientificName)) # 12360

# Filter species in the GBIF database
splist <- splist %>% filter(scientificName %in% all_sps)
length(unique(splist$scientificName)) # 4613

splist <- merge(splist, bioshifts[,c("sp_name_std","Eco")],
                by.x = "scientificName",
                by.y = "sp_name_std")

splist <- unique(splist)
length(unique(splist$scientificName)) # 4613

head(splist)
table(splist$class)

######################
# N occurrences GBIF

# Retrieve a summary table with N occurrence per species.
# We only retrieve N occurrences for records based:
# 1) "HUMAN_OBSERVATION" "OBSERVATION" "OCCURRENCE"
# 2) from which there was coordinates recorded
# 3) Registered after 1901 


ncores = parallelly::availableCores()

cl <- makeCluster(ncores)
clusterExport(cl, c("splist","basisOfRecord","temporal_range_env_data"))

N_OCC <- pbsapply(1:nrow(splist), function(i) {
    
    code_i <- gsub("GBIF:","",splist$db_code[i])
    tr_i <- temporal_range_env_data(splist$Eco[i])
    tr_i <- paste(tr_i[1],tr_i[2],sep = ",")
    
    records <- sapply(basisOfRecord, function(b) {
        rgbif::occ_count(taxonKey = code_i,
                         hasGeospatialIssue = FALSE,
                         basisOfRecord = b,
                         hasCoordinate = TRUE, 
                         occurrenceStatus = "PRESENT",
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

