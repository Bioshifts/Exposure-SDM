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

# source("R/my_functions.R")
# source("R/round-up.R")
source("R/settings.R")
source("R/Find_Sci_Names.R")
source("R/Clean_names.R")

# Load raw bioshifts
# Load v3
bioshifts <- read.csv(here::here(Bioshifts_dir, Bioshifts_DB_v3), header = T)

bioshifts <- bioshifts %>%
    filter(Type == "LAT") %>% # Use only latitudinal shifts
    filter(Eco == "Ter" & (MIDPOINT_firstperiod_v3 >= (temporal_range_env_data("Ter")[1] + n_yr_bioclimatic)) | # Shifts Marine or Terrestrial + within time period of the environmental data
               (Eco == "Mar" & (MIDPOINT_firstperiod_v3 >= (temporal_range_env_data("Mar")[1] + n_yr_bioclimatic)))) 

bioshifts$sp_reported_name_v1 <- gsub("_"," ",bioshifts$sp_reported_name_v1)
bioshifts$sp_reported_name_v1 <- Clean_Names(bioshifts$sp_reported_name_v1)
bioshifts$sp_reported_name_v1[which(bioshifts$sp_reported_name_v1=="NA")] <- NA

bioshifts$sp_reported_name_v2 <- gsub("_"," ",bioshifts$sp_reported_name_v2)
bioshifts$sp_reported_name_v2 <- Clean_Names(bioshifts$sp_reported_name_v2)
bioshifts$sp_reported_name_v2[which(bioshifts$sp_reported_name_v2=="NA")] <- NA

# N species
bioshifts$sp_name_v3 <- gsub("_"," ",bioshifts$sp_name_v3)
bioshifts$sp_name_v3 <- Clean_Names(bioshifts$sp_name_v3)
all_sps <- unique(bioshifts$sp_name_v3)
length(all_sps) # 5555


##############################
# Fix species names
mycols <- c("original_name","scientificName","kingdom","phylum","class","order","family","db_code")
sps_accepted_names <- data.frame(matrix(ncol = length(mycols), nrow = length(unique(all_sps))))
names(sps_accepted_names) <- mycols

sps_accepted_names$original_name <- unique(all_sps)

# retrieve sp names
sp_names_found <- Find_Sci_Names(sp_name = sps_accepted_names$original_name)
# took ~6 min on my personal computer


# ----------------
#  Summary 
# ----------------
# N taxa:
#     5555
# N taxa found:
#     |db   |    N|
#     |:----|----:|
#     |GBIF | 5525|
#     |ITIS |    1|
#     |NCBI |   21|
#     N taxa not found:
#     8

all(sp_names_found$requested_name %in% sps_accepted_names$original_name)

for(i in 1:length(sp_names_found$requested_name)){
    tofill <- unique(which(sps_accepted_names$original_name == sp_names_found$requested_name[i]))
    sps_accepted_names$scientificName[tofill] <- sp_names_found$scientificName[i]
    sps_accepted_names$kingdom[tofill] <- sp_names_found$kingdom[i]
    sps_accepted_names$phylum[tofill] <-sp_names_found$phylum[i]
    sps_accepted_names$class[tofill] <- sp_names_found$class[i]
    sps_accepted_names$order[tofill] <- sp_names_found$order[i]
    sps_accepted_names$family[tofill] <- sp_names_found$family[i]
    sps_accepted_names$db[tofill] <- sp_names_found$db[i]
    sps_accepted_names$db_code[tofill] <- sp_names_found$db_code[i]
}


splist <- sps_accepted_names

## *Classify species*
# Get taxonomic info for species missing
for(i in 1:nrow(splist)){ cat("\r",i, "from",nrow(splist))
    tmp <- splist[i,]
    
    if(any(is.na(tmp))){
        tmp_sp <- tmp$original_name
        missing_cols <- names(tmp[,c(2,3:7)])
        missing_cols <- missing_cols[which(is.na(tmp[,missing_cols]))]
        try_taxa <- try(wt_wikispecies(tmp_sp),silent = TRUE)
        while(class(try_taxa) == "try-error"){
            try_taxa <- try(wt_wikispecies(tmp_sp),silent = TRUE)
        }
        try_taxa <- try_taxa$classification
        if(!all(is.na(try_taxa$rank))){
            try_taxa <- data.frame(t(try_taxa))
            colnames(try_taxa) <- try_taxa[1,]
            try_taxa <- try_taxa[2,]
            try_taxa$kingdom <- ifelse(is.null(try_taxa$Regnum),NA,try_taxa$Regnum)
            try_taxa$phylum <- ifelse(is.null(try_taxa$Phylum),NA,try_taxa$Phylum)
            try_taxa$class <-ifelse(is.null(try_taxa$Classis),NA,try_taxa$Classis)
            try_taxa$order <- ifelse(is.null(try_taxa$Ordo),NA,try_taxa$Ordo)
            if("scientificName" %in% missing_cols){
                try_taxa$scientificName <- tmp_sp
            }
            if("family" %in% missing_cols){
                try_taxa2 <- wt_wikispecies_search(tmp_sp)
                try_taxa2 <- try_taxa2$query$search$snippet[1]
                try_taxa2 <- strsplit(try_taxa2,":")[[1]][2]
                try_taxa2 <- strsplit(try_taxa2," ")[[1]][2]
                try_taxa$family <- try_taxa2
            }
            tmp[,missing_cols] <- try_taxa[,missing_cols]
            splist[i,names(tmp)] <- tmp
            if("scientificName" %in% missing_cols){
                splist$db[i] <- "wiki"
                splist$db_code[i] <- wt_data_id(tmp_sp)[1]
            }
        } 
    }
}

# using itis
for(i in 1:nrow(splist)){ cat("\r",i, "from",nrow(splist))
    
    tmp <- splist[i,]
    
    if(any(is.na(tmp))){
        
        tmp_sp <- tmp$original_name
        missing_cols <- names(tmp[,c(2,3:7)])
        missing_cols <- missing_cols[which(is.na(tmp[,missing_cols]))]
        try_taxa <- taxize::tax_name(tmp_sp, get = missing_cols,messages = FALSE)
        if(length(missing_cols)>1){
            found_cols <- apply(try_taxa[,missing_cols],2,is.na)
            found_cols <- missing_cols[!found_cols] 
        } else {
            if(!is.na(try_taxa[,missing_cols])){
                found_cols <- missing_cols
            }
        }
        if(length(found_cols)>0){
            tmp[,found_cols] <- try_taxa[1,found_cols]
            splist[i,names(tmp)] <- tmp
        }
    }
}


# using NCBI
for(i in 1:nrow(splist)){ cat("\r",i, "from",nrow(splist))
    
    tmp <- splist[i,]
    
    if(any(is.na(tmp))){
        
        tmp_sp <- tmp$original_name
        missing_cols <- names(tmp[,c(2,3:7)])
        missing_cols <- missing_cols[which(is.na(tmp[,missing_cols]))]
        try_taxa <- taxize::tax_name(tmp_sp, get = missing_cols,messages = FALSE, db="ncbi")
        if(length(missing_cols)>1){
            found_cols <- apply(try_taxa[,missing_cols],2,is.na)
            found_cols <- missing_cols[!found_cols] 
        } else {
            if(!is.na(try_taxa[,missing_cols])){
                found_cols <- missing_cols
            }
        }
        if(length(found_cols)>0){
            tmp[,found_cols] <- try_taxa[1,found_cols]
            splist[i,names(tmp)] <- tmp
        }
    }
}

splist2 <- splist

# add original names
splist$sp_reported_name_v1 <- NA
splist$sp_reported_name_v2 <- NA

for(i in 1:nrow(splist)){ cat("\r",i, "from",nrow(splist))
    sp_i <- splist$original_name[i]
    v3_pos <- which(bioshifts$sp_name_v3 == sp_i)
    
    splist$sp_reported_name_v1[i] <- paste(na.omit(unique(bioshifts$sp_reported_name_v1[v3_pos])),collapse = ", ")
    splist$sp_reported_name_v2[i] <- paste(na.omit(unique(bioshifts$sp_reported_name_v2[v3_pos])),collapse = ", ")
}
splist$sp_reported_name_v1[which(splist$sp_reported_name_v1=="")] <- NA
splist$sp_reported_name_v2[which(splist$sp_reported_name_v2=="")] <- NA

any(is.na(splist$sp_reported_name_v1) & is.na(splist$sp_reported_name_v2))

all(splist$original_name %in% bioshifts$sp_name_v3)

## *Classify species*
splist$Group = NA
splist$kingdom[which(splist$class == "Phaeophyceae")] <- "Chromista"
splist$Group[which(splist$kingdom == "Chromista")] <- "Algae"
splist$Group[which(splist$phylum == "Rhodophyta")] <- "Algae"
splist$kingdom[which(splist$phylum == "Rhodophyta")] <- "Plantae"
splist$Group[which(splist$family == "Elminiidae")] <- "Barnacles"
splist$Group[which(splist$kingdom == "Bacteria")] <- "Bacteria"
splist$Group[which(splist$class == "Holothuroidea")] <- "Sea cucumber"
splist$Group[which(splist$class == "Aves")] <- "Bird"
splist$Group[which(splist$class == "Insecta")] <- "Insect"
splist$Group[which(splist$class == "Mammalia")] <- "Mammal"
splist$Group[which(splist$class == "Arachnida")] <- "Spider"
splist$kingdom[which(splist$kingdom == "Viridiplantae")] <- "Plantae"
splist$kingdom[which(splist$phylum == "Tracheophyta")] <- "Plantae"
splist$Group[which(splist$kingdom == "Plantae")] <- "Plant"
splist$Group[which(splist$class == "Hydrozoa")] <- "Hydrozoa"
splist$Group[which(splist$class == "Anthozoa")] <- "Sea anemones and corals"
splist$Group[which(splist$class == "Polychaeta")] <- "Polychaetes"
splist$Group[which(splist$phylum == "Mollusca")] <- "Molluscs"
splist$Group[which(splist$class == "Malacostraca")] <- "Crustacean"
splist$Group[which(splist$class == "Hexanauplia")] <- "Zooplankton"
splist$Group[which(splist$class == "Maxillopoda")] <- "Zooplankton"
splist$Group[which(splist$class == "Ostracoda")] <- "Zooplankton"
splist$Group[which(splist$class == "Branchiopoda")] <- "Crustacean"
splist$Group[which(splist$class == "Asteroidea")] <- "Starfish"
splist$Group[which(splist$class == "Ascidiacea")] <- "Ascidians tunicates and sea squirts"
splist$class[which(splist$class == "Actinopteri")] <- "Actinopterygii"
splist$Group[which(splist$class == "Actinopterygii")] <- "Fish"
splist$Group[which(splist$class == "Elasmobranchii")] <- "Fish"
splist$Group[which(splist$order == "Perciformes")] <- "Fish"
splist$Group[which(splist$class == "Chondrichthyes")] <- "Fish"
splist$Group[which(splist$class == "Holocephali")] <- "Fish"
splist$Group[which(splist$class == "Cephalaspidomorphi")] <- "Fish"
splist$Group[which(splist$class == "Echinoidea")] <- "Sea urchin"
splist$Group[which(splist$class == "Crinoidea")] <- "Crinoid"
splist$Group[which(splist$class == "Holothuroidea")] <- "Sea cucumber"
splist$Group[which(splist$class == "Reptilia")] <- "Reptile"
splist$Group[which(splist$order == "Squamata")] <- "Reptile"
splist$Group[which(splist$class == "Ophiuroidea")] <- "Brittle stars"
splist$Group[which(splist$class == "Chilopoda")] <- "Myriapoda"
splist$Group[which(splist$class == "Diplopoda")] <- "Myriapoda"
splist$Group[which(splist$class == "Amphibia")] <- "Amphibian"
splist$Group[which(splist$kingdom == "Fungi")] <- "Fungi"
splist$Group[which(splist$order == "Balanomorpha")] <- "Barnacles"
splist$Group[which(splist$phylum == "Nematoda")] <- "Nematodes"
splist$Group[which(splist$class == "Myxini")] <- "Hagfish"
splist$Group[which(splist$class == "Testudines")] <- "Turtle"
splist$Group[which(splist$class == "Teleostei")] <- "Fish"
splist$Group[which(splist$class == "Copepoda")] <- "Zooplankton"
splist$Group[which(splist$order == "Pleuronectiformes")] <- "Fish"


# Class species as Terrestrial, Marine or Freshwater
Terrestrials <- unique(bioshifts$sp_name_v3[which(bioshifts$Eco == "Ter")])
Marine <- unique(bioshifts$sp_name_v3[which(bioshifts$Eco == "Mar")])

# Freshwater fish
FFish <- unique(splist$original_name[which(splist$Group == "Fish")])
FFish <- FFish[which(FFish %in% Terrestrials)]

# Marine fish
MFish <- unique(splist$original_name[which(splist$Group == "Fish")])
MFish <- MFish[which(MFish %in% Marine)]

splist$Eco = NA
splist$Eco[which(splist$original_name %in% Terrestrials)] <- "Ter"
splist$Eco[which(splist$original_name %in% Marine)] <- "Mar"
splist$Eco[which(splist$original_name %in% MFish)] <- "Mar"
splist$Eco[which(splist$original_name %in% FFish)] <- "Aqua"

# save species list
write.csv(splist, "Data/Bioshifts/splist_v3_hamonized.csv", row.names = FALSE)
# splist <- read.csv("Data/Bioshifts/splist_v3_hamonized.csv")

######################
# N occurrences GBIF

# Retrieve a summary table with N occurrence per species.
# We only retrieve N occurrences for rEcords based:
# 1) "HUMAN_OBSERVATION" "OBSERVATION" "OCCURRENCE"
# 2) from which there was coordinates rEcorded
# 3) Registered after 1901 


ncores = parallelly::availableCores()

splist_GBIF <- splist %>% filter(db == "GBIF")

cl <- makeCluster(ncores)
clusterExport(cl, c("splist_GBIF","basisOfRecord","temporal_range_env_data"))

N_OCC <- pbsapply(1:nrow(splist_GBIF), function(i) {
    # N_OCC <- list()
    # for(i in 1:nrow(splist_GBIF)){
        
    code_i <- gsub("GBIF:","",splist_GBIF$db_code[i])
    tr_i <- temporal_range_env_data(splist_GBIF$Eco[i])
    tr_i <- paste(tr_i[1],tr_i[2],sep = ";")
    
    records <- sapply(basisOfRecord, function(b) {
        rgbif::occ_count(taxonKey = code_i,
                         hasGeospatialIssue = FALSE,
                         basisOfRecord = b,
                         year = tr_i) 
    })
    
    return(sum(records))
    
    # N_OCC[[i]] <- sum(records)
}
,
cl = cl)


stopCluster(cl)

# Took 9min my computer

splist_GBIF$N_OCC <- N_OCC

write.csv(splist_GBIF, "Data/N_OCC.csv", row.names = F)

N_OCC <- read.csv("Data/N_OCC.csv")

