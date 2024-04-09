# title: "Get GBIF data"
# author: "Brunno F Oliveira & Bioshifts group"


# Setup ----


rm(list=ls())
gc()

list.of.packages <- c("stringr", "crul", "dplyr", "here", "qs",
                      "parallel", "pbapply", 
                      "data.table", "ggplot2", "disk.frame",
                      "readxl", "knitr", 
                      "speciesgeocodeR", "CoordinateCleaner", "rgbif", "raster", "terra",
                      "sqldf","RSQLite","jsonlite")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)


########################
# set computer ----
computer = "matrics"
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}
work_dir <- getwd()

# Load functions
source("R/settings.R")
source("R/decompress_file.R")

decompress_file
range_env_data <- range(c(temporal_range_env_data("Ter"),temporal_range_env_data("Mar")))

# Get species list ----
N_OCC <- read.csv("Data/N_OCC.csv")
nrow(N_OCC)

N_OCC <- N_OCC %>% filter(N_OCC >= 30)
nrow(N_OCC)

# filter species I still do not have downloaded GBIF data
I_have_GBIF <- list.files("Data/GBIF_data")
I_have_GBIF <- gsub(".qs","",I_have_GBIF)
I_have_GBIF <- gsub("_"," ",I_have_GBIF)

N_OCC <- N_OCC %>% filter(!scientificName %in% I_have_GBIF)
nrow(N_OCC)

# GBIF requests ----

# Create requests for all taxa keys. The request is processed by GBIF. The processing of each request takes a while (GBIF website says it can take up to 15min). Once the request is processed, the file is ready for download under my user page at the GBIF website. The advantage of this method is reproducibility. It generates a DOI that can be cited and link for download that can be shared. Creating a single request can generate a very large file (File size was 300Gbs zipped == 600Gbs unziped). Thus, better submitting multiple requests. Here, we created divided species into chunks and loop through chunks for downloading data.

## Create queries


my_keys = N_OCC %>%
    group_by(db_code) %>%
    dplyr::summarise(N=sum(N_OCC))

my_keys$db_code <- gsub("GBIF:","",my_keys$db_code)

# split the keys to create N requests
reqs <- list()
# Max N occ per request
maxN <- 10^7
k = 1
n_obs <- 0
my_keys_n <- c()

for(i in 1:nrow(my_keys)){ cat(i, "from", nrow(my_keys), "\r")
    
    # get N Occ for key i
    ti <- my_keys$N[i]
    # sum N_occ for key N
    n_obs = n_obs + ti 
    my_keys_n <- c(my_keys_n, my_keys$db_code[i]) # store list of keys
    
    if(n_obs >= maxN) { # of sum N_occ >= maxN, save list of keys in chunk
        reqs[[k]] <- my_keys_n
        k = k + 1
        n_obs <- 0
        my_keys_n <- NULL
    }
}

if(length(reqs)==0){
    reqs[[1]] <- my_keys_n
}

# N keys
nrow(my_keys)
# Number of chunks
length(reqs)
# Number of keys per chunk
sapply(reqs, length)
# Avg of keys per chunk
mean(sapply(reqs, length))
# Max of keys per chunk
max(sapply(reqs, length))
# Min of keys per chunk
min(sapply(reqs, length))

# N obs per chunk
sapply(reqs, function(x){
    my_keys %>%
        filter(db_code %in% x) %>%
        summarise(sum(N)) %>%
        as.numeric()
})
# Distribution N obs per chunk
sapply(reqs, function(x){
    my_keys %>%
        filter(db_code %in% x) %>%
        summarise(sum(N)) %>%
        as.numeric()
}) %>% hist()
# Avg N obs per chunk
sapply(reqs, function(x){
    my_keys %>%
        filter(db_code %in% x) %>%
        summarise(sum(N)) %>%
        as.numeric()
}) %>% mean()
# max N obs per chunk
sapply(reqs, function(x){
    my_keys %>%
        filter(db_code %in% x) %>%
        summarise(sum(N)) %>%
        as.numeric()
}) %>% max()
# min N obs per chunk
sapply(reqs, function(x){
    my_keys %>%
        filter(db_code %in% x) %>%
        summarise(sum(N)) %>%
        as.numeric()
}) %>% min()


# Save my query arguments ----
query <- data.frame(predicate = c("taxonKey", 
                                  "basisOfRecord", 
                                  "year", 
                                  "year", 
                                  "hasCoordinate", 
                                  "hasGeospatialIssue", 
                                  "occurrenceStatus", 
                                  "creator", 
                                  "notification_address"),
                    type = c("pred_in", 
                             "pred_in", 
                             "pred_gte", 
                             "pred_lte", 
                             "pred", 
                             "pred", 
                             "pred", 
                             "pred", 
                             "pred"),
                    value = c(NA, 
                              paste(basisOfRecord, collapse = ", "), 
                              range_env_data[1], 
                              range_env_data[2], 
                              TRUE, 
                              FALSE, 
                              "PRESENT", 
                              "Brunno F. Oliveira", 
                              "brunno.oliveira@me.com"))

write.csv(query, "Data/my_query_to_GBIF.csv")

query <- read.csv("Data/my_query_to_GBIF.csv")

## Submit the requests ----


# Loop through chunks
requests <- data.frame()
k = 1
inforequest <- list()

for(i in 1:length(reqs)){ 
    
    cat("\r", "Requesting chunk", i, "from", length(reqs))
    
    if(k <= 3){ # GBIF allow 3 requests at a time
        inforequest[[k]] <- occ_download(
            # pred("taxonKey", 8417931),
            pred_in("speciesKey", reqs[[i]]),
            pred_in("basisOfRecord",
                    basisOfRecord),
            pred_gte("year", range_env_data[1]),
            pred_lte("year", range_env_data[2]),
            pred("hasCoordinate", TRUE),
            pred("hasGeospatialIssue", FALSE),
            pred("occurrenceStatus", "PRESENT"),
            format = "SIMPLE_CSV",
            user = Sys.getenv('GBIF_USER'), 
            pwd = Sys.getenv('GBIF_PWD'),
            email = "brunno.oliveira@me.com")
        
        tmp <- inforequest[[k]]
        
        # save information about request
        tmp.save <- attributes(tmp)
        tmp.save <- data.frame(download_key = tmp[1],
                               created = tmp.save$created,
                               download_link = tmp.save$downloadLink,
                               doi = tmp.save$doi,
                               citation = tmp.save$citation,
                               format = tmp.save$format,
                               user = tmp.save$user,
                               email = tmp.save$email)
        requests <- rbind(requests,tmp.save)
        
        k = k + 1 
        
    }
    
    if(k > 3){ # GBIF allows 3 requests at a time
        occ_download_wait(inforequest[[1]])
        occ_download_wait(inforequest[[2]])
        occ_download_wait(inforequest[[3]])
        
        k = 1
        inforequest <- list()
    }
}

write.csv(requests, "Data/requests.csv", row.names = F)



# Download requests at the HPC ----


requests <- read.csv("Data/requests.csv")
head(requests)

GBIF_zip_dir <- here::here(scratch_dir,"Data/GBIF")
if(!dir.exists(GBIF_zip_dir)){
    dir.create(GBIF_zip_dir,recursive = T)
}

ncores = parallelly::availableCores()

check <- "Error"
attempt <- 0
errors <- TRUE

keystogo <- requests$download_key

while( any(errors) ) {
    attempt <- attempt + 1
    
    cat("Attempt", attempt, "\n")
    
    check <- mclapply(1:length(keystogo), function(i){
        
        keystogo_i = as.character(keystogo[i])
        
        test = occ_download_meta(keystogo_i)
        test = test$status
        
        while(test == "RUNNING"){
            
            test = occ_download_meta(keystogo_i)
            test = test$status
            
            Sys.sleep(60)
            
        }
        
        check[[i]] <- try (
            {
                rgbif::occ_download_get(key = keystogo_i,
                                        path = GBIF_zip_dir,
                                        overwrite = TRUE)
                return(paste("OK", i))
            },
            silent = TRUE
        )
        
    }
    , mc.cores = ifelse(ncores>length(keystogo),length(keystogo),ncores))
    
    errors <- sapply(check, function(x) class(x)=="try-error")
    
    if(any(errors)){
        errors <- requests$download_key[which(errors)]
        keystogo <- errors
    }
    cat("Error on keys:", errors)
} 


# Test if everything downloaded correctly
# Compare file sizes of remote and local files
download_size <- function(url) as.numeric(httr::HEAD(url)$headers$`content-length`)

for(i in 1:length(keystogo)){ cat("\rChecking file", i, "from", length(keystogo))
    # local file size
    size_local <- file.size(here::here(GBIF_zip_dir, paste0(keystogo[i],".zip")))
    # remote file size
    size_remote <- download_size(requests$download_link[i])
    # download again if files dont have the same size
    if(!size_local == size_remote){
        cat("\rFiles don't have the same size\nDownloading file", i)
        rgbif::occ_download_get(key = keystogo[i],
                                path = GBIF_zip_dir,
                                overwrite = TRUE)
    }
}





# Filter and save species occurrences ----
GBIF_zip_dir <- here::here(scratch_dir,"Data/GBIF/")
if(!dir.exists(GBIF_zip_dir)){
    dir.create(GBIF_zip_dir,recursive = T)
}

myselection <- c("basisOfRecord", "speciesKey", "decimalLongitude", "decimalLatitude", "year", "month","species")
basisOfRecord_selected <- basisOfRecord

terrestrials <- N_OCC$scientificName[which(N_OCC$Eco == "Ter")]
marines <- N_OCC$scientificName[which(N_OCC$Eco == "Mar")]

# my terrestrial raster
rasdir <- "/storage/simple/projects/t_cesab/brunno/Exposure-SDM/Data"
ter.ras <- terra::rast(here::here(vars_dir("Ter"),"model_raster_ter_1km.tif"))
# my marine raster
mar.ras <- terra::rast(here::here(vars_dir("Mar"),"model_raster_mar.tif"))

# test
any(terrestrials %in% marines)

# zipped occ files
occs <- list.files(here::here(GBIF_zip_dir), pattern = ".zip")

# create dir to save sps occurrences
if(!dir.exists(occ_dir)){
    dir.create(occ_dir, recursive = T)
}

# create temp dir to decompress zip file
tmp.dir <- here::here(occ_dir,"tmp")
if(!dir.exists(tmp.dir)){
    dir.create(tmp.dir)
}


ncores = parallelly::availableCores()


test_if_work <- mclapply(1:length(occs), function(i){
    
    
    # decompress zip file
    zippedfile <- here::here(GBIF_zip_dir,occs[i])
    decompress_file(directory = tmp.dir, file = zippedfile)
    unzippedfile <- gsub(".zip",".csv",occs[i])
    unzippedfile <- gsub(".zip",".csv",here::here(tmp.dir,unzippedfile))
    
    # read in
    tmp <- fread(unzippedfile, select = myselection, nThread = 10)
    tmp <- na.omit(tmp) 
    
    # filter basisOfRecord
    tmp <- tmp %>% dplyr::filter(basisOfRecord %in% basisOfRecord_selected)
    tmp_species <- unique(tmp$species)
    
    # remove species I already have data
    # saved species
    saved_sps <- list.files(occ_dir,pattern = ".qs")
    saved_sps <- gsub(".qs","",saved_sps)
    saved_sps <- gsub("_"," ",saved_sps)
    
    if(any(tmp_species %in% saved_sps)){
        tmp <- tmp %>% dplyr::filter(!species %in% saved_sps)
        tmp_species <- unique(tmp$species)
    }
    
    if(length(tmp_species)>0){
        
        # subset terrestrials
        ter. <- tmp[which(tmp$species %in% terrestrials),]
        # get cells
        cells. <-  terra::extract(ter.ras,
                                  ter.[,c("decimalLongitude", "decimalLatitude")],
                                  cells=TRUE)
        ter.$cell <- cells.$cell
        ter.$layer <- cells.[,2]
        # remove cells falling in the ocean
        ter. <- ter.[which(ter.$layer==1),] # keep land / remove ocean
        ter. = ter.[,-"layer"]
        # remove occurrences outside the temporal range of the env data
        ter. <- ter. %>% filter(year >= temporal_range_env_data("Ter")[1]+1 + n_yr_bioclimatic # because we cannot calculate bioclimatics for the year 1
                                & year <= temporal_range_env_data("Ter")[2])
        
        # subset marines
        mar. <- tmp[which(tmp$species %in% marines),] 
        # get cells
        cells. <-  terra::extract(mar.ras, 
                                  mar.[,c("decimalLongitude", "decimalLatitude")],
                                  cells=TRUE)
        mar.$cell <- cells.$cell
        mar.$layer <- cells.[,2]
        # remove cells falling in land
        mar. <- mar.[which(mar.$layer==1),] # remove land / keep ocean
        mar. = mar.[,-"layer"]
        # remove occurrences outside the temporal range of the env data
        mar. <- mar. %>% filter(year > temporal_range_env_data("Mar")[1] + n_yr_bioclimatic # because we cannot calculate bioclimatics for the year 1
                                & year <= temporal_range_env_data("Mar")[2])
        
        # Group all
        sps. <- rbind(ter.,mar.)
        # remove temporal ter. and mar. data to save memory space
        rm(ter., mar.);gc()
        
        # Remove duplicates: same species in the same location and date
        rm <- duplicated(sps.[,c("scientificName", "cell", "year", "month")])
        if(any(rm)){
            sps. <- sps.[-which(rm),]
        }
        
        # Remove other potential issues with the package CoordinateCleaner
        sps. = CoordinateCleaner::clean_coordinates(
            x = sps.,
            lon = "decimalLongitude",
            lat = "decimalLatitude",
            tests = c("capitals", "centroids", "equal",
                      "gbif", "institutions", "zeros"),
            value = "clean")
        
        # keep species with more than 30 obs
        my_sps_i <- table(sps.$species)
        my_sps_i <- names(my_sps_i[which(my_sps_i >= 30)])
        sps. <- sps.[which(sps.$species %in% my_sps_i),]
        
        # save data per species
        if(length(my_sps_i)>1){
            for(j in 1:length(my_sps_i)){ 
                cat("\rsaving sps", j, "from", length(my_sps_i))
                tmp_sps <- subset(sps., species == my_sps_i[j])
                
                qs::qsave(tmp_sps, 
                          here::here(occ_dir, paste0(gsub(" ","_",my_sps_i[j]),".qs")))
            }
        }
        
        # delete unzipped file
        unlink(unzippedfile)
        
    }
    
    
}
, mc.cores = ncores)

# delete tmp dir used to decompress zipfiles
unlink(tmp.dir, recursive = TRUE)
unlink(GBIF_zip_dir, recursive = TRUE)


#####################
# how many species?

all_sps <- list.files(here::here(occ_dir), pattern = '.qs')
all_sps <- gsub("_"," ",all_sps)
all_sps <- gsub(".qs","",all_sps)

length(all_sps)
# 8515

N_OCC <- read.csv("Data/N_OCC.csv")
nrow(N_OCC)


all_sps[which(!all_sps %in% N_OCC$scientificName)]

N_OCC <- N_OCC %>% filter(scientificName %in% all_sps)

# how many clean occurrences?
all_sps_c <- list.files(here::here(occ_dir), pattern = '.qs')
all_sps_c <- pbsapply(all_sps_c, function(x){
    # read in
    sp_i <- qs::qread(here::here(occ_dir,x))
    # count
    nrow(sp_i)
})
# N occurrence per species
summary(all_sps_c)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 30     373    2441   39063   15682 3092867

# how many marine species?
length(which(marines %in% all_sps))
# 385

# how many terrestrial species?
length(which(terrestrials %in% all_sps))
# 6924

# original sp list after removing species with < 30 occs from GBIF
nrow(N_OCC)
# 10895

# % bioshift with occ

splist <- read.csv("Data/Bioshifts/splist.csv", header = T)
# remove duplicated sp_names
splist <- splist %>%
    dplyr::select(scientificName,kingdom,phylum,class,order,family,db) %>%
    filter(!duplicated(scientificName))

bio <- data.frame(table(splist$class))
names(bio) <- c("Class","Bioshifts_Freq")

gbifocc <- data.frame(table(N_OCC$class))
names(gbifocc) <- c("Class","GBIF_Freq")

bio <- merge(bio, gbifocc)
bio$percent = paste0(round((bio$GBIF_Freq/bio$Bioshifts_Freq)*100,2)," % (N = ", bio$GBIF_Freq,")")

melted <- reshape::melt(bio[,c(1:3)], id="Class")
melted <- merge(melted, bio[,c(1,4)])
melted$percent[which(melted$variable == "GBIF_Freq")] = NA

keep <- unique(melted$Class[which(melted$value > 20)])
toplot = melted[which(melted$Class %in% keep),]
toplot$Class <- factor(toplot$Class, levels = unique(toplot$Class[order(toplot$value)]))

ggplot(toplot, 
       aes(x = Class, y = value, fill = variable))+
    geom_bar(stat="identity",position = "identity", alpha = .5) +
    geom_text(aes(label=percent),
              stat='identity', size = 3, hjust = "inward")+
    coord_flip()+
    ylab("N species")+
    theme_classic()+
    theme(legend.position="bottom",
          legend.title = element_blank()) 

