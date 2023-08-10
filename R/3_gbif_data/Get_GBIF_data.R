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
                      "sqldf","RSQLite","jsonlite",
                      "meowR")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)






# set computer ----
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    work_dir <- getwd()
}



# Load functions ----

source("R/my_functions.R")

source("R/settings.R")

range_env_data <- range(c(temporal_range_env_data("Ter"),temporal_range_env_data("Mar")))



# Get species list ----



N_OCC <- read.csv("Data/n_occ.csv")

N_OCC <- N_OCC %>% filter(n_occ >= 30)



# GBIF requests ----

# Create requests for all taxa keys. The request is processed by GBIF. The processing of each request takes a while (GBIF website says it can take up to 15min). Once the request is processed, the file is ready for download under my user page at the GBIF website. The advantage of this method is reproducibility. It generates a DOI that can be cited and link for download that can be shared. Creating a single request can generate a very large file (File size was 300Gbs zipped == 600Gbs unziped). Thus, better submitting multiple requests. Here, we created divided species into chunks and loop through chunks for downloading data.

## Create queries



my_keys = N_OCC %>%
    group_by(db_code) %>%
    dplyr::summarise(N=sum(n_occ))

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


# # Find requests for these species
# sps_test = c("Eriochloa contracta", "Diplacus clevelandii", "Bloomeria clevelandii", "Astragalus mohavensis")
# ids_test <- N_OCC$db_code[which(N_OCC$scientificName %in% sps_test)]
# ids_test <- gsub("GBIF:","",ids_test)
# 
# reqs_sel <- sapply(reqs, function(x){
#     any(x %in% ids_test)
# })
# reqs_sel <- which(reqs_sel)

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
            pred_in("taxonKey", reqs[[i]]),
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

new.dir <- here::here(scratch_dir,"Data/GBIF")
if(!dir.exists(new.dir)){
    dir.create(new.dir,recursive = T)
}

ncores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))

check <- "Error"
attempt <- 0
errors <- TRUE

keystogo <- requests$download_key

while( any(errors) ) {
    attempt <- attempt + 1
    
    cat("Attempt", attempt, "\n")
    
    check <- mclapply(1:length(keystogo), function(i){
        
        check[[i]] <- try (
            {
                rgbif::occ_download_get(key = keystogo[i],
                                        path = new.dir,
                                        overwrite = TRUE)
                return(paste("OK", i))
            },
            silent = TRUE
        )
    }
    , mc.cores = ncores)
    
    errors <- sapply(check, function(x) class(x)=="try-error")
    
    if(any(errors)){
        errors <- requests$download_key[which(errors)]
        keystogo <- errors
    }
    cat("Error on keys:", errors)
} 


# Test if everythink downloaded correctly
# Compare file sizes of remote and local files
download_size <- function(url) as.numeric(httr::HEAD(url)$headers$`content-length`)

for(i in 1:length(keystogo)){ cat("\rChecking file", i, "from", length(keystogo))
    # local file size
    size_local <- file.size(here::here(new.dir, paste0(keystogo[i],".zip")))
    # remote file size
    size_remote <- download_size(requests$download_link[i])
    # download again if files dont have the same size
    if(!size_local == size_remote){
        cat("\rFiles don't have the same size\nDownloading file", i)
        rgbif::occ_download_get(key = keystogo[i],
                                path = new.dir,
                                overwrite = TRUE)
    }
}





# Filter and save species occurrences ----

myselection <- c("basisOfRecord", "speciesKey", "decimalLongitude", "decimalLatitude", "year", "month","species")
basisOfRecord_selected <- basisOfRecord

occ.dir <- here::here(scratch_dir,"Data/GBIF")

terrestrials <- N_OCC$scientificName[which(N_OCC$ECO == "T")]
marines <- N_OCC$scientificName[which(N_OCC$ECO == "M")]
aquatic <- N_OCC$scientificName[which(N_OCC$ECO == "A")]

# my terrestrial raster
rasdir <- "/storage/simple/projects/t_cesab/brunno/Exposure-SDM/Data"
ter.ras <- terra::rast(here::here(rasdir,"Land/cruts/model_raster_ter_1km.tif"))
# my marine raster
mar.ras <- terra::rast(here::here(rasdir,"Marine/model_raster_mar.tif"))
# my aquatic raster
# aqua.ras <- terra::rast(here::here(rasdir,"Aqua","model_raster_aqua.tif"))

# test
any(terrestrials %in% marines)
# test
any(aquatic %in% marines)

# zipped occ files
occs <- list.files(here::here(occ.dir), pattern = ".zip")

# create dir to save sps occurrences
sps.dir <- here::here(work_dir,"Data/GBIF_data")
if(!dir.exists(sps.dir)){
    dir.create(sps.dir, recursive = T)
}

# create temp dir to decompress zip file
tmp.dir <- here::here(occ.dir,"tmp")
if(!dir.exists(tmp.dir)){
    dir.create(tmp.dir)
}


ncores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))


test_if_work <- mclapply(1:length(occs), function(i){
    
    
    # decompress zip file
    zippedfile <- here::here(occ.dir,occs[i])
    decompress_file(directory = tmp.dir, file = zippedfile)
    unzippedfile <- gsub(".zip",".csv",occs[i])
    unzippedfile <- gsub(".zip",".csv",here::here(tmp.dir,unzippedfile))
    
    # read in
    tmp <- fread(unzippedfile, select = myselection, nThread = 1)
    tmp <- na.omit(tmp)
    
    # delete unzipped file
    unlink(unzippedfile)
    
    # filter basisOfRecord
    tmp <- tmp %>% dplyr::filter(basisOfRecord %in% basisOfRecord_selected)
    
    # saved species
    saved_sps <- list.files(sps.dir,pattern = ".qs")
    saved_sps <- gsub(".qs","",saved_sps)
    saved_sps <- gsub("_"," ",saved_sps)
    if(!any(tmp$species %in% saved_sps)){
        
        
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
        ter. <- ter. %>% filter(year >= temporal_range_env_data("Ter")[1] & year <= temporal_range_env_data("Ter")[2])
        
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
        mar. <- mar. %>% filter(year > temporal_range_env_data("Mar")[1])
        
        # Group all
        sps. <- rbind(ter.,mar.)
        
        # Remove duplicates: same species in the same location and date
        rm <- duplicated(sps.[,c("species", "cell", "year", "month")])
        if(any(rm)){
            sps. <- sps.[-which(rm),]
        }
        
        # Remove other potential issues
        sps. = CoordinateCleaner::clean_coordinates(x = sps.,
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
            for(j in 1:length(my_sps_i)){ cat("\rsaving sps", i, "from", length(my_sps_i))
                tmp_sps <- subset(sps., species == my_sps_i[j])
                
                qs::qsave(tmp_sps, here::here(sps.dir, paste0(gsub(" ","_",my_sps_i[j]),".qs")))
            }
        }
        
        
    }
    
    
}
, mc.cores = ncores)

# delete tmp dir used to decompress zipfiles
unlink(tmp.dir, recursive = TRUE)


# how many species?
all_sps <- list.files(here::here(sps.dir), pattern = '.qs')
all_sps <- gsub("_"," ",all_sps)
all_sps <- gsub(".qs","",all_sps)

length(all_sps)
# 7309

# how many clean occurrences?
all_sps_c <- list.files(here::here(sps.dir), pattern = '.qs')
all_sps_c <- pbsapply(all_sps_c, function(x){
    # read in
    sp_i <- qs::qread(here::here(sps.dir,x))
    # count
    nrow(sp_i)
})

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

N_OCC <- read.csv("Data/n_occ.csv")
N_OCC <- N_OCC %>% filter(scientificName %in% all_sps)


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

