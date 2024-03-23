# title: "Get species data for SDMs"
# author: "Brunno F Oliveira & Bioshifts group"


# Setup

rm(list=ls())
gc()

list.of.packages <- c("stringr", "rgbif", "parallel", "pbapply", "raster", "dplyr", "here", "tidyr",
                      "data.table", "ggplot2", "readxl", "knitr","googledrive","googleAuthR",
                      "sqldf","RSQLite", "scales")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)



########################
# set computer
computer = "personal"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
    work_dir <- getwd()
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
    work_dir <- getwd()
}
if(computer == "personal"){
    work_dir <- getwd()
}

# Load functions

source("R/my_functions.R")
source("R/round-up.R")
source("R/settings.R")

range_env_data <- range(c(temporal_range_env_data("Ter"),temporal_range_env_data("Mar")))


# *Load species list*

splist <- read.csv("Data/Bioshifts/splist.csv", header = T)
# remove duplicated sp_names
splist <- splist %>%
    dplyr::select(scientificName,kingdom,phylum,class,order,family,db,db_code) %>%
    filter(!duplicated(scientificName))



# Load raw bioshifts
# *While we dont have a v3, just add papers from v2 to v1*
    
# Load v1
biov1 <- read.csv("Data/Bioshifts/biov1_fixednames.csv", header = T)

biov1$sp_name_std_v1 <- gsub("_"," ",biov1$sp_name_std_v1)

biov1 <- biov1 %>%
    mutate(
        Type = case_when(
            Type=="HOR" ~ "LAT",
            TRUE ~ as.character(Type)),
        Species = sp_name_std_v1) %>%
    filter(Type == "LAT", # Use only latitudinal shifts
           (ECO == "T" | ECO == "M")) # Shifts Marine or Terrestrial

biov1 <- biov1 %>%
    filter(!is.na(sp_name_std_v1))

all(biov1$sp_name_std_v1 %in% splist$scientificName)


# Load v2 
biov2 <- read.csv("Data/Bioshifts/biov2_fixednames.csv", header = T)
biov2$Scientific.Name <- gsub(" ","_",biov2$Scientific.Name)
biov2$sp_name_std_v2 <- gsub("_"," ",biov2$sp_name_std_v2)

biov2 <- biov2  %>%
    mutate(
        Type = case_when(
            Dimension=="latitude" ~ "LAT",
            Dimension=="elevation" ~ "ELE",
            TRUE ~ as.character(Dimension)),
        group=Taxonomic.Group,
        ECO=case_when(
            Ecosystem.Type=="terrestrial" ~ "T",
            Ecosystem.Type=="marine" ~ "M",
            Ecosystem.Type=="aquatic" ~ "A",
            TRUE ~ as.character(Ecosystem.Type)),
        Phylum=phylum,
        Class=class,
        Order=order,
        Family=family,
        Genus=genus,
        Species=sp_name_std_v2) %>%
    filter(Type == "LAT", # Use only latitudinal shifts
           (ECO == "T" | ECO == "M")) # Shifts Marine or Terrestrial

biov2 <- biov2 %>%
    filter(!is.na(sp_name_std_v2))

all(biov2$biov2 %in% splist$scientificName)

# filter sp list
splist <- splist %>% filter((scientificName %in% biov1$sp_name_std_v1) | (scientificName %in% biov2$sp_name_std_v2))

## *Remove species identified to the genus level or cf.*



if(any(grep("sp[.]",biov1$sp_reported_name_v1))){
    biov1 <- biov1 %>% filter(!grepl("sp[.]",sp_reported_name_v1))
}
if(any(grep("sp[.]",biov1$sp_name_std_v1))){
    biov1 <- biov1 %>% filter(!grepl("sp[.]",sp_name_std_v1))
}
if(any(grep("cf[.]",biov1$sp_reported_name_v1))){
    biov1 <- biov1 %>% filter(!grepl("cf[.]",sp_reported_name_v1))
}
if(any(grep("cf[.]",biov1$sp_name_std_v1))){
    biov1 <- biov1 %>% filter(!grepl("cf[.]",sp_name_std_v1))
}

if(any(grep("sp[.]",biov2$sp_reported_name_v2))){
    biov2 <- biov2 %>% filter(!grepl("sp[.]",sp_reported_name_v2))
}
if(any(grep("sp[.]",biov2$sp_name_std_v2))){
    biov2 <- biov2 %>% filter(!grepl("sp[.]",sp_name_std_v2))
}
if(any(grep("cf[.]",biov2$sp_reported_name_v2))){
    biov2 <- biov2 %>% filter(!grepl("cf[.]",sp_reported_name_v2))
}
if(any(grep("cf[.]",biov2$sp_name_std_v2))){
    biov2 <- biov2 %>% filter(!grepl("cf[.]",sp_name_std_v2))
}


if(any(grep("sp[.]",splist$scientificName))){
    splist <- splist %>% filter(!grepl("sp[.]",scientificName))
}
if(any(grep("sp[.]",splist$scientificName))){
    splist <- splist %>% filter(!grepl("sp[.]",scientificName))
}
if(any(grep("cf[.]",splist$scientificName))){
    splist <- splist %>% filter(!grepl("cf[.]",scientificName))
}
if(any(grep("cf[.]",splist$scientificName))){
    splist <- splist %>% filter(!grepl("cf[.]",scientificName))
}



## *Classify species*



# Class species as Terrestrial, Marine or Freshwater
Terv1 <- unique(biov1$sp_name_std_v1[which(biov1$ECO == "T")])
Terv2 <- unique(biov2$sp_name_std_v2[which(biov2$ECO == "T")])
Terrestrials = unique(c(Terv1,Terv2))

Mar1 <- unique(biov1$sp_name_std_v1[which(biov1$ECO == "M")])
Mar2 <- unique(biov2$sp_name_std_v2[which(biov2$ECO == "M")])
Marine = unique(c(Mar1,Mar2))

Aquatic = unique(biov2$sp_name_std_v2[which(biov2$ECO == "A")])

# Freshwater fish
FFishv1 <- unique(biov1$sp_name_std_v1[(biov1$Class == "Actinopterygii" | biov1$Class == "Cephalaspidomorphi") & biov1$ECO == "T"])
FFishv2 <- unique(biov2$sp_name_std_v2[biov2$group == "Fish" & 
                                           (biov2$ECO == "T" | biov2$ECO == "A")])
FFish = unique(c(FFishv1,FFishv2))

# Marine fish
MFishv1 <- unique(biov1$sp_name_std_v1[biov1$Class == "Actinopterygii" & biov1$ECO == "M"])
MFishv2 <- unique(biov2$sp_name_std_v2[biov2$group == "Fish" & biov2$ECO == "M"])
MFish = unique(c(MFishv1,MFishv2))

splist$ECO = NA
splist$ECO[which(splist$scientificName %in% Terrestrials)] <- "T"
splist$ECO[which(splist$scientificName %in% Marine)] <- "M"
splist$ECO[which(splist$scientificName %in% Aquatic)] <- "A"
splist$ECO[which(splist$scientificName %in% MFish)] <- "M"
splist$ECO[which(splist$scientificName %in% FFish)] <- "A"

biov1$ECO = NA
biov1$ECO[which(biov1$sp_name_std_v1 %in% Terrestrials)] <- "T"
biov1$ECO[which(biov1$sp_name_std_v1 %in% Marine)] <- "M"
biov1$ECO[which(biov1$sp_name_std_v1 %in% Aquatic)] <- "A"
biov1$ECO[which(biov1$sp_name_std_v1 %in% MFish)] <- "M"
biov1$ECO[which(biov1$sp_name_std_v1 %in% FFish)] <- "A"

biov2$ECO = NA
biov2$ECO[which(biov2$sp_name_std_v2 %in% Terrestrials)] <- "T"
biov2$ECO[which(biov2$sp_name_std_v2 %in% Marine)] <- "M"
biov2$ECO[which(biov2$sp_name_std_v2 %in% Aquatic)] <- "A"
biov2$ECO[which(biov2$sp_name_std_v2 %in% MFish)] <- "M"
biov2$ECO[which(biov2$sp_name_std_v2 %in% FFish)] <- "A"


splist$Group = NA
splist$kingdom[which(splist$class == "Phaeophyceae")] <- "Chromista"
splist$Group[which(splist$kingdom == "Chromista")] <- "Chromista"

splist$Group[which(splist$phylum == "Rhodophyta")] <- "Seaweed"
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
splist$Group[which(splist$class == "Hexanauplia")] <- "Crustacean"
splist$Group[which(splist$class == "Maxillopoda")] <- "Crustacean"
splist$Group[which(splist$class == "Ostracoda")] <- "Crustacean"
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
splist$Group[which(splist$class == "Chilopoda")] <- "Centipedes"
splist$Group[which(splist$class == "Diplopoda")] <- "Millipedes"
splist$Group[which(splist$class == "Amphibia")] <- "Amphibian"
splist$Group[which(splist$kingdom == "Fungi")] <- "Fungi"
splist$Group[which(splist$order == "Balanomorpha")] <- "Barnacles"
splist$Group[which(splist$phylum == "Nematoda")] <- "Nematodes"
splist$Group[which(splist$class == "Myxini")] <- "Hagfish"

######################################
biov1 <- merge(biov1[,-which(names(biov1) %in% c("Group","ECO"))], 
               splist[,c("Group","ECO","scientificName")], 
               by.x = "sp_name_std_v1", by.y = "scientificName", 
               all.x = T)

biov2 <- merge(biov2[,-which(names(biov2) %in% c("Group","ECO"))], splist[,c("Group","ECO","scientificName")], 
               by.x = "sp_name_std_v2", by.y = "scientificName", 
               all.x = T)

table(biov1$ECO)
table(biov2$ECO)
table(splist$ECO)

all(biov1$sp_name_std_v1 %in% splist$scientificName)

any(duplicated(splist$scientificName))


## Filter fresh water fish and marine birds


#remove freshwater fishes
splist <- splist[-which((splist$class == "Actinopterygii" | splist$Group == "Fish") & splist$ECO=="A"),]

#remove marine birds
splist <- splist[-which(splist$class == "Aves" & splist$ECO=="M"),]



# Subset

## Filter species with GBIF data


# splist
splist <- splist[which(splist$db == "gbif"),]



## Filter shifts within the temporal range of the environmental data
## work with v1 only as we dont have v3 yet
## remove fresh water shifts

biov1 <- biov1 %>%
    dplyr::filter(
        # For terrestrials
        ((ECO == "T") & 
             # start > beginning env data
             (START >= temporal_range_env_data("Ter")[1])) |
            # For marine
            ((ECO == "M") & 
                 # start > beginning env data
                 (START >= temporal_range_env_data("Mar")[1]))
    )

# splist
splist <- splist %>% filter(scientificName %in% unique(biov1$sp_name_std_v1))

table(splist$ECO)

table(biov1$ECO)



# N occurrences GBIF

# Retrieve a summary table with N occurrence per species.
# We only retrieve N occurrences for records based:
# 1) "HUMAN_OBSERVATION" "OBSERVATION" "OCCURRENCE"
# 2) from which there was coordinates recorded
# 3) Registered after 1901 


ncores = parallelly::availableCores()

cl <- makeCluster(ncores)
clusterExport(cl, c("splist","basisOfRecord","range_env_data"))

N_OCC <- pbsapply(1:length(splist$db_code), function(i) {
    
    code <- gsub("GBIF:","",splist$db_code[i])
    
    records <- sapply(basisOfRecord, function(b) {
        rgbif::occ_count(taxonKey = code,
                         hasGeospatialIssue = FALSE,
                         basisOfRecord = b,
                         from = range_env_data[1], 
                         to = range_env_data[2]) 
    })
    
    sum(records)
}
,
cl = cl)

stopCluster(cl)

# Took 9min FRB
# Took 3min MESO

length(N_OCC) == length(splist$scientificName)

N_OCC <- data.frame(splist, n_occ = N_OCC)

write.csv(N_OCC, "Data/n_occ.csv", row.names = F)

N_OCC <- read.csv("Data/n_occ.csv")


# N species
length(unique(N_OCC$scientificName))

tmp <- N_OCC %>%
    filter(n_occ > 30) 
length(unique(tmp$scientificName))

sort(table(tmp$Group))

tmp %>%
    group_by(ECO, Group) %>%
    tally() 
sort(table(tmp$ECO))

# N shifts
biov1 %>%
    filter(sp_name_std_v1 %in% tmp$scientificName) %>%
    nrow()
tmp2 <- splist %>%
    filter(scientificName %in% tmp$scientificName) 
sort(table(tmp2$Group))
sort(table(tmp2$ECO))



# Stats




## N occurrences by group

# Occ by Class
ggplot(N_OCC, aes(x = Group, y = n_occ))+
    geom_col()+
    theme_classic()+
    # coord_flip()+
    ylab("N occurrences")+xlab("Class")+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    facet_wrap(.~ECO, scales = "free")



## N species by group

toplot <- N_OCC %>%
    filter(n_occ > 30) %>%
    group_by(kingdom, Group, ECO) %>%
    dplyr::summarise(n_occ = length(unique(scientificName)))

# Occ by Class
ggplot(toplot, aes(x = Group, y = n_occ))+
    geom_col()+
    theme_classic()+
    # coord_flip()+
    ylab("N species")+xlab("Class")+
    theme(axis.text.x = element_text(angle = 45, hjust=1))+
    facet_wrap(.~ECO, scales = "free")


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
bio$percent = paste(round((bio$GBIF_Freq/bio$Bioshifts_Freq)*100,2),"%")

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



