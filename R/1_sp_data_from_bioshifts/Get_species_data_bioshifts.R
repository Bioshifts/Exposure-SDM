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



# Load parameters


# set computer
computer = "muse"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
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
    
    

# Load v3
biov3 <- read.csv("Data/Bioshifts/bioshifts_v3_raw.csv")
biov3 <- biov3 %>%
    filter(!is.na(sp_name)) %>%
    mutate(sp_name=gsub("_"," ",sp_name))
# add taxonomy
tojoin <- splist %>%
    dplyr::select(scientificName, kingdom, phylum, class, order, family) %>%
    mutate(sp_name = scientificName) %>%
    filter(!duplicated(sp_name))

any(!biov3$sp_name %in% tojoin$sp_name)

biov3 <- merge(biov3, tojoin, by = "sp_name", all.x = TRUE)

# Load v1
biov1 <- read.csv("Data/Bioshifts/biov1_fixednames.csv", header = T)
# Fix references in biov1
biov1$ID <- gsub("[^0-9.-]", "", biov1$Article_ID)
biov1$ID <- as.numeric(biov1$ID)

biov1$sp_name_std_v1 <- gsub("_"," ",biov1$sp_name_std_v1)

biov1 <- biov1 %>%
    mutate(
        Type = case_when(
            Type=="HOR" ~ "LAT",
            TRUE ~ as.character(Type)),
        Species = sp_name_std_v1, 
        # Standardize shift measures to km/year
        SHIFT = case_when(
            UNIT == "m/year" ~ SHIFT/1000, 
            UNIT == "km/year" ~ SHIFT, 
            TRUE ~ NA_real_),
        Trend_unit = Unit.trend,
        # Standardize trend to km 
        Trend.std = SHIFT*DUR,
        Trend.std = case_when(
            Type == "LAT" ~ round_up(Trend.std,3), 
            Type == "ELE" ~ round_up(Trend.std,6), 
            TRUE ~ NA_real_),
        Azi = ifelse(is.na(Azimuth),0,1))

biov1 <- biov1 %>%
    filter(!is.na(sp_name_std_v1))

all(biov1$sp_name_std_v1 %in% splist$scientificName[which(splist$v1==1)])


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
        ID=Paper.ID,
        Study_ID=Study.Period,
        Type=Type,
        Param=Parameter,
        Trend=Numeric.Change,
        Trend_unit=Unit,
        DUR=as.numeric(Total.Years.Studied),
        START=Start.Year,
        END=End.Year,
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
    mutate(
        Type = case_when(
            Dimension=="latitude" ~ "LAT",
            Dimension=="elevation" ~ "ELE",
            TRUE ~ as.character(Dimension)),
        Param = case_when(
            Parameter=="leading edge" ~ "LE",
            Parameter=="maximum/optimum" ~ "O",
            Parameter=="trailing edge" ~ "TE",
            Parameter=="mean" ~ "O",
            TRUE ~ as.character(Parameter)),
        # Standardize shift measures to km/year
        SHIFT = case_when( 
            Trend_unit == "ft" ~ Trend/3280.8 /Total.Years.Studied,
            Trend_unit == 'degrees lat' ~ Trend*111/Total.Years.Studied,
            Trend_unit =='degrees lat/year' ~ Trend*111,
            Trend_unit=='degree latitude per decade'~ Trend*111/10,
            Trend_unit=='km' ~ Trend/Total.Years.Studied,
            Trend_unit=='km/decade' ~ Trend/10,
            Trend_unit=='km/year' ~ Trend,
            Trend_unit=='m' ~ Trend/1000/Total.Years.Studied,
            Trend_unit=='m/year' ~ Trend/1000,
            Trend_unit=='m/decade'~ Trend/1000/10,
            TRUE ~ NA_real_),
        # Standardize trend to km 
        Trend.std = SHIFT*DUR,
        Trend.std = case_when(
            Type == "LAT" ~ round_up(Trend.std,3), 
            Type == "ELE" ~ round_up(Trend.std,6), 
            TRUE ~ NA_real_))

biov2 <- biov2 %>%
    filter(!is.na(sp_name_std_v2))

all(biov2$biov2 %in% splist$scientificName[which(splist$v2==1)])



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


if(any(grep("sp[.]",biov3$sp_reported_name_v3))){
    biov3 <- biov3 %>% filter(!grepl("sp[.]",sp_reported_name_v3))
}
if(any(grep("sp[.]",biov3$sp_name_std_v3))){
    biov3 <- biov3 %>% filter(!grepl("sp[.]",sp_name_std_v3))
}
if(any(grep("cf[.]",biov3$sp_reported_name_v3))){
    biov3 <- biov3 %>% filter(!grepl("cf[.]",sp_reported_name_v3))
}
if(any(grep("cf[.]",biov3$sp_name_std_v3))){
    biov3 <- biov3 %>% filter(!grepl("cf[.]",sp_name_std_v3))
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
MFishv2 <- unique(biov2$sp_name_std_v2[biov2$group == "Fish" & biov2$ECO == "marine"])
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

biov3$ECO = NA
biov3$ECO[which(biov3$sp_name_std_v3 %in% Terrestrials)] <- "T"
biov3$ECO[which(biov3$sp_name_std_v3 %in% Marine)] <- "M"
biov3$ECO[which(biov3$sp_name_std_v3 %in% Aquatic)] <- "A"
biov3$ECO[which(biov3$sp_name_std_v3 %in% MFish)] <- "M"
biov3$ECO[which(biov3$sp_name_std_v3 %in% FFish)] <- "A"

splist$Group = NA
splist$Group[which(splist$class == "Phaeophyceae")] <- "Chromista"
splist$kingdom[which(splist$class == "Phaeophyceae")] <- "Chromista"
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

biov3 <- merge(biov3[,-which(names(biov3) %in% c("Group","ECO"))], splist[,c("Group","ECO","scientificName")], 
               by.x = "sp_name", by.y = "scientificName", 
               all.x = T)

table(biov1$ECO)
table(biov2$ECO)
table(biov3$ECO)
table(splist$ECO)

all(biov1$sp_name_std_v1 %in% biov3$scientificName)



# *Populate v3 columns*



# Trend
biov3$Trend <- as.numeric(biov3$Trend_v1)
pos <- which(is.na(biov3$Trend))
biov3$Trend[pos] <- biov3$Trend_v2[pos]

# Unit trend
biov3$Unit.trend <- biov3$Unit.trend_v1
pos <- which(is.na(biov3$Unit.trend))
biov3$Unit.trend[pos] <- biov3$Unit.trend_v2[pos]

# Start
biov3$START <- biov3$START_v1
pos <- which(is.na(biov3$START))
biov3$START[pos] <- biov3$START_v2[pos]

# End
biov3$END <- biov3$END_v1
pos <- which(is.na(biov3$END))
biov3$END[pos] <- biov3$END_v2[pos]

# Duration
biov3$DUR <- biov3$END-biov3$START+1

# Trend km (standardized)
unique(biov3$Unit.trend)
biov3 <- biov3 %>%
    mutate(
        Trend_km = as.numeric(Trend),
        DUR = as.numeric(DUR),
        Trend_km = case_when(
            Unit.trend == "degree" ~ Trend_km*111,
            Unit.trend == "degrees lat" ~ Trend_km*111,
            Unit.trend =='degree/year' ~ (Trend_km*111)*DUR,
            Unit.trend=='degrees lat/year'~ (Trend_km*111)*DUR,
            Unit.trend=='degree/decade'~ (Trend_km*111/10)*DUR,
            Unit.trend=='degree latitude per decade'~ (Trend_km*111/10)*DUR,
            Unit.trend=='km' ~ Trend_km,
            Unit.trend=='km/decade' ~ (Trend_km/10)*DUR,
            Unit.trend=='km/year' ~ Trend_km*DUR,
            Unit.trend=='m' ~ Trend_km/1000,
            Unit.trend=='m/year' ~ (Trend_km/1000)*DUR,
            Unit.trend=='m/decade'~ (Trend_km/1000/10)*DUR
        )
    )

# SHIFT (shift values come standardized)
biov3$SHIFT <- biov3$SHIFT_v1
pos <- which(is.na(biov3$SHIFT))
biov3$SHIFT[pos] <- biov3$SHIFT_v2[pos]

# remove categorical shifts
biov3 <- biov3 %>%
    mutate(SHIFT=as.numeric(SHIFT)) %>%
    filter(!is.na(SHIFT))



## Filter LAT / ELE shifts


biov3 <- biov3 %>%
    filter(Type %in% c("ELE","LAT"))  # Use LAT ELE shifts

# splist
splist <- splist %>% filter(scientificName %in% biov3$sp_name)



## Filter fresh water fish and marine birds


#remove freshwater fishes
biov3 <- biov3[-which((biov3$class == "Actinopterygii" | biov3$Group == "Fish") & biov3$ECO=="A"),]

#remove marine birds
biov3 <- biov3[-which(biov3$class == "Aves" & biov3$ECO=="M"),]



# Subset

## Filter species with GBIF data


# splist
splist <- splist[which(splist$db == "gbif"),]

biov3 <- biov3 %>% filter(sp_name %in% splist$scientificName)





## Filter shifts within the temporal range of the environmental data


biov3 <- biov3 %>%
    dplyr::filter(
        # For terrestrials
        ((ECO == "T" | ECO == "A") & 
             # start > beginning env data
             (START >= get_temporal_range_env_data("Ter")[1])) |
            # For marine
            ((ECO == "M") & 
                 # start > beginning env data
                 (START >= get_temporal_range_env_data("Mar")[1]))
    )

# splist
bios <- unique(biov3$sp_name)
splist <- splist %>% filter(scientificName %in% bios)

table(splist$ECO)

table(biov3$ECO)
table(biov1$ECO)
table(biov2$ECO)




# N occurrences GBIF

# Retrieve a summary table with N occurrence per species.
# We only retrieve N occurrences for records based:
#     1) "HUMAN_OBSERVATION", "LITERATURE", "LIVING_SPECIMEN", "OBSERVATION", "PRESERVED_SPECIMEN", "OCCURRENCE"
# 2) from which there was coordinates recorded
# 3) Registered after 1979 (limit is 2022)


ncores = as.numeric(Sys.getenv("SLURM_CPUS_ON_NODE"))

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

# Took 4min FRB
# Took 3min MESO

length(N_OCC) == length(splist$scientificName)

N_OCC <- data.frame(splist, n_occ = N_OCC)

write.csv(N_OCC, "Data/n_occ.csv", row.names = F)

N_OCC <- read.csv("Data/n_occ.csv")

length(unique(N_OCC$scientificName[which(N_OCC$n_occ>30)]))

# N species
tmp <- N_OCC %>%
    filter(n_occ > 30) 
length(unique(tmp$scientificName))
sort(table(tmp$Group))
tmp %>%
    group_by(ECO, Group) %>%
    tally() 
sort(table(tmp$ECO))

# N shifts
biov3 %>%
    filter(sp_name %in% tmp$scientificName) %>%
    nrow()
tmp2 <- biov3 %>%
    filter(sp_name %in% tmp$scientificName) 
sort(table(tmp2$Group))
sort(table(tmp2$ECO))



# Stats




## N occurrences by group

toplot <- N_OCC %>%
    dplyr::filter(n_occ > 30) %>%
    group_by(kingdom, Group) %>%
    dplyr::summarize(n_occ = sum(n_occ))

toplot$class <- factor(toplot$Group, levels = toplot$Group[order(toplot$n_occ)])

# Occ by Class
ggplot(toplot, aes(x = Group, y = n_occ))+
    geom_col()+
    theme_classic()+
    # coord_flip()+
    ylab("N occurrences")+xlab("Class")+
    theme(axis.text.x = element_text(angle = 45, hjust=1))



## N species by group

toplot <- N_OCC %>%
    filter(n_occ > 30) %>%
    group_by(kingdom, Group) %>%
    dplyr::summarise(n_occ = length(unique(scientificName)))

toplot$class <- factor(toplot$Group, levels = toplot$Group[order(toplot$n_occ)])

# Occ by Class
ggplot(toplot, aes(x = Group, y = n_occ))+
    geom_col()+
    theme_classic()+
    # coord_flip()+
    ylab("N species")+xlab("Class")+
    theme(axis.text.x = element_text(angle = 45, hjust=1))



## N species by ECO

toplot <- N_OCC %>%
    group_by(ECO, Group) %>%
    dplyr::summarise(n_sps = length(unique(scientificName)))
toplot$n_sps <- as.numeric(toplot$n_sps)

# Occ by ECO
ggplot(toplot, aes(x = ECO, y = n_sps))+
    geom_col(fill = "white", color = "black")+
    theme_classic()+
    ylab("N species")+xlab("Habitat")+
    facet_grid(scales = "free", space = "free", cols = vars(Group))+
    theme(strip.text.x = element_text(angle = 90))



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



