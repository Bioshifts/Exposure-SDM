
# Supplementary analyses 
library(here)
library(dplyr)
library(readr)
library(stringr)

# set computer
computer = "matrics"

if(computer == "muse"){
    setwd("/storage/simple/projects/t_cesab/brunno/Exposure-SDM")
}
if(computer == "matrics"){
    setwd("/users/boliveira/Exposure-SDM")
}
work_dir <- getwd()

# Get the list of species used in analyses
bioshifts <- read.csv(here(work_dir,"Data/Output/Bioshifts_merge_Exposure_all.csv"))
# fixing Locomotion_mode
bioshifts$Locomotion_mode = as.character(str_to_sentence(bioshifts$Locomotion_mode))

pos <- grep("Flying",bioshifts$Locomotion_mode)
bioshifts$Locomotion_mode[pos] <- "Fly"
pos <- grep("wimming",bioshifts$Locomotion_mode)
bioshifts$Locomotion_mode[pos] <- "Swim"
pos <- grep("alk",bioshifts$Locomotion_mode)
bioshifts$Locomotion_mode[pos] <- "Run/Walk"
pos <- grep("unning",bioshifts$Locomotion_mode)
bioshifts$Locomotion_mode[pos] <- "Run/Walk"
pos <- grep("lithering",bioshifts$Locomotion_mode)
bioshifts$Locomotion_mode[pos] <- "Slither/Crawl"
pos <- grep("rawling",bioshifts$Locomotion_mode)
bioshifts$Locomotion_mode[pos] <- "Slither/Crawl"
pos <- grep("Burrowing",bioshifts$Locomotion_mode)
bioshifts$Locomotion_mode[pos] <- "Sessile"

bioshifts$Locomotion_mode <- factor(
    bioshifts$Locomotion_mode,
    levels = c("Fly", "Run/Walk", "Swim", "Slither/Crawl", "Planktonic", "Sessile"))

# fixing Category
table(bioshifts$Category)
# remove synthesis because there is a only a single study
bioshifts <- bioshifts |> filter(!Category == "Synthesis")
# combine census-resurvey and survey-resurvey 
pos <- grep("Census-Resurvey",bioshifts$Category)
bioshifts$Category[pos] <- "Survey-Resurvey"
table(bioshifts$Category)
bioshifts$Category <- factor(bioshifts$Category, 
                             levels = c("TimeSeries", "CensusPeriods", "Survey-Resurvey"))

# make numeric as.numeric
pos <- which(sapply(bioshifts, is.numeric))
tmp <- lapply(bioshifts[pos], as.numeric)
tmp <- do.call(cbind,tmp)
bioshifts[pos] <- tmp

# remove one outlier species
bioshifts <- bioshifts |>
    filter(!(Eco == "Ter" & bvel_lat < -15))

# remove swimming terrestrials
bioshifts <- bioshifts |>
    filter(!(Eco == "Ter" & Locomotion_mode == "Swim"))

# select classes with more than X observations per edge
N_obs_class = 5
class_param_select <- bioshifts |>
    group_by(class, Param) |>
    tally() |>
    filter(n >= N_obs_class) |>
    mutate(class_param = paste(class,Param))

bioshifts <- bioshifts |> 
    mutate(class_param = paste(class,Param)) |> 
    filter(class_param %in% class_param_select$class_param) |>
    select(!class_param)

length(unique(bioshifts$sp_name_std))

# Select 10% of the terrestrial and marine species with >10.000 occurrence records
N_OCC <- read_csv("Data/sdms_occ.csv")
N_OCC <- N_OCC %>%
    dplyr::filter(Species %in% bioshifts$sp_name_std)
length(unique(N_OCC$Species))

table(N_OCC$ECO)

# Sample 10% of species of supplementary analyses

# Avoid species located in very large study areas
acceptable_sps <- bioshifts %>%
    dplyr::filter(Areakm2 <= 1848511.7) %>%
    dplyr::select(sp_name_std) %>%
    pull() %>%
    unique()
    
set.seed(123)  # for reproducibility

# 1) Compute original group sizes
eco_sizes <- N_OCC %>%
    dplyr::filter(ECO %in% c("Mar", "Ter")) %>%
    count(ECO, name = "N_total")

# 2) Filter eligible species and sample based on original sizes
N_OCC_sub <- N_OCC %>%
    filter(N_occ < 10000, # select species < 10.000 occurrences
           Species %in% acceptable_sps) %>% # Filter out very large areas
    left_join(eco_sizes, by = "ECO") %>%
    group_by(ECO) %>%
    group_modify(~ {
        n_target <- ceiling(0.1 * unique(.x$N_total))
        slice_sample(.x, n = min(nrow(.x), n_target))
    }) %>%
    ungroup()

table(N_OCC_sub$ECO)

write.csv(N_OCC_sub,
          "Data/Output/sps_sub_supp_analyses.csv",
          row.names = FALSE)


