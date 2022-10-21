
get_temporal_range_env_data <- function(realm){
    if(realm == "Ter"){c(1980, 2019)} else {c(1982, 2021)}
}
get_min_year_layers <- function(realm){
    if(realm == "Mar") {1982} else {1980}
}
get_myvars <- function(realm){
    if(realm == "Ter"){c("tasmax", "tasmin", "tas", "pr")} else {"SST"}
}
get_varsdir <- function(realm){
    if(realm == "Ter"){"/media/seagate/boliveira/Land"} else {"/media/seagate/boliveira/Marine"}
}
get_n_tiles <- function(realm){
    if(realm == "Ter"){40} else {10}
}

#################

n_yr_bioclimatic <- 1

temporal_range_env_data <- list(Ter = c(1980, 2019),
                                Mar = c(1982, 2021))

min_year_layers = list(Ter = 1982,
                       Mar = 1980)

limit_recs = 20000

myvars <- list(Ter = c(1980, 2019),
               Mar = c(1982, 2021))

varsdir <- list(Ter = "/media/seagate/boliveira/Land",
                Mar = "/media/seagate/boliveira/Marine")

n_tiles <- list(Ter = 40,
                Mar = 10)

