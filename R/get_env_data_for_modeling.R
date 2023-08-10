get_env_data_for_modeling <- function(realm, myvars, varsdir){
    env_data <- list()
    for(j in 1:length(myvars)){
        tmp <- list.files(here::here(varsdir,myvars[j]),pattern = ".tif", full.names = T)
        names(tmp) <- cleanNamesFiles(tmp)
        env_data[[j]] <- tmp
    }
    env_data <- unlist(env_data)
    return(env_d