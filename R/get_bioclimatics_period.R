get_bioclimatics_period <- function(p_1,p_2,n_yr_bioclimatic,BA,n_cores,env_data){
    ### Get bioclimatics for the BA at time periods
    # time periods for model i
    p_start_i <- (p_1-n_yr_bioclimatic):p_1
    p_end_i <- (p_2-n_yr_bioclimatic):p_2
    
    ### start
    start_period <- lapply(p_start_i, function(x) paste(x, sprintf("%02d", 1:12), sep = "_"))
    start_period <- do.call(c,start_period)
    start_period <- paste(start_period,collapse = "|")
    # load env data for model i
    env_data_start <- env_data[grep(start_period,env_data)]
    env_data_start_names <- strsplit(env_data_start,"[/]")
    env_data_start_names <- sapply(env_data_start_names, function(x){
        tmp <- x[length(x)]
        gsub(".tif","",tmp)
    })
    env_data_start <- terra::rast(env_data_start)
    names(env_data_start) <- env_data_start_names
    ### end
    end_period <- lapply(p_end_i, function(x) paste(x, sprintf("%02d", 1:12), sep = "_"))
    end_period <- do.call(c,end_period)
    end_period <- paste(end_period,collapse = "|")
    # load env data for model i
    env_data_end <- env_data[grep(end_period,env_data)]
    env_data_end_names <- strsplit(env_data_end,"[/]")
    env_data_end_names <- sapply(env_data_end_names, function(x){
        tmp <- x[length(x)]
        gsub(".tif","",tmp)
    })
    env_data_end <- terra::rast(env_data_end)
    names(env_data_end) <- env_data_end_names
    
    # Mask data to the BG area
    cells_BA <- cells(BA,1)
    cells_BA <- cells_BA$BA
    BG_data <- mclapply(list(env_data_start,env_data_end), 
                        function(x) x[cells_BA],
                        mc.cores = 2)
    BG_data_start <- BG_data[[1]]
    BG_data_end <- BG_data[[2]]
    BG_xys = xyFromCell(env_data_start[[1]],cells_BA)
    rm(BG_data)
    # set chunks for faster extraction
    x=1:length(BG_data_start[,1])
    chunks <- split(x, cut_number(x, n_cores))
    # calculate bioclimatics
    ### start
    if(realm == "Mar"){
        bio_back_start <- mclapply(chunks, function(x) {
            data.frame(bioclimatics_ocean(BG_data_start[x,], 
                                          myvars = myvars, 
                                          n_yr_bioclimatic = n_yr_bioclimatic))
        }, mc.cores = n_cores)
    } else {
        bio_back_start <- mclapply(chunks, function(x) {
            data.frame(bioclimatics_land(BG_data_start[x,], 
                                         myvars = myvars, 
                                         n_yr_bioclimatic = n_yr_bioclimatic))
        }, mc.cores = n_cores)
    }
    bio_back_start <- rbindlist(bio_back_start)
    bio_back_start <- data.frame(bio_back_start)
    
    ### end
    if(realm == "Mar"){
        bio_back_end <- mclapply(chunks, function(x) {
            data.frame(bioclimatics_ocean(BG_data_end[x,], 
                                          myvars = myvars, 
                                          n_yr_bioclimatic = n_yr_bioclimatic))
        }, mc.cores = n_cores)
    } else {
        bio_back_end <- mclapply(chunks, function(x) {
            data.frame(bioclimatics_land(BG_data_end[x,], 
                                         myvars = myvars, 
                                         n_yr_bioclimatic = n_yr_bioclimatic))
        }, mc.cores = n_cores)
    }
    bio_back_end <- rbindlist(bio_back_end)
    bio_back_end <- data.frame(bio_back_end)
    
    return(list(back_start = bio_back_start, back_end = bio_back_end, xy = BG_xys))
    
}
