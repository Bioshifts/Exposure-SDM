# from a distance object to data.frame
dist.to.df <- function(d){
    size <- attr(d, "Size")
    return(
        data.frame(
            subset(expand.grid(row=2:size, col=1:(size-1)), row > col),
            distance=as.numeric(d),
            row.names = NULL
        )
    )
}

# Markdown table
nice_table <- function(x){
    DT::datatable(
        x,
        rownames = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = 100, 
                       scrollY = "400px",
                       scrollX = T,
                       dom = 'Blfrtip',
                       buttons = c('csv', 'excel'),
                       lengthMenu = list(c(10,25,50,-1),
                                         c(10,25,50,"All"))))
}

# submit system job
# args = spgoing
# code_dir= here::here("R/get_data_sps_sdms.R")
# Rout_file = here::here(jobs.dir, paste0(i,"_",gsub(" ","_",spgoing), ".Rout"))
systemjob <- function(args,code_dir,Rout_file){
    system(paste("R CMD BATCH --vanilla", 
                 dQuote(paste("--args", args), q = "C"), 
                 code_dir,
                 Rout_file))
}

slurmjobsingularity <- function(job_args, slurm_args, Rscript_file){
    
    slurm_cmd <- paste("srun",
                       "-N", slurm_args$n_nodes, 
                       "-n", slurm_args$n_task_per_node, 
                       "-c", slurm_args$n_cores, 
                       paste0("--time=", slurm_args$time), 
                       paste0("--mem=", slurm_args$mem), 
                       "-J", slurm_args$job_name,
                       "singularity exec --disable-cache $HOME/work_t_cesab/brunno/brunnospatial.sif",
                       "Rscript", Rscript_file, job_args)
    
    system(slurm_cmd)
}

# clean var names from raster files
cleanNamesFiles <- function(orinames, tile = F){
    # rename layers for easy reading
    if(any(grepl("ph|o2",orinames))){
        pos <- grep("ph|o2",orinames)
        names.tmp <- sapply(orinames[pos], function(x){
            tmp = strsplit(x, "_")[[1]]
            if(tile){
                return(paste(tmp[-4:-5],collapse = "_"))
            }else{
                return(paste(tmp[-4],collapse = "_"))
            }
        })
        orinames[pos] <- names.tmp
    }
    if(!any(grepl("ph|o2",orinames))){
        pos <- which(!grepl("ph|o2",orinames))
        names.tmp <- sapply(orinames[pos], function(x){
            tmp = strsplit(x, "_")[[1]]
            if(tile){
                return(paste(tmp[-4],collapse = "_"))
            }else{
                return(paste(tmp,collapse = "_"))
            }
        })
        orinames[pos] <- names.tmp
    }
    return(orinames)
}

cleanNamesLayers <- function(orinames, tile = F){
    # rename layers for easy reading
    if(any(grepl("ph|o2",orinames))){
        pos <- grep("ph|o2",orinames)
        names.tmp <- sapply(orinames[pos], function(x){
            tmp = strsplit(x, "_")[[1]]
            tmp = tmp[-4]
            return(paste(tmp,collapse = "_"))
        })
        orinames[pos] <- names.tmp
    }
    if(any(grepl("SST",orinames))){
        pos <- grep("SST",orinames)
        names.tmp <- sapply(orinames[pos], function(x){
            tmp = strsplit(x, "_")[[1]]
            tmp1 = tmp[1] # varname
            tmp2 = substr(tmp[2],1,4) # year
            tmp3 = substr(tmp[2],5,6) # month
            if(tile){
                tmp4 = tmp[3] # tile
                return(paste(c(tmp1,tmp2,tmp3,tmp4),collapse = "_"))
            }else{
                return(paste(c(tmp1,tmp2,tmp3),collapse = "_"))
            }
        })
        orinames[pos] <- names.tmp
    }
    if(any(grepl("pr|tas",orinames))){
        pos <- grep("pr|tas",orinames)
        names.tmp <- sapply(orinames[pos], function(x){
            tmp = strsplit(x, "_")[[1]]
            if(tile){
                return(paste(c(tmp[2],tmp[4],tmp[3],tmp[6]),collapse = "_"))
            }else{
                return(paste(c(tmp[2],tmp[4],tmp[3]),collapse = "_"))
            }
        })
        orinames[pos] <- names.tmp
    }
    return(orinames)
}

##################################
# Get time range
# enter dates (data frame with 2 columns - year and month) and N years needed to get a time range of all months/years before N years
GetTimeRange <- function(dates, n_years){
    dates <- dates[, c("year","month")]
    uniquedates <- paste(dates,collapse = '_')
    dates2go.tmp <- 
        lapply(1:nrow(dates), function(i){
            mydate <- dates[i, c("year","month")]
            mydate <- paste(mydate[,1] - n_years, 
                            mydate[,2], 
                            "01", 
                            sep = "/")
            mydate <- as.Date(mydate)
            daterange <- seq(mydate, by = "month", length.out = n_years*12)
            format(daterange, "%Y_%m")
        })
}

##################################
# calculate bioclimatic variables from time series
bioclimatics_ocean_simple <- function(env_data){
    x = env_data
    # annual average (equivalent to mean annual temperature)
    amean <- apply(x, 1, mean, na.rm = T) 
    # seasonality (equivalent to temperature seasonality)
    asd <- apply(x, 1, sd, na.rm = T) 
    # minimum (equivalent to max Temperature of Coldest Quarter
    amin <- apply(x, 1, min, na.rm = T)
    # maximum (equivalent to min Temperature of Warmest Quarter 
    amax <- apply(x, 1, max, na.rm = T)
    
    return(data.frame(mean = amean, sd = asd, min = amin, max = amax))
}

# calculate bioclimatic variables from time series
bioclimatics_ocean <- function(env_data){
    tmp <- bioclimatics_ocean_inside(env_data)
    return(tmp)
}
bioclimatics_ocean_inside <- function(x){
    
    window <- function(x)  { 
        lng <- length(x)
        x <- c(x,  x[1:3])
        m <- matrix(ncol=3, nrow=lng)
        for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
        apply(m, MARGIN=1, FUN=sum)
    }
    
    # annual average (equivalent to mean annual temperature)
    amean <- apply(x, 1, mean, na.rm = T) 
    
    tmp <- t(apply(x, 1, window)) / 3
    
    # seasonality (equivalent to temperature seasonality)
    asd <- apply(x, 1, sd, na.rm = T) 
    # minimum (equivalent to Mean Temperature of Coldest Quarter
    amin <- apply(tmp, 1, min, na.rm = T)
    # maximum (equivalent to Mean Temperature of Warmest Quarter 
    amax <- apply(tmp, 1, max, na.rm = T)
    
    return(data.frame(mean = amean, sd = asd, min = amin, max = amax))
}

##################################
# calculate bioclimatic (terrestrial like >> tas, tasmim, tasmax, pr) variables from timeseries
bioclimatics_land_simple <- function(env_data){
    tmin <- env_data[,grep("tmin",names(env_data))]
    tmax <- env_data[,grep("tmax",names(env_data))]
    pr <- env_data[,grep("prec",names(env_data))]
    tas <- (tmin + tmax) / 2
    
    # mean annual temperature
    amean <- apply(tas, 1, mean, na.rm = T) 
    # temperature seasonality
    asd <- apply(tas, 1, sd, na.rm = T) 
    # Max Temperature of Warmest Quarter 
    amax <- apply(tmax, 1, max, na.rm = T)
    # Min Temperature of Coldest Quarter
    amin <- apply(tmin, 1, min, na.rm = T)
    
    # mean annual precipitation
    pmean <- apply(pr, 1, mean, na.rm = T) 
    # precipitation seasonality
    psd <- apply(pr, 1, sd, na.rm = T) 
    # Precipitation of Wettest Quarter 
    pmin <- apply(pr, 1, max, na.rm = T)
    # Precipitation of Driest Quarter 
    pmax <- apply(pr, 1, min, na.rm = T)
    
    # return bioclimatics
    return(data.frame(mat = amean, seat = asd, mint = amin, maxt = amax,
                      map = pmean, seap = psd, minp = pmin, maxp = pmax))
}

# calculate bioclimatic (terrestrial like >> tas, tasmim, tasmax, pr) variables from timeseries
bioclimatics_land <- function(env_data){
    resp <- bioclimatics_land_inside(env_data)
    return(resp)
}
bioclimatics_land_inside <- function(env_data){
    
    window <- function(x)  { 
        lng <- length(x)
        x <- c(x,  x[1:3])
        m <- matrix(ncol=3, nrow=lng)
        for (i in 1:3) { m[,i] <- x[i:(lng+i-1)] }
        apply(m, MARGIN=1, FUN=sum)
    }
    
    tmin <- env_data[,grep("tasmin",names(env_data))]
    tmax <- env_data[,grep("tasmax",names(env_data))]
    pr <- env_data[,grep("pr",names(env_data))]
    tas <- (tmin + tmax) / 2
    
    tmp <- t(apply(tas, 1, window)) / 3
    
    # mean annual temperature
    amean <- apply(tas, 1, mean, na.rm = T) 
    # temperature seasonality
    asd <- apply(tas, 1, sd, na.rm = T) 
    # Mean Temperature of Warmest Quarter 
    amax <- apply(tmp, 1, max, na.rm = T)
    # Mean Temperature of Coldest Quarter
    amin <- apply(tmp, 1, min, na.rm = T)
    
    # precip by quarter (3 months)		
    wet <- t(apply(pr, 1, window))
    
    # mean annual precipitation
    pmean <- apply(pr, 1, mean, na.rm = T) 
    # precipitation seasonality
    psd <- apply(pr, 1, sd, na.rm = T) 
    # Precipitation of Wettest Quarter 
    pmin <- apply(wet, 1, max, na.rm = T)
    # Precipitation of Driest Quarter 
    pmax <- apply(wet, 1, min, na.rm = T)
    
    # return bioclimatics
    return(data.frame(mat = amean, seat = asd, mint = amin, maxt = amax,
                      map = pmean, seap = psd, minp = pmin, maxp = pmax))
}

##################################
# get data from each cells in time-series data

# env_var <- env_vars[[1]]
env_time_cell <- function(dates, env_data, temporalrange, mc = 10){
    # get time range for each cell
    dates$month <- sprintf("%02d", dates$month)
    dates$id <- sapply(1:nrow(dates), function(i) {paste(dates$year[i], dates$month[i], "01", sep = "/")})
    dates2go <- dates[-which(duplicated(dates$id)),]
    timerange <- mclapply(1:nrow(dates2go), function(i){
        # get temporal range
        mydate <- dates2go[i, c("year","month")]
        mydate <- paste(mydate[,1] - temporalrange, 
                        mydate[,2], 
                        "01", 
                        sep = "/")
        mydate <- as.Date(mydate)
        daterange <- seq(mydate, by = "month", length.out = temporalrange*12)
        daterange <- sapply(daterange, function(x) {
            x <- strsplit(as.character(x), "-")[[1]]
            paste(x[1],x[2],sep = "_")
        })
        return(daterange)
    }, mc.cores = mc)
    # get dates from env var
    dates_from_env <- names(env_data)
    dates_from_env <- sapply(dates_from_env, function(x){
        tmp <- strsplit(x,"_")[[1]][-1] 
        paste(tmp[1],tmp[2],sep = "_")
    })
    # get positions of each time period for variables
    pos <- lapply(timerange, function(x){
        which(dates_from_env %in% x)
    })
    # get climate from time periods at each cell
    climtime <- mclapply(1:length(pos), function(y){
        as.numeric(env_data[y, pos[[y]]])
    },mc.cores = mc)
    # feed
    names(climtime) <- dates2go$id
    myList <- lapply(dates$id, function(x) climtime[[x]])
    return(myList)
}

##################################
# calculates bioclimatic variables from raw climate data
# allows specifying time period for calculate of bioclimatic variables
# temporal range in years. Represents the period from which raw variables will be averaged for the calculus of bioclimatics
# occ is a data.frame object with at least the colnames: cell, year, month
# myvars is the list of raw variables that will be used for calculating bioclimatics
# varsdir is the directory when variables are stored
# environment should be "T" (terrestrial), "M" (marine), "A" (aquatic). The way bioclimatic variables are calculated depend on the environment

# occ = spoccur[1:5,]
# myvars = c("tasmax", "tasmin", "tas", "pr")
# varsdir = "/media/seagate/boliveira/Land"
# temporalrange = 2
# env_type = "T"

# occ = spoccur[1:5,]
# myvars = c("ph", "o2", "SST")
# varsdir = "/media/seagate/boliveira/Marine"
# temporalrange = 2
# env_type = "M"

minorThreat <- function(occ, myvars, varsdir, crs = "+proj=longlat +datum=WGS84 +no_defs", temporalrange = 2, env_type = NULL, limit = NULL){
    if(is.null(env_type)){
        stop("env_type type should be specified.
             The way bioclimatic variables are calculated depends on the environment type.")
    }
    if(!env_type %in% c("T", "M", "A")){
        stop("env_type should be either: 'T' (terrestrial), 'M' (marine) or 'A' (aquatic).
             The way bioclimatic variables are calculated depends on the environment type.")
    }
    # load climate variables
    files <- list.files(here::here(varsdir), full.names = T, recursive = T, pattern = "tif$")
    files <- files[-grep("raw|model",files)]
    # Select layers within the time-frame of interest
    timeframe <- paste(c((min(occ$year)-2): max(occ$year)),collapse = "|")
    pos <- grep(timeframe, files)
    files <- files[pos]
    # get data for each env variable
    cellstogo <- occ$cell
    coords <- st_as_sf(SpatialPoints(occ[,c("decimalLongitude", "decimalLatitude")]))
    poly <- st_buffer(coords, dist = 1)
    st_crs(poly) <- crs
    env_cells <- list()
    for(i in 1:length(myvars)){
        files.tmp <- files[grep(myvars[i],files)]
        layers.tmp <- read_stars(files.tmp, proxy = TRUE, crs = crs)
        names(layers.tmp) <- cleanNames(layers.tmp)
        # crop star object for faster extraction
        layers.tmp <- layers.tmp[poly]
        # get env data from each cell
        env_cells.tmp <- st_extract(layers.tmp, st_coordinates(coords))
        env_cells.tmp$cell <- cellstogo
        env_cells[[i]] <- env_cells.tmp
    }
    names(env_cells) <- myvars
    # get env data from cells at dates
    dates = occ[,c("year","month")]
    envtimecells <- lapply(1:length(env_cells), function(i) {
        env_time_cell(dates = dates,
                      env_cells = env_cells[[i]],
                      temporalrange = temporalrange)
    })
    names(envtimecells) <- myvars
    # calc bioclimatics for each cell
    if(env_type == "M"){
        bios <- lapply(envtimecells, function(x) {
            tmp <- lapply(x, function(i){
                bioclimatics(i)
            })
            return(do.call(rbind,tmp))
        })
        # fix bioclimatic names
        for(i in 1:length(bios)) { names(bios[[i]]) <- paste(myvars[i],names(bios[[i]]),sep = "_")}
        names(bios) <- NULL
        bios <- do.call(cbind,bios)
    } else {
        bios <- bioclimatics_land(envtimecells)
    }
    bios <- cbind(cell=cellstogo,bios)
    if(is.null(limit)){
        return(bios)
    } else {
        if(nrow(bios) >= limit){
            bios <- bios[sample(1:nrow(bios),limit),]
        }
        return(bios)
    }
}

# Clean coordinates
drugfree <- function(x, inverse = FALSE, my.mask){
    original_cols <- names(x)
    # Clean
    ## get cells
    cells. <-  terra::extract(my.mask,x[,c("decimalLongitude", "decimalLatitude")],cells=TRUE)
    x$cell <- cells.$cell
    x$layer <- cells.[,2]
    # remove/keep cells falling in the ocean/Land
    rm <- x$layer
    if(inverse){
        x <- x[which(rm==0),] # remove land / keep ocean
    } else{
        x <- x[which(rm==1),] # keep land / remove ocean
    }
    x = x[,-"layer"]
    ## Remove duplicates: One coordinate per species/cellID/month/year
    ## get combination species/cellID/month/year
    rm <- duplicated(x[,c("species", "cell", "year", "month")])
    if(any(rm)){
        x <- x[-which(rm),]
    }
    ## Remove other potential issues
    x = CoordinateCleaner::clean_coordinates(x = x,
                                             lon = "decimalLongitude",
                                             lat = "decimalLatitude",
                                             tests = c("capitals", "centroids", "equal",
                                                       "gbif", "institutions", "zeros"),
                                             value = "clean")
    return(x)
}

##################################
# Unzip large file
decompress_file <- function(directory, file, .file_cache = FALSE) {
    
    if (.file_cache == TRUE) {
        print("decompression skipped")
    } else {
        
        # Set working directory for decompression
        # simplifies unzip directory location behavior
        wd <- getwd()
        setwd(directory)
        
        # Run decompression
        decompression <-
            system2("unzip",
                    args = c("-o", # include override flag
                             file),
                    stdout = TRUE)
        
        # uncomment to delete archive once decompressed
        # file.remove(file) 
        
        # Reset working directory
        setwd(wd); rm(wd)
        
        # Test for success criteria
        # change the search depending on 
        # your implementation
        if (grepl("Warning message", tail(decompression, 1))) {
            print(decompression)
        }
    }
}

##################################
# Get Copernicus Global ocean biogeochemistry hindcast
GetCopernicusCMEMS <- function(user, password, variables, nc_path, years, months, depths = 1, keep.original = FALSE, over.write = TRUE){
    
    if(.Platform$OS.type == "windows"){
        stop("This function only works on Unix system")
    }
    
    new.dir <- here::here(nc_path)
    if(!dir.exists(new.dir)){
        dir.create(new.dir,recursive = T)
    }
    
    for(y in 1:length(years)){
        for(m in 1:length(months)){
            
            filename <- paste0("mercatorfreebiorys2v4_global_mean_",years[y],months[m],".nc")
            path.file <- file.path(nc_path, filename)
            if(!over.write){
                over.write <- !file.exists(path.file)
            }
            if(over.write){
                system(paste0("wget -q --user=", user, " --password=", password,
                              " -O ", path.file,
                              " ftp://my.cmems-du.eu/Core/GLOBAL_MULTIYEAR_BGC_001_029/cmems_mod_glo_bgc_my_0.25_P1M-m/",
                              years[y], "/", filename))
                
                for(v in 1:length(variables)){
                    for(d in 1:length(depths)){
                        myfile.tmp <- brick(path.file, varname=variables[v])[[depths[d]]]
                        writeRaster(myfile.tmp, 
                                    file.path(nc_path, paste0(paste(variables[v],years[y],months[m],depths[d],sep = "_"), ".tif")),
                                    format="GTiff",
                                    overwrite=TRUE)
                    }
                }
                
                if(!keep.original){
                    unlink(path.file)
                }
            }
        }
    }
}

##################################
# Get Copernicus Global SST
GetCopernicusSST <- function(user, password, nc_path, years, months, model.raster = NULL){
    
    if(.Platform$OS.type == "windows"){
        stop("This function only works on Unix system")
    }
    
    new.dir <- here::here(nc_path)
    if(!dir.exists(new.dir)){
        dir.create(new.dir,recursive = T)
    }
    
    tmp.dir <- here::here(new.dir,"tmp")
    if(!dir.exists(tmp.dir)){
        dir.create(tmp.dir,recursive = T)
    }
    
    for(y in 1:length(years)){
        for(m in 10:length(months)){
            # list files
            files <- system(paste0("wget --user=", user, " --password=", password,
                                   " --no-remove-listing --spider",
                                   " ftp://my.cmems-du.eu/Core/SST_GLO_SST_L4_REP_OBSERVATIONS_010_011/METOFFICE-GLO-SST-L4-REP-OBS-SST/",
                                   years[y], "/", months[m], "/"))
            files <- read.table(".listing")
            files <- files$V9[-1:-2]
            
            path.files <- sapply(files, function(x) here::here(tmp.dir, x))
            
            unlink(".listing")
            
            for(f in 1:length(files)){
                
                cat("\r Downloading year", years[y], "month", months[m], "day", f)
                
                system(paste0("wget -q --user=", user, " --password=", password,
                              " -O ", path.files[f],
                              " ftp://my.cmems-du.eu/Core/SST_GLO_SST_L4_REP_OBSERVATIONS_010_011/METOFFICE-GLO-SST-L4-REP-OBS-SST/",
                              years[y], "/", months[m], "/", files[f]))
                
            }
            
            cat("\nResampling\n")
            
            cl <- makeCluster(detectCores()-10)
            clusterExport(cl, c("path.files", "model.raster"))
            
            pbsapply(path.files, FUN = function(x){
                myfile.tmp <- raster::raster(x)
                myfile.tmp <- raster::resample(myfile.tmp, model.raster)
                
                raster::writeRaster(myfile.tmp, 
                                    paste0(x, ".tif"),
                                    format="GTiff",
                                    overwrite=TRUE)
                
                unlink(x)
            },
            cl = cl)
            
            stopCluster(cl)
            
            cat("\nStacking and saving")
            
            myraster.files <- list.files(tmp.dir, full.names = TRUE)
            mystack <- stack(myraster.files)
            myraster <- raster::calc(mystack, mean)
            
            writeRaster(myraster, 
                        paste0(here::here(new.dir, paste0(paste("SST",years[y],months[m],sep="_"),".tif"))),
                        format="GTiff",
                        overwrite=TRUE)
            
            unlink(myraster.files)
        }
    }
    unlink(tmp.dir, recursive = T)
}

##################################
# First letter upper case
.simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep = "", collapse = " ")
}

# Fix species names
fix.spnames <- function(x){
    # remove letters with accent
    tmp <- iconv(x, "latin1", "ASCII", "")
    # remove var and subspecies
    tmp <- gsub("_x_","_", tmp, ignore.case = T)
    tmp <- gsub("_ssp_","_", tmp, ignore.case = T)
    tmp <- gsub("_ssp.","_", tmp, ignore.case = T)
    tmp <- gsub("_subsp_","_", tmp, ignore.case = T)
    tmp <- gsub("_subsp.","_", tmp, ignore.case = T)
    tmp <- gsub("_var_","_", tmp, ignore.case = T)
    tmp <- gsub("_var[.]","_", tmp, ignore.case = T)
    tmp <- gsub("__","_", tmp, ignore.case = T)
    
    tmp = strsplit(tmp, "_")[[1]]
    
    if(length(tmp) == 1){
        tmp = 0
    } else {
        # remove special characters
        tmp = sapply(tmp, function(i) str_replace_all(i, "[^[:alnum:]]", ""))
        if(any(
            sapply(tmp, function(i) any(is.na(i) | grepl("NA", i) | i == "sp" | i == "x" | i == "X")))){
            tmp = 0
        } else {
            tmp <- sapply(tmp, tolower)
            tmp[1] <- .simpleCap(tmp[1])
            tmp = paste(tmp,collapse = "_")
        }
    }
    return(tmp)
}

# delete duplicated SDM generated from biomod
delete_duplicated_models <- function(realm,species){
    # models
    files_sdms <- list.files(
        here::here(sdm_dir(realm),species,gsub("_",".",species),"models"), 
        full.names = TRUE)  
    
    if(length(files_sdms) > 1){ # if there are > 1 model, keep the most recent one
        del <- order(file.info(files_sdms)$ctime,decreasing = TRUE)[-1]
        
        unlink(files_sdms[del],recursive = TRUE)
    }
    
    # link to models
    files_sdms_model <- list.files(
        here::here(sdm_dir(realm),species,gsub("_",".",species)), 
        pattern = "models.out",
        full.names = TRUE)  
    
    if(length(files_sdms)==0){ # if there are mo models, delete all
        unlink(files_sdms_model)
    }
    
    if(length(files_sdms_model) > 2){ # if there are > 2 (simple + ensemble)
        
        which_sdm <- strsplit(files_sdms,split = "/")[[1]]
        which_sdm <- which_sdm[length(which_sdm)]
        
        del <- files_sdms_model[!grepl(which_sdm, files_sdms_model)]
        
        unlink(del)
        
    }
}