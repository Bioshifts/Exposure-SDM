
###############################################
slurm_job_singularity <- function(jobdir, logdir, sptogo, args,
                                  N_Nodes = 1, tasks_per_core = 1, cores = 1, time = "24:00:00", memory = "8G", partition = "normal",
                                  singularity_image,Rscript_file){
    
    # Start writing to this file
    sink(here::here(jobdir,paste0(sptogo,'.sh')))
    
    # the basic job submission script is a bash script
    cat("#!/bin/bash\n")
    
    cat("#SBATCH -N",N_Nodes,"\n")
    cat("#SBATCH -n",tasks_per_core,"\n")
    cat("#SBATCH -c",cores,"\n")
    cat("#SBATCH --partition",partition,"\n")
    cat("#SBATCH --mem=",memory,"\n", sep="")
    cat("#SBATCH --time=",time,"\n", sep="")
    
    cat("#SBATCH --job-name=",sptogo,"\n", sep="")
    cat("#SBATCH --output=",here::here(logdir,paste0(sptogo,".out")),"\n", sep="")
    cat("#SBATCH --error=",here::here(logdir,paste0(sptogo,".err")),"\n", sep="")
    # cat("#SBATCH --mail-type=ALL\n")
    # cat("#SBATCH --mail-user=brunno.oliveira@fondationbiodiversite.fr\n")
    
    cat(paste0("IMG_DIR='",singularity_image,"'\n"))
    
    cat("module purge\n")
    cat("module load singularity\n")
    
    cat("singularity exec --disable-cache $IMG_DIR Rscript",Rscript_file, args,"\n", sep=" ")
    
    # Close the sink!
    sink()
    
    # Submit to run on cluster
    system(paste("sbatch", here::here(jobdir, paste0(sptogo,'.sh'))))
    
}


###############################################
# create temporal pseudo-absences
# BA = raster of the background area to draw pseudo records
# env_range = range of environmental data to generate pseudo records vector of 
# sp_occ = species presence data. Should include the columns year, month, and cell numbers. Used to remove from the pseudo-absence data the same cells within the same dates as in the sp_occ.
create_temporal_pseudo_absences <- function(BA, env_range, sp_occ){
    # 1) get cells at the BA
    cells_BA <- terra::cells(BA, 1)[[1]]
    # 2) get dates
    env_range[1] <- env_range[1] + n_yr_bioclimatic # + n_yr_bioclimatic to be able to calculate bioclimatics for the previous n_yr_bioclimatic
    all_dates <- format(
        seq.Date(from = as.Date(paste0(env_range[1],"/01/01")), 
                 to = as.Date(paste0(env_range[2],"/12/01")), 
                 by = "month"), 
        "%m_%Y")
    # 3) sample random cells and dates
    random_cells <- sample(cells_BA, size = nrow(sp_occ), replace = TRUE)
    random_dates <- sample(all_dates, size = nrow(sp_occ), replace = TRUE)
    PA_cell_dates <- paste(random_cells,random_dates)
    # 4) check for the existence of any combination of cell dates from PA in sp_occ
    occ_cell_dates <- paste(sp_occ$cell, paste(sp_occ$month,sp_occ$year, sep="_"))
    run = 0
    while(any(PA_cell_dates %in% occ_cell_dates)){
        run <- run + 1
        # resample
        random_cells <- sample(cells_BA, size = nrow(sp_occ), replace = TRUE)
        random_dates <- sample(all_dates, size = nrow(sp_occ), replace = TRUE)
        PA_cell_dates <- paste(random_cells,random_dates)
        if(run > 100){
            all_possible_cells_dates <- expand.grid(cells_BA, all_dates)
            all_possible_cells_dates <- paste(all_possible_cells_dates[,1],all_possible_cells_dates[,2])
            all_possible_cells_dates <- all_possible_cells_dates[-which(all_possible_cells_dates %in% occ_cell_dates)]
            PA_cell_dates <- sample(all_possible_cells_dates, size = nrow(sp_occ), replace = TRUE)
            tmp <- strsplit(PA_cell_dates, " ")
            random_cells <- as.numeric(sapply(tmp, function(x) x[1]))
            random_dates <- sapply(tmp, function(x) x[2])
        }
    }
    # Create a PA dataset
    tmp <- lapply(random_dates, function(x) strsplit(x, "_")[[1]])
    PA_years <- sapply(tmp, function(x) x[2])
    PA_months <- sapply(tmp, function(x) x[1])
    
    PA_occ <- data.frame(year = PA_years, 
                         month = PA_months,
                         species = sp_occ$species[1],
                         cell = random_cells,
                         pa = 0)
    PA_xy <- terra::xyFromCell(BA, random_cells)    
    PA_occ <- cbind(PA_occ, PA_xy) 
    return(PA_occ)
}


###############################################
# make bioshifts v3 workable
# change colnames to reflect names in v1
bioshifts_fix_columns <- function(bioshifts){
    require(dplyr)
    bioshifts %>%
        dplyr::mutate(sp_name_std = gsub(" ","_",sp_name_checked),
                      Start = MIDPOINT_firstperiod,
                      End = MIDPOINT_secondperiod,
                      ShiftRate = case_when(
                          UnitRate == "m/year" ~ Rate/1000, 
                          UnitRate == "km/year" ~ Rate, 
                          TRUE ~ NA_real_),
                      UnitRate = "km/year") %>%
        dplyr::select(ID,
                      Type,
                      Param,
                      sp_name_std,
                      Start, 
                      End,
                      
                      Grain_size, 
                      N_periodes, 
                      Sampling, 
                      Category, 
                      Obs_type, 
                      Uncertainty_distribution,
                      Duration,
                      
                      ShiftRate,
                      UnitRate,
                      Areakm2,
                      LatExtentk,
                      Eco,
                      kingdom,
                      phylum,
                      class,
                      family)
}

###############################################
# subset bioshifts for use in sdms:
# 1) Latitude shifts
# 2) Terrestrial or marine shifts
# 3) Range shifts within the time-period that matches the environmental data
bioshifts_sdms_selection <- function(x){
    
    x %>%
        bioshifts_fix_columns %>%
        # filter latitude shifts
        filter(Type == "LAT") %>%
        # filter terrestrial or marine shifts
        filter(Eco == "Ter" | Eco == "Mar") %>%
        # Shifts marine or terrestrial within time period of the environmental data
        filter(Eco == "Ter" & (Start >= (temporal_range_env_data("Ter")[1] + n_yr_bioclimatic)) | 
                   (Eco == "Mar" & (Start >= (temporal_range_env_data("Mar")[1] + n_yr_bioclimatic)))) %>%
        # Remove marine birds
        filter(!(Eco == "Mar" & class == "Aves")) %>%
        # Remove freshwater fish
        filter(!(Eco == "Ter" & class == "Actinopterygi*")) %>%
        # Remove fungi, bacteria and mosses
        filter(!(kingdom == "Fungi" | kingdom == "Bacteria" | phylum == "Bryophyta"))
    
}

###############################################
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

###############################################
# Markdown table
nice_table <- function(x){
    DT::datatable(
        x,
        rownames = FALSE,
        extensions = 'Buttons',
        options = list(pageLength = 100, 
                       scrollY = "400px",
                       scrollX = T,
                       dom = 'Bft',
                       buttons = c('csv'),
                       lengthMenu = list(c(10,25,50,-1),
                                         c(10,25,50,"All"))))
}

###############################################
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



###############################################
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

###############################################
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

# calculate bioclimatics from a given date
# enter the date and calculate bioclimatics for the previous N years
# date_i = date to calculate the bioclimatics for the previous N years
# realm = realm to find the environmental data and function associated realm level bioclimatics
# sp_occ = to find cells and dates bioclimatics will be attached
# n_yr_bioclimatic = N previous years to calculate bioclimatics
bioclimatics_from_date <- function(all_layers,date_i,realm,sp_occ,n_yr_bioclimatic,n_cores=NULL){
    
    # subset occ data for date i
    sub_occ <- sp_occ %>% dplyr::filter(date == date_i)
    
    # get spatvect from data.frame
    occ_vect <- vect(data.frame(sub_occ[,c("x","y")]),geom=c("x","y"))
    
    # get range of dates to get the bioclimatics
    periods_i <- paste0("01_",date_i)
    periods_i <- as.Date(periods_i,"%d_%m_%Y")
    periods_i <- format(
        seq.Date(from = periods_i-(365*n_yr_bioclimatic), 
                 to = periods_i, 
                 by = "month"), 
        "%m_%Y")[1:12]
    
    
    # subset env layers for this period
    layers_i_pos <- unique(grep(paste(periods_i,collapse = "|"),all_layers))
    layers_i <- all_layers[layers_i_pos]
    layers_i <- rast(layers_i)
    
    # extract data
    terra::window(layers_i) <- terra::ext(terra::buffer(occ_vect,1))
    
    if(is.null(n_cores)){
        layers_i <- terra::extract(layers_i, occ_vect)
    } else {
        layers_i <- parallel::mclapply(1:nrow(occ_vect), 
                                       function(x) { terra::extract(layers_i, occ_vect[x]) }, 
                                       mc.cores = n_cores)
        layers_i <- data.frame(data.table::rbindlist(layers_i))
    }
    
    # calculate bioclimatics
    if(realm == "Ter"){
        bios_i <- bioclimatics_land_simple(layers_i)
    }
    if(realm == "Mar"){
        bios_i <- bioclimatics_ocean_simple(layers_i)
    }
    
    bios_i$ID <- paste(sub_occ$cell,sub_occ$date,sep = "_")
    
    return(bios_i)
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
    
    models_folder <- here::here(sdm_dir(realm), species, gsub("_",".",species),"models")
    
    existing_models <- list.dirs(models_folder, full.names = TRUE, recursive = FALSE)
    
    if(length(existing_models)==0){ # if there are no models, delete everything
        unlink(models_folder, recursive = TRUE)
        unlink(list.files(
            here::here(sdm_dir(realm),species,gsub("_",".",species)), 
            pattern = "models.out",
            full.names = TRUE) , recursive = TRUE)
    } else {
        existing_models_ids <- sapply(existing_models, function(x){
            tmp <- strsplit(x,"/")[[1]]
            tmp[length(tmp)]
        })
        goog_model_id <- sort(as.numeric(existing_models_ids), decreasing = TRUE)[1]
        
        models_folder_del <- existing_models[grep(goog_model_id,existing_models,invert = TRUE)]
        unlink(models_folder_del, recursive = TRUE)
        
        models_del <- list.files(
            here::here(sdm_dir(realm),species,gsub("_",".",species)), 
            pattern = "models.out",
            full.names = TRUE)
        models_del <- models_del[grep(goog_model_id,models_del,invert = TRUE)]
        unlink(models_del , recursive = TRUE)
        
    }
}

# fix biomod directory
apply_gsub_to_s4 <- function(obj, pattern, replacement) {
    # Get the slot names of the S4 object
    # obj = model_sp
    slot_names <- slotNames(obj)
    
    # Loop through each slot and apply gsub if the slot contains character data
    for (i in 1:length(slot_names)) {
        slot_value <- slot(obj, slot_names[i])
        if (is.character(slot_value)) {
            slot(obj, slot_names[i]) <- gsub(pattern, replacement, slot_value)
        } else {
            slot_value <- slotNames(slot(obj, slot_names[i]))
            if(length(slot_value)>1){
                for (j in 1:length(slot_value)) {
                    slot_value_j <- slot(slot(obj, slot_names[i]),slot_value[j])
                    if (is.character(slot_value_j)) {
                        slot(slot(obj, slot_names[i]),slot_value[j]) <- 
                            gsub(pattern, replacement, slot(slot(obj, slot_names[i]),slot_value[j]))
                    }
                }  
            }else{
                slot_value <- slot(obj, slot_names[i])
                if (is.character(slot_value)) {
                    slot(obj, slot_names[i]) <- 
                        gsub(pattern, replacement, slot(obj, slot_names[i]))
                }
            }
        }
    }
    
    return(obj)
}

get_CV_file <- function(realm,sptogo,type){
    if(type=="all"){
        CV_type <- "_CV_ens_all"
    }
    if(type=="RUN"){
        CV_type <- "_CV_ens_RUN"
    }
    # CV file address
    CV_file_add <- here::here(sdm_dir(realm),sptogo,paste0(sptogo,CV_type,".csv"))
    
    if(!file.exists(CV_file_add)){
        
        # Load full model
        output_dir <- here::here(sdm_dir(realm), sptogo)
        
        files_sdms <- list.files(
            here::here(output_dir,gsub("_",".",sptogo)), 
            full.names = TRUE,
            pattern = "models.out")  
        
        files_sdms <- files_sdms[-grep("ensemble",files_sdms)]
        if(files_sdms > 1){
            delete_duplicated_models(realm = realm, species = sptogo)
            
            files_sdms <- list.files(
                here::here(output_dir,gsub("_",".",sptogo)), 
                full.names = TRUE,
                pattern = "models.out")  
            files_sdms <- files_sdms[-grep("ensemble",files_sdms)]
        }
        model_sp <- get(load(files_sdms))
        
        model_sp <- apply_gsub_to_s4(model_sp, pattern = "/lustre/oliveirab", replacement = "/scratch/boliveira")
        
        models_folder <- here::here(sdm_dir(realm), sptogo, gsub("_",".",sptogo),"models")
        existing_model <- list.dirs(models_folder, full.names = TRUE, recursive = FALSE)
        existing_model <- strsplit(existing_model,"/")[[1]]
        existing_model <- existing_model[length(existing_model)]
        current_model <- model_sp@modeling.id
        if(!current_model==existing_model){
            model_sp <- apply_gsub_to_s4(model_sp, pattern = current_model, replacement = existing_model)
        }
        
        if(type=="all"){
            # calculate ensemble
            ens_model_sp <- BIOMOD_EnsembleModeling(
                model_sp, 
                chosen.models = 'all',
                em.by = 'all',
                metric.select = "TSS",
                metric.eval = "TSS")
        }
        if(type=="RUN"){
            # calculate ensemble
            ens_model_sp <- BIOMOD_EnsembleModeling(
                model_sp, 
                metric.select = "TSS",
                metric.eval = "TSS")
        }
        
        
        # get evaluations from ensemble model
        CV_file = get_evaluations(ens_model_sp)
        
        write.csv(CV_file, 
                  CV_file_add,
                  row.names = FALSE)
        rm(ens_model_sp)
    } else {
        CV_file <- read.csv(CV_file_add)
    }
    
}

## Select the best partition for job submission
ter_partitions <- c("normal-amd","normal","bigmem-amd","bigmem")
limits <- data.frame(partition = ter_partitions,
                     max_mem_node = c(250, 125, 1000, 500),
                     max_mem_user = c(2000, 2000, 2000, 2000),
                     max_cpu = c(512, 448, 128, 112))

select_partition <- function(request_mem, request_cpu, limits){
    test <- select_partition_inset(request_mem = request_mem, request_cpu = request_cpu, limits = limits)
    if(nrow(test)==0){
        cat("\n\rThere are no enough resources for running submitting job.\nWaiting for resources to become available...\n")
    }
    while(nrow(test)==0){
        
        
        Sys.sleep(10)
        test <- select_partition_inset(request_mem, request_cpu, limits)
        
    } 
    return(as.character(test$partition[1]))
}
select_partition_inset <- function(request_mem, request_cpu, limits){
    limits <- limits %>% filter((max_mem_node > request_mem))
    
    test_partition <- system("squeue --format='%.50P' --me", intern = TRUE)[-1]
    if(length(test_partition)==0){
        return(limits[1,])
    }else{
        test_partition <- gsub(" ","",test_partition)
        using <- data.frame(table(test_partition))
        names(using)[1] <- "partition"
        using <- merge(using, limits, all = TRUE)
        using[is.na(using)] <- 0
        using$using_mem <- request_mem * using$Freq
        using$using_cpu <- request_cpu * using$Freq
        
        tmp <- using %>% filter((max_mem_user > using$using_mem) & (max_cpu > using$using_cpu))
        return(tmp)
    }
}

# check if species has sdms for all possible years
check_if_has_sdms_for_all_shifts <- function(all_sps,bioshifts){
    
    require("pbapply")
    
    my_sdms <- pbapply::pblapply(1:nrow(all_sps), function(i){
        
        # sp_i <- "Centropristis_striata"
        # sp_i_realm <- "Mar"
        sp_i <- gsub(" ","_",all_sps$sps[i])
        sp_i_realm <- all_sps$realm[i]
        
        # test if have sdm for sp_i
        test <- list.files(
            here::here(sdm_dir(sp_i_realm),sp_i,gsub("_",".",sp_i)), 
            pattern = "ensemble.models.out")
        
        test <- length(test) > 0
        
        if(test){ # if yes, check if have projections for all possible shifts
            
            # all possible shifts for sp_i?
            shift_info <- filter(bioshifts, sp_name_std == sp_i)
            shift_info <- select(shift_info, c(ID, Start, End))
            shift_info <- unique(shift_info)
            
            ID_i <- unique(shift_info$ID)
            # for each ID_i, look if has projections for all years
            tmp <- data.frame(sps=sp_i,
                              realm=sp_i_realm,
                              ID=ID_i)
            
            tmp$I_have_sdms <- sapply(ID_i, function(x){
                shift_info_i <- shift_info[which(shift_info$ID==x),]
                years_ID_i <- round(shift_info_i$Start,0):round(shift_info_i$End,0)
                # check
                # Focus on SA for now
                sdms_sp_i <- list.files(here(sdm_dir(sp_i_realm),sp_i,gsub("_",".",sp_i)),pattern = "SA")
                # get ensemble models
                sdms_sp_i_ens <- sdms_sp_i[grep(" ens",sdms_sp_i)]
                # get ID_i models
                sdms_sp_i_ens <- sdms_sp_i_ens[grep(x,sdms_sp_i_ens)]
                # check if all years exist
                sdms_sp_i_ens <- all(sapply(years_ID_i, function(x){any(grepl(x,sdms_sp_i_ens))}))
                sdms_sp_i_ens
            })
            
            
        } else {
            tmp <- data.frame(sps = sp_i, realm = sp_i_realm, ID = NA, I_have_sdms = FALSE)
        }
        return(tmp)
    })
    my_sdms <- do.call(rbind, my_sdms)
    return(my_sdms)
}