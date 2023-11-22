# load ecoregions
get_ecoregions <- function(realm, sptogo, PresAbs, mask.ras, return.shp = TRUE, return.raster = TRUE, varsdir, check_if_exists = TRUE, output_dir){
    
    if(all(!return.shp,!return.raster)){
        stop("Use return.shp = TRUE or return.raster = TRUE.\nAt least one type of file should be returned.")
    }
    
    shptosave <- here::here(output_dir,paste0(sptogo,"_BA.shp"))
    rastertosave <- here::here(output_dir,paste0(sptogo,"_BA.tiff"))
    
    if(return.shp){
        
        if(check_if_exists){
            test <- file.exists(shptosave)
            if(test){ # if file exists load the BA
                ecoreg_shp <- terra::vect(shptosave)
            } else { # if not, create the BA
                ecoreg_shp <- get_BA_shp(varsdir, realm, PresAbs)
                terra::writeVector(ecoreg_shp,shptosave,overwrite=TRUE)
            }
        } else {
            ecoreg_shp <- get_BA_shp(varsdir, realm, PresAbs)
            terra::writeVector(ecoreg_shp,shptosave,overwrite=TRUE)
        }
        
        
    } else {
        ecoreg_shp = NULL
    }
    
    
    if(return.raster){
        
        if(check_if_exists){
            test <- file.exists(rastertosave)
            if(test){ # if file exists load the BA
                ecoreg_raster <- terra::rast(rastertosave)
            } else { # if not, create the BA
                ecoreg_raster <- get_BA_raster(varsdir, realm, PresAbs)
                terra::writeRaster(ecoreg_raster,rastertosave)
            }
        } else {
            ecoreg_raster <- get_BA_raster(varsdir, realm, PresAbs)
            terra::writeRaster(ecoreg_raster,rastertosave)
        }
        
        
    } else {
        ecoreg_raster = NULL
    }
    
    return(list(raster_file = ecoreg_raster, shape_file = ecoreg_shp))
    
}
# Inside functions
get_wwf_shp <- function(varsdir){
    # check if shape was downloaded and if not download it
    if(dir.exists(here::here(varsdir,"WWF_ecoregions"))){
        ecoreg_shp <- terra::vect(here::here(varsdir,"WWF_ecoregions/official/wwf_ecorealms.shp"))
        return(ecoreg_shp)
    } else {
        cat("Downloading WWF shape file")
        ecoreg_shp <- speciesgeocodeR::WWFload(here::here(varsdir))
        ecoreg_shp <- terra::vect(ecoreg_shp)
        # Create the EcoRealm layer
        # but first fix an issue >> there is conflict with the code NA for Neartic. Fix to avoid confusion with real.NA.
        tmp <- as.character(ecoreg_shp$REALM)
        tmp[which(ecoreg_shp$REALM=="NA")] <- "NEA" 
        tmp <- paste(ecoreg_shp$ECO_ID, tmp)
        tmp[grep("NA",tmp)] <- NA 
        ecoreg_shp$EcoRealm <- tmp
        # remove Antarctica
        ecoreg_shp <- terra::subset(ecoreg_shp, ecoreg_shp$REALM != "AN")
        writeVector(ecoreg_shp, here::here(varsdir,"WWF_ecoregions/official/wwf_ecorealms.shp"), overwrite = TRUE)
        return(ecoreg_shp)
    }
}
get_wwf_raster <- function(varsdir){
    # check if raster file exists, if not create and save
    if(file.exists(here::here(varsdir,"WWF_ecoregions/official/EcoRealm.tif"))){
        ecoreg_raster <- terra::rast(here::here(varsdir,"WWF_ecoregions/official/EcoRealm.tif"))
        return(ecoreg_raster)
    } else {
        ecoreg_shp <- get_wwf_shp(varsdir)
        ecoreg_raster <- terra::rasterize(ecoreg_shp, mask.ras, field = "EcoRealm")
        # save
        writeRaster(ecoreg_raster, here::here(varsdir,"WWF_ecoregions/official/EcoRealm.tif"), overwrite = TRUE)
        return(ecoreg_raster)
    }
}
get_meow_shp <- function(varsdir){
    # create dir to save ecoregion 
    if(!dir.exists(here::here(varsdir,"MEOW"))){
        dir.create(here::here(varsdir,"MEOW"),recursive = TRUE)
    }
    
    # check if shape was downloaded and if not download it
    if(file.exists(here::here(varsdir,"MEOW","MEOW.zip"))){
        ecoreg_shp <- terra::vect(here::here(varsdir,"MEOW", "DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp"))
        return(ecoreg_shp)
    } else {
        cat("Downloading MEOW shape file")
        download.file(url = "https://datadownload-production.s3.amazonaws.com/WCMC036_MEOW_PPOW_2007_2012_v1.zip",
                      destfile = here::here(varsdir,"MEOW","MEOW.zip"))
        ## Unzip
        decompress_file(directory = here::here(varsdir,"MEOW"), 
                        file = here::here(varsdir,"MEOW","MEOW.zip"))
        ecoreg_shp <- terra::vect(here::here(varsdir,"MEOW", "DataPack-14_001_WCMC036_MEOW_PPOW_2007_2012_v1/01_Data/WCMC-036-MEOW-PPOW-2007-2012-NoCoast.shp"))
        return(ecoreg_shp)
    }
}
get_meow_raster <- function(varsdir){
    # create dir to save ecoregion 
    if(!dir.exists(here::here(varsdir,"MEOW"))){
        dir.create(here::here(varsdir,"MEOW"),recursive = TRUE)
    }
    # check if raster file exists, if not create and save
    if(file.exists(here::here(varsdir,"MEOW/MEOW.tif"))){
        ecoreg_raster <- terra::rast(here::here(varsdir,"MEOW/MEOW.tif"))
        return(ecoreg_raster)
    } else {
        ecoreg_shp <- get_meow_shp(varsdir)
        ecoreg_raster <- terra::rasterize(ecoreg_shp, mask.ras, field = "PROVINC")
        # Removes cells without environmental data
        ecoreg_raster[mask.ras==0] <- NA
        # save
        writeRaster(ecoreg_raster, here::here(varsdir,"MEOW/MEOW.tif"), overwrite = TRUE)
        return(ecoreg_raster)
    }
}
get_BA_shp <- function(varsdir,realm,PresAbs){
    if(realm=="Ter"){
        ecoreg_shp <- get_wwf_shp(varsdir)
        # At which ecoregions/realms the species occur?
        ecoreg <- terra::extract(ecoreg_shp, PresAbs)
        ecoreg <- na.omit(unique(ecoreg$EcoRealm))
        # Subset shape file to the ecoregions
        ecoreg_shp <- terra::subset(ecoreg_shp, ecoreg_shp$EcoRealm %in% ecoreg)
    }
    if(realm=="Mar"){
        ecoreg_shp <- get_meow_shp(varsdir)
        # At which ecoregions/realms the species occur?
        ecoreg <- terra::extract(ecoreg_shp, PresAbs)
        ecoreg <- na.omit(unique(ecoreg$PROVINC))
        # Subset shape file to the ecoregions
        ecoreg_shp <- terra::subset(ecoreg_shp, ecoreg_shp$PROVINC %in% ecoreg)
    }
    return(ecoreg_shp)
}
get_BA_raster <- function(varsdir,realm,PresAbs){
    if(realm=="Ter"){
        ecoreg_raster <- get_wwf_raster(varsdir)
        # Get a pres/abs BA
        tmpcells <- terra::cellFromXY(ecoreg_raster, xy=data.frame(PresAbs))
        # which ecoregions/realms the species occur?
        ecoreg <- ecoreg_raster[tmpcells]
        ecoreg <- as.character(ecoreg$EcoRealm)
        ecoreg <- na.omit(unique(ecoreg))
        # Crop rasters to the ecoregions/realms the species occur
        ecoreg_raster <- ecoreg_raster %in% ecoreg
    }
    if(realm=="Mar"){
        ecoreg_raster <- get_meow_raster(varsdir)
        # Get a pres/abs BA
        tmpcells <- terra::cellFromXY(ecoreg_raster, xy=data.frame(PresAbs))
        # which ecoregions/realms the species occur?
        ecoreg <- ecoreg_raster[tmpcells]
        ecoreg <- as.character(ecoreg$PROVINC)
        ecoreg <- na.omit(unique(ecoreg))
        # Crop rasters to the ecoregions/realms the species occur
        ecoreg_raster <- ecoreg_raster %in% ecoreg
    }
    # from TRUE/FALSE to 1/0 
    ecoreg_raster <- ifel(ecoreg_raster,1,0)
    names(ecoreg_raster) <- "BA"
    return(ecoreg_raster)
}