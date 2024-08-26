# Get file size before downloading
DownloadSizeChelsa <- function(myvar, period){
    
    download_size <- function(url) as.numeric(httr::HEAD(url)$headers$`content-length`)
    
    yr <- as.numeric(strsplit(period,"_")[[1]][2])
    mnth <- strsplit(period,"_")[[1]][1]
    mnth_fx <- as.numeric(strsplit(period,"_")[[1]][1])
    
    if(yr < 1980){ # CHELSA CRUTS
        link <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_cruts/"
        link <- paste0(link,myvar,"/")
        ver_file <- "V.1.0.tif"
        mydataset <- "CHELSAcruts_"
        myfile_link <- paste0(mydataset, paste(myvar, mnth_fx, yr, ver_file, sep = "_"))
    } else { # CHELSA
        link <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/"
        if(myvar == "tmax"){myvar_new <- "tasmax"}
        if(myvar == "tmin"){myvar_new <- "tasmin"}
        if(myvar == "prec"){myvar_new <- "pr"}
        link <- paste0(link,myvar_new,"/")
        ver_file <- "V.2.1.tif"
        mydataset <- "CHELSA_"
        myfile_link <- paste0(mydataset, paste(myvar_new, mnth, yr, ver_file, sep = "_"))
    }
    
    mylink <- paste0(link, myfile_link)
    
    download_size(mylink)
}

# Get Chelsa data
getChelsa <- function(myvar, 
                      period, 
                      dest.folder, 
                      silent = TRUE, # avoids errors to stop the download in a loop
                      quiet = FALSE, # show progress of downloads
                      over.write = FALSE){
    
    if(!is.logical(quiet)){
        cat("quiet should be TRUE or FALSE")
        stop()
    }
    
    new.dir <- here::here(dest.folder, myvar)
    if(!dir.exists(new.dir)){
        dir.create(new.dir,recursive = T)
    }
    
    yr <- as.numeric(strsplit(period,"_")[[1]][2])
    mnth <- strsplit(period,"_")[[1]][1]
    mnth_fx <- as.numeric(strsplit(period,"_")[[1]][1])
    
    link <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V2/GLOBAL/monthly/"
    if(myvar == "tmax"){myvar_new <- "tasmax"}
    if(myvar == "tmin"){myvar_new <- "tasmin"}
    if(myvar == "prec"){myvar_new <- "pr"}
    link <- paste0(link,myvar_new,"/")
    ver_file <- "V.2.1.tif"
    mydataset <- "CHELSA_"
    myfile_save <- paste0(paste(myvar, mnth, yr, sep = "_"),".tif")
    myfile_link <- paste0(mydataset, paste(myvar_new, mnth, yr, ver_file, sep = "_"))
    
    mylink <- paste0(link, myfile_link)
    savefile <- here::here(dest.folder, myvar, myfile_save)
    
    if(!over.write){
        over.write <- !file.exists(savefile)
    }
    if(over.write){
        if(silent){
            try(download.file(url = mylink,
                              destfile = savefile,
                              quiet = quiet,
                              cacheOK = FALSE),
                silent = TRUE)
        } else {
            download.file(url = mylink,
                          destfile = savefile,
                          quiet = quiet,
                          cacheOK = FALSE)
        }
        
    }
    
}

getChelsa_cruts <- function(myvar, 
                            period, 
                            dest.folder, 
                            silent = TRUE, # avoids errors to stop the download in a loop
                            quiet = FALSE, # show progress of downloads
                            over.write = FALSE){
    
    if(!is.logical(quiet)){
        cat("quiet should be TRUE or FALSE")
        stop()
    }
    
    new.dir <- here::here(dest.folder, myvar)
    if(!dir.exists(new.dir)){
        dir.create(new.dir,recursive = T)
    }
    
    yr <- as.numeric(strsplit(period,"_")[[1]][2])
    mnth <- strsplit(period,"_")[[1]][1]
    mnth_fx <- as.numeric(strsplit(period,"_")[[1]][1])
    
    link <- "https://os.zhdk.cloud.switch.ch/envicloud/chelsa/chelsa_V1/chelsa_cruts/"
    link <- paste0(link,myvar,"/")
    ver_file <- "V.1.0.tif"
    mydataset <- "CHELSAcruts_"
    myfile_save <- paste0(paste(myvar, mnth, yr, sep = "_"),".tif")
    myfile_link <- paste0(mydataset, paste(myvar, mnth_fx, yr, ver_file, sep = "_"))
    
    mylink <- paste0(link, myfile_link)
    savefile <- here::here(dest.folder, myvar, myfile_save)
    
    if(!over.write){
        over.write <- !file.exists(savefile)
    }
    if(over.write){
        if(silent){
            try(download.file(url = mylink,
                              destfile = savefile,
                              quiet = quiet,
                              cacheOK = FALSE),
                silent = TRUE)
        } else {
            download.file(url = mylink,
                          destfile = savefile,
                          quiet = quiet,
                          cacheOK = FALSE)
        }
        
    }
    
}