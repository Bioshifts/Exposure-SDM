mask_bios_BG <- function(bioclimatics_BG, output_BG){
    bioclimatics_BG <- lapply(bioclimatics_BG, function(x) {
        my_file <- strsplit(terra::sources(x),"/")[[1]]
        my_file <- here::here(output_BG,my_file[length(my_file)])
        # check if exists
        if(file.exists(my_file)){
            # if exists load in to see if it works
            test <- try(terra::rast(my_file), silent = TRUE)
            if(class(test) == "try-error"){
                # if error then mask and save
                terra::window(x) <- terra::ext(BA$shape_file)
                terra::mask(x, BA$shape_file, filename = my_file, overwrite = TRUE)
            } else {
                # if works, do nothing
                return(test)
            }
        } else {
            # if file does not exists, then mask
            terra::window(x) <- terra::ext(BA$shape_file)
            terra::mask(x, BA$shape_file, filename = my_file)
        }
    })
    names_bioclimatics_BG <- sapply(bioclimatics_BG, function(x){
        tmp <- strsplit(terra::sources(x),"/")[[1]]
        tmp <- tmp[length(tmp)]
        tmp <- strsplit(tmp,"_")[[1]]
        tmp <- tmp[length(tmp)]
        gsub(".tif","",tmp)
    })
    names(bioclimatics_BG) <- names_bioclimatics_BG
    return(bioclimatics_BG)
}

