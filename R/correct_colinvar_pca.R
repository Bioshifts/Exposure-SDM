# based on correct_colinvar {flexsdm}
options(warn=2)

correct_colinvar_pca <- function(
        sptogo,
        bios_BG=NULL, 
        bios_SA=NULL, 
        bios_PresAbs=NULL, 
        perc = 0.95, 
        env_cap = 10^4,
        output_dir = "",
        which_bioclimatics_BG = "all",
        check_if_PCA_model_exists = TRUE,
        N_cpus = NULL){
    
    # Calculate PCA
    if(check_if_PCA_model_exists){
        p <- try(readRDS(here::here(output_dir,"PCA_model.RDS")), silent = TRUE)
        coefficients <- try(read.csv(here::here(output_dir,"PCA_coefficients.csv")), silent = TRUE)
        cumulative_variance <- try(read.csv(here::here(output_dir,"PCA_cumulative_variance.csv")), silent = TRUE)
        
        error_test <- any(c(class(p)=="try-error",
                            class(coefficients)=="try-error",
                            class(cumulative_variance)=="try-error")) |
            any(c(is.null(p),
                  is.null(coefficients),
                  is.null(cumulative_variance)))
    } else {
        error_test = FALSE
    }
    
    if(error_test | !check_if_PCA_model_exists){
        
        cat("Extracting background env data\n")
        non_na_cells <- length(!is.na(values(bioclimatics_BG[[1]][[1]])))
        if(non_na_cells > env_cap){
            env_layer <- terra::spatSample(bioclimatics_BG[[1]], 
                                           size = env_cap, 
                                           na.rm=TRUE)
        } else {
            env_layer <- as.data.frame(bioclimatics_BG[[1]])
        }
        
        p <- stats::prcomp(env_layer, retx = TRUE, scale. = TRUE, center = TRUE)
        means <- p$center
        stds <- p$scale
        cof <- p$rotation
        cvar <- summary(p)$importance["Cumulative Proportion", ]
        naxis <- Position(function(x) {
            x >= perc
        }, cvar)
        cvar <- data.frame(cvar)
        
        coefficients = data.frame(cof) %>% 
            dplyr::tibble(variable = rownames(.), .)
        
        cumulative_variance = dplyr::tibble(PC = 1:nrow(cvar), 
                                            cvar)
        
        p <- stats::prcomp(env_layer, retx = TRUE, scale. = TRUE, center = TRUE, rank. = naxis)
    }
    
    
    # Predict to the PresAbs
    if (!is.null(bios_PresAbs)) {
        bios_PresAbs_pred <- predict(p, bios_PresAbs)
        bios_PresAbs <- data.table(data.frame(bios_PresAbs[,c("x","y","pa")], bios_PresAbs_pred))
        
        if(!output_dir==""){
            qs::qsave(bios_PresAbs,
                      here::here(output_dir,paste0(sptogo,"_PresAbs_PC.qs")))
        } 
        bios_PresAbs <- bios_PresAbs
    } else {
        bios_PresAbs <- NULL
    }
    
    
    # Predict to the BG
    if (!is.null(bios_BG)) {
        BG_PC_dir <- here::here(output_dir,"BG_PC")
        if(!dir.exists(BG_PC_dir)){
            dir.create(BG_PC_dir)
        }
        which_ones <- ifelse(which_bioclimatics_BG=="all",1:length(bios_BG),which_bioclimatics_BG)
        bios_BG <- bios_BG[which_ones]
        tmp <- lapply(bios_BG, function(x){
            my_file <- strsplit(terra::sources(x),"/")[[1]]
            my_file <- my_file[length(my_file)]
            test <- try(rast(here::here(BG_PC_dir,my_file)),silent = TRUE)
            if(class(test)=="try-error"){
                terra::predict(x, 
                               model = p, 
                               filename = here::here(BG_PC_dir,my_file),
                               overwrite = TRUE)
            } else {
                test
            }
        })
        bios_BG_PC <- lapply(bios_BG, function(x){
            my_file <- strsplit(terra::sources(x),"/")[[1]]
            my_file <- my_file[length(my_file)]
            rast(here::here(BG_PC_dir,my_file))
        })
        
        
    } else {
        bios_BG_PC <- NULL
    }
    
    
    # Predict to the SA
    if (!is.null(bios_SA)) {
        SA_PC_dir <- here::here(output_dir,"SA_PC")
        if(!dir.exists(SA_PC_dir)){
            dir.create(SA_PC_dir)
        }
        tmp <- lapply(bios_SA, function(area){
            parallel::mclapply(area, function(x){
                my_file <- strsplit(terra::sources(x),"/")[[1]]
                my_file <- my_file[length(my_file)]
                test <- try(rast(here::here(SA_PC_dir,my_file)),silent = TRUE)
                if(class(test)=="try-error"){
                    terra::predict(x, 
                                   model = p, 
                                   filename = here::here(SA_PC_dir,my_file), 
                                   overwrite = TRUE)
                } 
            }, mc.cores = 10)
        })
        bios_SA_PC <- lapply(bios_SA, function(y){
            lapply(y, function(x){
                my_file <- strsplit(terra::sources(x),"/")[[1]]
                my_file <- my_file[length(my_file)]
                rast(here::here(SA_PC_dir,my_file))
            })
        })
        
        
    } else {
        bios_SA_PC <- NULL
    }
    
    
    result <- list(coefficients = coefficients, 
                   cumulative_variance = cumulative_variance,
                   bios_PresAbs = bios_PresAbs,
                   bios_BG = bios_BG_PC,
                   bios_SA = bios_SA_PC,
                   PCA_model = p)
    return(result)
    
}

PCA_env <- function(
        sptogo,
        bioclimatics_BG,
        bioclimatics_SA,
        which_bioclimatics_BG = "all",
        env_cap = 10^4, # Maximum size of dataset used to calculate the PCA
        output_dir,
        shift_info,
        PresAbs,
        check_if_PCA_model_exists = TRUE,
        check_if_data_exists = TRUE,
        N_cpus = NULL
){
    
    if(check_if_data_exists){
        
        # 1st check if has data
        ## PCA model
        if(file.exists(here::here(output_dir,"PCA_model.RDS"))){
            my_test_PCAmodel <- NULL
            cat("PCA model exists\n")
        } else{
            cat("Calculating PCA model\n")
            my_test_PCAmodel <- NA
        }
        ## all bios BG
        all_bios_BG <- list.files(here::here(output_dir,"BG"))
        ## all bios BG PC
        all_bios_BG_PC <- list.files(here::here(output_dir,"BG_PC"))
        ## test if all bios BG PC open
        all_bios_BG_PC_work <- sapply(list.files(here::here(output_dir,"BG_PC"), full.names = TRUE),
                                      function(x){
                                          test <- try(terra::rast(x))
                                          !class(test)=="try-error"
                                      })
        
        ## all bios SA
        all_bios_SA <- sapply(1:nrow(shift_info), function(i) {
            paste0(shift_info$ID[i], " bios ", round(shift_info$Start[i],0):round(shift_info$End[i],0), ".tif")
        })
        all_bios_SA <- as.vector(unlist(all_bios_SA))
        ## all bios SA PC 
        all_bios_SA_PC <- list.files(here::here(output_dir,"SA_PC"))
        ## test if all bios SA PC open
        all_bios_SA_PC_work <- sapply(list.files(here::here(output_dir,"SA_PC"), full.names = TRUE),
                                      function(x){
                                          test <- try(terra::rast(x))
                                          !class(test)=="try-error"
                                      })
        
        ### My test
        # PresAbs
        if(file.exists(here::here(output_dir,paste0(sptogo,"_PresAbs_PC.qs")))){
            my_test_PresAbs <- NULL
            cat("PCs exist for PresAbs\n")
        } else{
            my_test_PresAbs <- PresAbs
            cat("Calculating PCs for PresAbs\n")
        }
        
        # BG
        if(all(all_bios_BG %in% all_bios_BG_PC) & all(all_bios_BG_PC_work)){
            my_test_BG <- NULL
            cat("PCs exist for BG\n")
        } else {
            my_test_BG <- bioclimatics_BG
            cat("Calculating PCs for BG\n")
        }
        
        # SA
        if(all(all_bios_SA %in% all_bios_SA_PC) & all(all_bios_SA_PC_work)){
            my_test_SA <- NULL
            cat("PCs exist for SA\n")
        } else {
            my_test_SA <- bioclimatics_SA
            cat("Calculating PCs for SA\n")
        }
        
        
        test <- any(!c(is.null(my_test_PCAmodel),is.null(my_test_BG),is.null(my_test_SA),is.null(my_test_PresAbs)))
        
    } else {
        
        my_test_PCAmodel <- NA
        my_test_BG <- bioclimatics_BG
        my_test_SA <- bioclimatics_SA
        my_test_PresAbs <- PresAbs
        
        test <- any(!c(is.null(my_test_PCAmodel),is.null(my_test_BG),is.null(my_test_SA),is.null(my_test_PresAbs))) 
        
    }
    
    #########################
    # Get PCs
    if(test){
        
        cat("Calculating PCs\n")
        new_data <- correct_colinvar_pca(
            sptogo = sptogo,
            bios_BG = my_test_BG, 
            bios_SA = my_test_SA, 
            bios_PresAbs = my_test_PresAbs,
            env_cap = env_cap,
            which_bioclimatics_BG = which_bioclimatics_BG,
            output_dir = output_dir,
            check_if_PCA_model_exists = check_if_PCA_model_exists)
        
    } else {
        new_data <- list(bios_BG = my_test_BG, 
                         bios_SA = my_test_SA, 
                         bios_PresAbs = my_test_PresAbs,
                         p = readRDS(here::here(output_dir,"PCA_model.RDS")),
                         coefficients = read.csv(here::here(output_dir,"PCA_coefficients.csv")),
                         cumulative_variance = read.csv(here::here(output_dir,"PCA_cumulative_variance.csv")))
        
    }
    
    
    #########################
    # Update data
    # BG
    if(is.null(new_data$bios_BG)){
        new_data$bios_BG <- lapply(all_bios_BG_PC, function(x){
            terra::rast(here::here(output_dir,"BG_PC",x))
        })
        names_bioclimatics_BG <- sapply(new_data$bios_BG, function(x){
            tmp <- strsplit(terra::sources(x),"/")[[1]]
            tmp <- tmp[length(tmp)]
            tmp <- strsplit(tmp,"_")[[1]]
            tmp <- tmp[length(tmp)]
            gsub(".tif","",tmp)
        })
        names(new_data$bios_BG) <- names_bioclimatics_BG
    } 
    
    # SA
    if(is.null(new_data$bios_SA)){
        new_data$bios_SA <- lapply(shift_info$ID, function(x) {
            tmp_names <- list.files(here::here(output_dir,"SA_PC"), pattern = x)
            tmp_names <- gsub(paste0(c(x,"bios",".tif"," "),collapse = "|"),"",tmp_names)
            tmp <- list.files(here::here(output_dir,"SA_PC"), pattern = x, full.names = TRUE)
            tmp <- lapply(tmp, terra::rast)
            names(tmp) <- tmp_names
            return(tmp)
        })
        names(new_data$bios_SA) <- shift_info$ID
    } 
    
    # PresAbs
    if(is.null(my_test_PresAbs)){
        new_data$bios_PresAbs <- qs::qread(here::here(output_dir,paste0(sptogo,"_PresAbs_PC.qs")))
    }
    
    
    #########################
    # return results
    return(new_data)    
    
}



options(warn=1)