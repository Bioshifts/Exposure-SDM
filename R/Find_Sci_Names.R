# search species names in GBIF backbone taxonomy
# can provide kingdom information to help with searches (optional)
get_rgbif <- function(sp_name, kingdom=NULL, return_species_rank=TRUE, return_found_only = TRUE){
    
    cl <- makeCluster(detectCores()-2)
    clusterExport(cl, c("sp_name","kingdom"), envir = environment())
    
    if(!is.null(kingdom)){
        tmp <- pblapply(1:length(sp_name), function(i){
            x = rgbif::name_backbone(
                name = sp_name[i],
                kingdom = kingdom[i])
            
            if(x$matchType=="HIGHERRANK"){
                x2 = rgbif::name_backbone(
                    name = sp_name[i],
                    kingdom = kingdom[i],
                    verbose = TRUE)
                
                if(any(x2$matchType=="EXACT")){
                    x2 = x2[which(x2$matchType=="EXACT"),]
                    x2 = x2[1,]
                    x = x2
                }
            }
            if(any(x$matchType=="EXACT")){
                x = x[which(x$matchType=="EXACT"),]
                x = x[1,]
            }
            
            x$requested_name = sp_name[i]
            
            return(x)
        }, cl = cl) 
    } else {
        tmp <- pblapply(1:length(sp_name), function(i){
            x = rgbif::name_backbone(
                name = sp_name[i])
            
            if(x$matchType=="HIGHERRANK"){
                x2 = rgbif::name_backbone(
                    name = sp_name[i],
                    verbose = TRUE)
                
                if(any(x2$matchType=="EXACT")){
                    x2 = x2[which(x2$matchType=="EXACT"),]
                    x2 = x2[1,]
                    x = x2
                }
            }
            if(any(x$matchType=="EXACT")){
                x = x[which(x$matchType=="EXACT"),]
                x = x[1,]
            }
            
            x$requested_name = sp_name[i]
            
            return(x)
        }, cl = cl) 
    }
    
    stopCluster(cl)
    
    tmp <- rbindlist(tmp, fill = TRUE)
    
    # fix db_code
    if(!is.null(tmp$acceptedUsageKey)){
        tmp$db_code <- tmp$acceptedUsageKey
    } else {
        tmp$db_code <- tmp$usageKey
    }
    tmp$db_code[which(is.na(tmp$db_code))] <- tmp$usageKey[which(is.na(tmp$db_code))]
    tmp$db_code <- paste0("GBIF:",tmp$db_code)
    
    tmp$scientificName = tmp$species
    tmp$db <- "GBIF"
    
    if(return_species_rank){
        # Remove which rank is not species
        pos <- which(!tmp$rank=="SPECIES")
        tmp[pos,
            c("kingdom", "phylum", "class", "order", "family", "scientificName", "db_code", "db")] <- NA
    }
    
    if(return_found_only){
        found_position <- !is.na(tmp$db_code)
        tmp <- tmp[which(found_position),]
    }
    
    # organize columns
    tmp <- 
        tmp[,
            c("kingdom", "phylum", "class", "order", "family", "scientificName", "db_code", "db","requested_name")]
    
    return(tmp)
    
}

get_gbif_bdc <- function(sp_name, kingdom=NULL, return_species_rank=TRUE, return_found_only = TRUE){
    
    
    
    if(!is.null(kingdom)){
        
        kingdom_i <- unique(kingdom)
        tmp <- list()
        
        for(i in 1:length(kingdom_i)){
            
            sps_to_go <- sp_name[which(kingdom == kingdom_i[i])]
            
            tmp[[i]] <- 
                bdc::bdc_query_names_taxadb(sci_name = sps_to_go,
                                            rank_name = kingdom_i[i],
                                            rank = "kingdom",
                                            replace_synonyms = TRUE,
                                            suggest_names = FALSE,
                                            db = "gbif",
                                            parallel = TRUE,
                                            ncores = detectCores()-2)
            
            tmp <- rbindlist(tmp)
        }
    } else {
        
        tmp <- 
            bdc::bdc_query_names_taxadb(sci_name = sp_name,
                                        replace_synonyms = TRUE,
                                        suggest_names = FALSE,
                                        db = "gbif",
                                        parallel = TRUE,
                                        ncores = detectCores()-2)
        
    }
    
    
    tmp$requested_name <- tmp$original_search
    tmp$db_code <- tmp$taxonID
    tmp$db <- "GBIF"
    
    if(return_species_rank){
        # Remove which rank is not species
        pos <- which(!tmp$taxonRank=="species" | is.na(tmp$taxonRank))
        tmp[pos,
            c("kingdom", "phylum", "class", "order", "family", "scientificName", "db_code", "db")] <- NA
    }
    
    if(return_found_only){
        found_position <- !is.na(tmp$db_code) & !is.na(tmp$acceptedUsageKey)
        tmp <- tmp[which(found_position),]
    }
    
    # organize columns
    tmp <- 
        tmp[,
            c("kingdom", "phylum", "class", "order", "family", "scientificName", "db_code", "db","requested_name")]
    
    return(tmp)
    
}

get_taxadb <- function(sp_name, providers = c("gbif","itis","ncbi"), return_species_rank=TRUE, return_found_only = TRUE){
    
    sp_name = sp_name
    
    new_search <- list()
    
    for(i in 1:length(providers)){
        cat("\nSearching in", providers[i])
        tmp <- taxadb::filter_name(sp_name,
                                   provider = providers[i])
        tmp$requested_name <- tmp$scientificName
        tmp$db <- providers[i]
        if(providers[i] %in% c("col","itis","gbif")){
            tmp$scientificName <- paste(tmp$genus,tmp$specificEpithet)
        } else {
            tmp$scientificName <- tmp$specificEpithet
        }
        new_search[[i]] <- tmp
    }
    
    new_search <- rbindlist(new_search,fill = TRUE)
    new_search <- new_search[-which(duplicated(new_search$requested_name)),]
    new_search$db <- toupper(new_search$db)
    new_search$db_code <- new_search$acceptedNameUsageID
    
    # organize columns
    new_search <- 
        new_search[,
                   c("kingdom", "phylum", "class", "order", "family", "scientificName", "db_code", "db","requested_name")]
    
    return(new_search)
    
}


#' Harmonizing taxon names based on multiple sources giving priority to GBIF.
#'
#' @param sp_name Vector of species names
#' @param kingdom optional. Vector of kingdom names to help resolve species names
#' @param providers list of taxonomy providers for performing additional searches
#' @param return_species_rank logical. Should return only taxa identified to the species level? Default = TRUE.
#' @param return_found_only logical. Should return only taxa with found accepted names? Default = TRUE.
#'
#' @return data.frame
#' 
Find_Sci_Names <- function(sp_name, kingdom=NULL, providers = c("gbif","itis","ncbi"), return_species_rank=TRUE, return_found_only = TRUE){
    
    tofind <- data.frame(matrix(nrow = length(sp_name), ncol = 8))
    names(tofind) = c("scientificName", "kingdom", "phylum", "class", "order", "family", "db", "db_code")
    
    tofind <- data.frame(species = sp_name, tofind)
    
    tofind <- tofind %>%
        mutate(across(everything(), as.character))
    
    if(!is.null(kingdom)){
        tofind$kingdom <- kingdom
    }
    
    # using rgbif package
    cat("\n\n----------------\nSearching", length(sp_name), "species using rgbif package",
        "\n----------------\n\n")
    
    rgbif_names <- get_rgbif(sp_name = sp_name, 
                             kingdom = kingdom,
                             return_species_rank = return_species_rank,
                             return_found_only = return_found_only)
    
    rgbif_names$method <- "rgbif"
    
    ## Feed
    for(i in 1:nrow(rgbif_names)){
        tofill <- unique(which(tofind$species == rgbif_names$requested_name[i]))
        tofind$scientificName[tofill] <- rgbif_names$scientificName[i]
        tofind$kingdom[tofill] <- rgbif_names$kingdom[i]
        tofind$phylum[tofill] <-rgbif_names$phylum[i]
        tofind$class[tofill] <- rgbif_names$class[i]
        tofind$order[tofill] <- rgbif_names$order[i]
        tofind$family[tofill] <- rgbif_names$family[i]
        tofind$db[tofill] <- rgbif_names$db[i]
        tofind$db_code[tofill] <- rgbif_names$db_code[i]
    }
    
    #############################
    # using bdc package
    
    ## Found vs not found 
    not_found_position <- is.na(tofind$scientificName)
    tofind_ <- tofind[not_found_position,]
    
    cat("\n\n----------------\nSearching", length(tofind_$species), "species using bdc package",
        "\n----------------\n\n")
    
    if(nrow(tofind_)>0){
        if(is.null(kingdom)){
            gbif_names_bdc <- get_gbif_bdc(sp_name = tofind_$species, 
                                           return_species_rank = return_species_rank,
                                           return_found_only = return_found_only)
            
            gbif_names_bdc$method <- "bdc"
        } else {
            gbif_names_bdc <- get_gbif_bdc(sp_name = tofind_$species, 
                                           kingdom = tofind_$kingdom,
                                           return_species_rank = return_species_rank,
                                           return_found_only = return_found_only)
            
            gbif_names_bdc$method <- "bdc"
        }
        
        
    } else {
        gbif_names_bdc <- data.frame(matrix(ncol = ncol(rgbif_names)))
        colnames(gbif_names_bdc) <- colnames(rgbif_names)
    }
    
    
    ## Feed
    for(i in 1:nrow(gbif_names_bdc)){
        tofill <- unique(which(tofind$species == gbif_names_bdc$requested_name[i]))
        tofind$scientificName[tofill] <- gbif_names_bdc$scientificName[i]
        tofind$kingdom[tofill] <- gbif_names_bdc$kingdom[i]
        tofind$phylum[tofill] <-gbif_names_bdc$phylum[i]
        tofind$class[tofill] <- gbif_names_bdc$class[i]
        tofind$order[tofill] <- gbif_names_bdc$order[i]
        tofind$family[tofill] <- gbif_names_bdc$family[i]
        tofind$db[tofill] <- gbif_names_bdc$db[i]
        tofind$db_code[tofill] <- gbif_names_bdc$db_code[i]
    }
    
    #############################
    # Use taxadb
    
    # Looking through other databases
    ## Found vs not found 
    not_found_position <- is.na(tofind$scientificName)
    tofind_ <- tofind[not_found_position,]
    
    cat("\n\n----------------\nSearching", length(tofind_$species), "species using taxadb package",
        "\n----------------\n\n")
    
    if(nrow(tofind_)>0){
        
        names_others <- get_taxadb(sp_name = tofind_$species,
                                   providers = providers,
                                   return_species_rank = return_species_rank,
                                   return_found_only = return_found_only)
        
        names_others$method <- "taxadb"
    } else {
        names_others <- data.frame(matrix(ncol = ncol(rgbif_names)))
        colnames(names_others) <- colnames(rgbif_names)
    }
    
    
    ## Feed
    for(i in 1:nrow(names_others)){
        tofill <- unique(which(tofind$species == names_others$requested_name[i]))
        tofind$scientificName[tofill] <- names_others$scientificName[i]
        tofind$kingdom[tofill] <- names_others$kingdom[i]
        tofind$phylum[tofill] <-names_others$phylum[i]
        tofind$class[tofill] <- names_others$class[i]
        tofind$order[tofill] <- names_others$order[i]
        tofind$family[tofill] <- names_others$family[i]
        tofind$db[tofill] <- names_others$db[i]
        tofind$db_code[tofill] <- names_others$db_code[i]
    }
    
    #############################
    ## Group all
    all_found <- rbindlist(list(data.frame(rgbif_names),
                                data.frame(gbif_names_bdc),
                                data.frame(names_others)), 
                           fill = TRUE)
    
    if(any(is.na(all_found$scientificName))){
        all_found <- all_found[-which(is.na(all_found$scientificName)),]
    }
    
    
    #############################
    # GBIF priority
    cat("\n\n----------------\nPrioritizing GBIF\n----------------\n")
    
    # Use species names found in other databases, and look for names in GBIF backbone
    ## Found vs not found 
    not_found_position <- !all_found$db == "GBIF"
    tofind_ <- all_found[which(not_found_position),]
    
    if(nrow(tofind_)>0){
        
        gbif_names2 <- get_rgbif(sp_name = tofind_$scientificName, 
                                 kingdom = tofind_$kingdom,
                                 return_species_rank = return_species_rank,
                                 return_found_only = return_found_only)
        
        if(nrow(gbif_names2)>0){
            ## Feed
            for(i in 1:nrow(gbif_names2)){
                tofill <- unique(which(tofind$species == gbif_names2$requested_name[i]))
                tofind$scientificName[tofill] <- gbif_names2$scientificName[i]
                tofind$kingdom[tofill] <- gbif_names2$kingdom[i]
                tofind$phylum[tofill] <-gbif_names2$phylum[i]
                tofind$class[tofill] <- gbif_names2$class[i]
                tofind$order[tofill] <- gbif_names2$order[i]
                tofind$family[tofill] <- gbif_names2$family[i]
                tofind$db[tofill] <- gbif_names2$db[i]
                tofind$db_code[tofill] <- gbif_names2$db_code[i]
            }
        }
    }
    
    # Summary
    res <- data.frame(table(all_found$db))
    names(res) <- c("db","N")
    cat("\n\n----------------\n Summary \n----------------\n",
        "N taxa:",length(sp_name),
        "N taxa found:",knitr::kable(res),
        "N taxa not found:", length(sp_name)-nrow(all_found),
        sep = "\n")
    
    return(all_found)
}

