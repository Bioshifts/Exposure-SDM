
range_shift <- function(sui_t1,
                        sui_t2,
                        var_get = "LAT", 
                        info, 
                        get_mean = FALSE,
                        probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)){
    
    if(terra::ncell(sui_t1) > 100000){
        sui_t1 <- terra::spatSample(sui_t1, 100000, na.rm = TRUE, xy = TRUE)
        sui_t2 <- terra::spatSample(sui_t2, 100000, na.rm = TRUE, xy = TRUE)
    } else {
        sui_t1 <- terra::as.data.frame(sui_t1, na.rm = TRUE, xy = TRUE)
        sui_t2 <- terra::as.data.frame(sui_t2, na.rm = TRUE, xy = TRUE)
    }
    
    # calc quantile weighted type for SA in t1
    if(var_get=="LAT"){
        var2go <- sui_t1$y
    }
    wti1 <- wtd.quantile(var2go, sui_t1[,3], probs = probs) 
    names(wti1) <- gsub(" ","",names(wti1))
    wti1 <- as.data.frame(wti1)
    wti1$quant <- rownames(wti1)
    if(get_mean){
        # calc mean type for t1
        tmp <- data.frame(as.numeric(wtd.mean(var2go, sui_t1[,3])), "mean")
        names(tmp) <- names(wti1)
        wti1 <- rbind(wti1, tmp)
    }
    rownames(wti1) <- NULL
    
    # calc quantile weighted type for SA in t2
    wti2 <- wtd.quantile(var2go, sui_t2[,3], probs = probs) 
    names(wti2) <- gsub(" ","",names(wti2))
    wti2 <- as.data.frame(wti2)
    wti2$quant <- rownames(wti2)
    if(get_mean){
        # calc mean type for t2
        tmp <- data.frame(as.numeric(wtd.mean(var2go, sui_t2[,3])), "mean")
        names(tmp) <- names(wti2)
        wti2 <- rbind(wti2, tmp)
    }
    rownames(wti2) <- NULL
    
    # calc shift t1 t2
    shift_i <- data.frame(quant = wti1$quant, w_t1 = wti1$wti1, w_t2 = wti2$wti2)
    # when shifts are in the south hemisphere, use absolute values...
    south <- which(shift_i$w_t1 < 0 & shift_i$w_t2 < 0)
    shift_i$shift[south] <- abs(shift_i$w_t2[south]) - abs(shift_i$w_t1[south])
    # when shifts are in the north hemisphere, or cross hemispheres, use raw values...
    north <- which(!shift_i$w_t1 < 0 & !shift_i$w_t2 < 0)
    shift_i$shift[north] <- shift_i$w_t2[north] - shift_i$w_t1[north]
    
    # add info
    shift_i$ID <- info$ID
    shift_i$Type <- var_get
    shift_i$START <- info$START
    shift_i$END <- info$END
    shift_i$Species <- info$Species
    
    return(shift_i)
}

range01raster <- function(x){
    mm <- terra::minmax(x)
    tmp <- (x-mm[1])/(mm[2]-mm[1])
    return(tmp)
}

exposure_shift <- function(back_t1, back_t2, sui_t1, sui_t2){
    
    # calc mean type for t1
    mp <- data.frame(as.numeric(wtd.mean(var_t, sui_t1)), "mean")
    names(tmp) <- names(wti1)
    wti1 <- rbind(wti1, tmp)
    rownames(wti1) <- NULL
    
    # calc mean type for t1
    tmp <- data.frame(as.numeric(wtd.mean(var_t, sui_t2)), "mean")
    names(tmp) <- names(wti2)
    wti2 <- rbind(wti2, tmp)
    rownames(wti2) <- NULL
    
    # calc shift t1 t2
    shift_i <- data.frame(quant = wti1$quant, wti1 = wti1$wti1, wti2 = wti2$wti2)
    shift_i$shift <- shift_i$wti2 - shift_i$wti1
    
    return(shift_i)
}


sens_shift <- function(SA, BG, type, probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),variable=NULL){
    
    # get type SA
    tSA <- terra::as.data.frame(SA, xy = TRUE)
    if(type[i] == "LAT"){ tSA <- tSA$y }
    if(type[i] == "ELE"){ 
        # crop to the SA 
    }
    # get type BG
    tBG <- terra::as.data.frame(BG, xy = TRUE)
    if(type[i] == "LAT"){ tBG <- tBG$y }
    if(type[i] == "ELE"){ 
        # crop to the BG 
    }
    
    # difference in type
    tSA <- quantile(tSA, probs)
    tGB <- quantile(tBG, probs)
    return(data.frame(SA = quantile(tSA, probs), BG = quantile(tBG, probs), Quantile = names(quantile(tBG, probs))))
}


range_pos <- function(SA,BG){
    
    MaxSAt1 = SA$wti1[which(SA$quant=="95%")]
    MinSAt1 = SA$wti1[which(SA$quant=="5%")]
    RangeSAt1 = MaxSAt1-MinSAt1
    MaxBGt1 = BG$wti1[which(SA$quant=="95%")]
    MinBGt1 = BG$wti1[which(SA$quant=="5%")]
    RangeBGt1 = MaxBGt1-MinBGt1
    Fillingt1 = RangeSAt1/RangeBGt1
    
    MaxSAt2 = SA$wti2[which(SA$quant=="95%")]
    MinSAt2 = SA$wti2[which(SA$quant=="5%")]
    RangeSAt2 = MaxSAt2-MinSAt2
    MaxBGt2 = BG$wti2[which(SA$quant=="95%")]
    MinBGt2 = BG$wti2[which(SA$quant=="5%")]
    RangeBGt2 = MaxBGt2-MinBGt2
    Fillingt2 = RangeSAt2/RangeBGt2
    
    return(data.frame(MaxSAt1, MinSAt1, RangeSAt1, MaxBGt1, MinBGt1, RangeBGt1, Fillingt1,
                      MaxSAt2, MinSAt2, RangeSAt2, MaxBGt2, MinBGt2, RangeBGt2, Fillingt2))
}


