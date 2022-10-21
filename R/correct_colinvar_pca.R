# based on correct_colinvar {flexsdm}

correct_colinvar_pca <- function(env_layer, proj, PresAbs, perc = 0.95){
    p1 <- terra::as.data.frame(env_layer, xy = FALSE, na.rm = TRUE)
    p <- stats::prcomp(p1, retx = TRUE, scale. = TRUE, center = TRUE)
    means <- p$center
    stds <- p$scale
    cof <- p$rotation
    cvar <- summary(p)$importance["Cumulative Proportion", ]
    naxis <- Position(function(x) {
        x >= perc
    }, cvar)
    cvar <- data.frame(cvar)
    p <- stats::prcomp(p1, retx = TRUE, scale. = TRUE, center = TRUE, rank. = naxis)
    env_layer <- terra::predict(env_layer, p)
    PresAbs <- predict(p, PresAbs)
    
    result <- list(PresAbs = PresAbs,
                   env_layer = env_layer, 
                   coefficients = data.frame(cof) %>% 
                       dplyr::tibble(variable = rownames(.), .), 
                   cumulative_variance = dplyr::tibble(PC = 1:nrow(cvar), 
                                                       cvar))
    
    if (!is.null(proj)) {
        scen <- terra::predict(proj, p)
        result$proj <- scen
    }
    
    return(result)
}
