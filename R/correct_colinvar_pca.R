# based on correct_colinvar {flexsdm}

correct_colinvar_pca <- function(env_layer, proj, PresAbs, perc = 0.95){
    
    p <- stats::prcomp(env_layer, retx = TRUE, scale. = TRUE, center = TRUE)
    means <- p$center
    stds <- p$scale
    cof <- p$rotation
    cvar <- summary(p)$importance["Cumulative Proportion", ]
    naxis <- Position(function(x) {
        x >= perc
    }, cvar)
    cvar <- data.frame(cvar)
    
    p <- stats::prcomp(env_layer, retx = TRUE, scale. = TRUE, center = TRUE, rank. = naxis)
    
    PresAbs <- predict(p, PresAbs)
    
    result <- list(PresAbs = PresAbs,
                   coefficients = data.frame(cof) %>% 
                       dplyr::tibble(variable = rownames(.), .), 
                   cumulative_variance = dplyr::tibble(PC = 1:nrow(cvar), 
                                                       cvar))
    
    if (!is.null(proj)) {
        scen <- lapply(proj, terra::predict, model = p)
            
        result$proj <- scen
    }
    
    return(result)
    
}