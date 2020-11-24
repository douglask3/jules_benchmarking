openObs <- function(obs, Layers, scale, modEG, ..., layer = NULL) { 
    
    tempFile = paste(c('temp/', filename.noPath(obs, TRUE),
                     range(Layers), ..., scale, '.nc'),
                     collapse = '-')
    cat("\nopening:", obs)
    cat("\n\tinto:", tempFile)
    
    if (file.exists(tempFile)) obs = brick(tempFile)
    else {
        obs = brick(obs)
        obs = obs[[Layers]]
        obs = datConvert(obs, scale, modEG = modEG[[1]], ...)

        #obs[is.na(modEG)] = NaN
        #browser()
        writeRaster(obs, file = tempFile, overwrite = TRUE)
    }
    if (!is.null(layer)) obs = obs[[layer]]
    if (nlayers(obs) == 1) obs = obs[[1]]
    return(obs)
}
