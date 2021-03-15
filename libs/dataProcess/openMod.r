openMod <- function(mod, dir, varName, years, modScale, ..., layer = NULL) {
    tempFile = paste(c('temp/', paste(strsplit(mod, '/')[[1]], collapse = '_'),
                     varName, range(years), ..., modScale, '.nc'),
                     collapse = '-')
    cat("\nopening:", mod)
    cat("\n\tinto:", tempFile, "\n")
    
    if (file.exists(tempFile)) dat = brick(tempFile)
    else {
        files = list.files(paste0(dir, mod, '/'), full.names = TRUE)
        fyr = substr(files, nchar(files)-6, nchar(files)-3)
   
        files = files[apply( sapply(fyr, '==', years), 2, any)]
        
        varNames = strsplit(varName, ';')[[1]]
        openVar <- function(varName) {
            dat = layer.apply(files, process.jules.file, 1, varName)
            
            dat = datConvert(dat, modScale, ...)
            dat = convert_pacific_centric_2_regular(dat)
        }
        dati = lapply(varNames, openVar)
        dat = dati[[1]]
        if (length(dati) > 1) for (d in dat[-1]) dat = dat + d
        
        dat = writeRaster(dat, file = tempFile, overwrite = TRUE)
    }
    if (!is.null(layer)) dat = dat[[layer]]
    if (nlayers(dat) == 1) dat = dat[[1]]
    
    return(dat)
}
