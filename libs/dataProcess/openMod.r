openMod <- function(mod, dir, varName, years, modScale, ..., stream = NULL, levels = 1, layer = NULL) {
    
    if (!exists("extent") || is.null(extent)) extent = c(-180, 180, -90, 90)
    if (is.list(levels)) tLayers = paste(sapply(levels, paste0, collapse = '_'), collapse = '-')
        else tLayers = range(levels)
    tempFile = paste(c('temp/', paste(strsplit(mod, '/')[[1]], collapse = '_'),
                     varName, range(years), tLayers, ..., extent, modScale, '.nc'),
                     collapse = '-')
    cat("\nopening:", mod)
    cat("\n\tinto:", tempFile, "\n")
    
    if (file.exists(tempFile)) dat = brick(tempFile)
    else {
        files = list.files(paste0(dir, '/', mod, '/'), full.names = TRUE)
        if (!is.null(stream)) files = files[grepl(stream, files)]
        files = files[!grepl('rcp2p6', files)]
        files = files[!grepl('moourl', files)]
        
        

        fyr = substr(files, nchar(files)-6, nchar(files)-3)
        if (length(years) == 1) {
            if (class(years) == "list") years = years[[1]]
            years = c(years, years) 
        }
        files = files[apply( sapply(fyr, '==', years), 2, any)]

        if (sum(grepl('2006', files)) == 2) 
            files = files[!(grepl('2006', files) & !grepl('rcp', files))]
        
        varNames = strsplit(varName, ';')[[1]]
        openVar <- function(varName) {
            openLvs <- function(level) {
                dat = layer.apply(files, process.jules.file, level, varName)
            
                dat = datConvert(dat, modScale, ...)
                dat = convert_pacific_centric_2_regular(dat)
                if (any(extent != c(-180, 180, -90, 90))) dat = crop(dat, extent)
                return(dat)
            }
            if (is.list(levels)) dat = layer.apply(levels, openLvs) else dat = openLvs(levels)
              
            return(dat)
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
