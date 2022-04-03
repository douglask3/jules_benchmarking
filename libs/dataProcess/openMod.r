openMod <- function(mod, dir, varName, years, modScale, ..., fill = NULL,
                    levels = 1, layer = NULL, fileID = 'ilamb', fileID2 = '',
                    temp_dir_move = NULL) {
    
    if (!exists("extent") || is.null(extent) || class(extent) == "standardGeneric")
        extent = c(-180, 180, -90, 90)
    if (is.list(levels)) tLayers = paste(sapply(levels, paste0, collapse = '_'), collapse = '-')
        else tLayers = range(levels)
    
    tempFile = paste(c('temp/', filename.noPath(dir), paste(strsplit(mod, '/')[[1]], collapse = '_'),
                     varName, range(years), tLayers, ..., fileID, fileID2, extent, modScale, '.nc'),
                     collapse = '-')
    cat("\nopening:", mod)
    cat("\n\tinto:", tempFile, "\n")
    
    if (file.exists(tempFile)) dat = brick(tempFile)
    else {
        dir = paste0(dir, mod, '/')
        listJfiles <- function(dir, years) {
            files = list.files(dir, full.names = TRUE)
            if (fileID != '') files = files[grepl(fileID, files)]
            fyr = substr(files, nchar(files)-6, nchar(files)-3)
            if(length(years == 1)) years = c(years, years) 
            files = files[apply( sapply(fyr, '==', years), 2, any)] 
        }
        if (!is.null(temp_dir_move)) {
            temp_dir_move = paste0(temp_dir_move, '/', mod, '/')
            makeDir(temp_dir_move)
            if (length(list.files(temp_dir_move)) == 0) {
                files = listJfiles(dir, 1950:9999)
                 
                file.copy(files, temp_dir_move)
            }
            
            dir = temp_dir_move
        }
        files = listJfiles(dir, years)
        
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
        if (!is.null(fill)) dat[!is.na(dat)] = fill
        dat = writeRaster(dat, file = tempFile, overwrite = TRUE)
    }
    if (!is.null(layer)) dat = dat[[layer]]
    if (nlayers(dat) == 1) dat = dat[[1]]
    
    return(dat)
}
