openModsFromDir <- function(dirs, years, ModLevels,
                            extent = NULL, varname = 'frac', scale = 1) {
    openMods <- function(dir, mods) {
        
        mod = lapply(mods, openMod, dir, varname, years, scale, extent = extent,
                     levels = ModLevels, cover = TRUE)
       return(mod)
    }
    mods = paste0('/', c("GFDL-ESM2M", "HADGEM2-ES", "IPSL-CM5A-LR", "MIROC5"), '/')   
    #mods = lapply(dirs, function(dir) models)
    mod = mapply(openMods, dirs, mods, SIMPLIFY = FALSE)
    return(list(mods, mod))
}
