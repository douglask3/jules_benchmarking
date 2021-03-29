openModsFromDir <- function(dirs, years, ModLevels,
                            extent = NULL, varname = 'frac', scale = 1) {
    openMods <- function(dir, mods) {
    
        mod = lapply(mods, openMod, dir, varname, years, scale, extent = extent,
                     levels = ModLevels, cover = TRUE)
       return(mod)
    }
    
    mods = lapply(dirs, function(dir) list.dirs(dir, full.names = FALSE)[-1])
    mod = mapply(openMods, dirs, mods, SIMPLIFY = FALSE)
    return(list(mods, mod))
}
