openModsFromDir <- function(dirs, years, ModLevels,
                            extent = NULL, varname = 'frac', scale = 1, ...) {
    openMods <- function(dir, mods) {
    
        mod = lapply(mods, openMod, dir, varname, years, scale, extent = extent,
                     levels = ModLevels, cover = TRUE, exName = tail(strsplit(dir, '/')[[1]], 1), ...)
       return(mod)
    }
    
    mods = lapply(dirs, function(dir) list.dirs(dir, full.names = FALSE)[-1])
    mod = mapply(openMods, dirs, mods[[1]], SIMPLIFY = FALSE)
    return(list(mods, mod))
}