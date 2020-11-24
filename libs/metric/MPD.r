MPD_only <- function(x, y, w = NULL) 
    out = MPDonly(x, y)[[1]]
    
null.MPD_only <- function(x, n = 5) 
    null.FUN(x, FUN = MPDonly, n = 5)

