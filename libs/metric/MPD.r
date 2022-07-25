MPD_only <- function(x, y, w = NULL) {
    if (is.raster(x)) { 
        mask = !is.na(x + y)
        x = matrix(x[mask]); y = matrix(y[mask]); 
        if (!is.null(w)) w = matrix(w[mask])
    }
    
    out = MPDonly(x, y, w)
}
null.MPD_only <- function(x, n = 5) 
    null.FUN(x, FUN = MPDonly, n = 5)

