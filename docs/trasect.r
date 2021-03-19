source("cfg.r")

dir = "../jules_outputs/RUNC20C_u-by276_isimip_0p5deg_origsoil_dailytrif_fire/hadgem2-es/"

extent = c(-20, 55, -35, 33)
extent = c(19.75, 20.75, -35, 33)
ModLevels = list(1:5, 12:13, c(6:11), 16)

obsLayers = list(1:2, 5, c(3:4), 8)
cols = c('#a6611a','#dfc27d','#80cdc1','#018571')
obs =  "../fireMIPbenchmarking/data/benchmarkData/vegfrac_refLC_refCW.nc"

mod = openMod('', dir, 'frac', 2001:2005, 1, extent = extent,  levels = ModLevels, cover = TRUE)
obs = openObs(obs, obsLayers, 1, modEG = mod, cover = TRUE); cat("\n")

lats = yFromCell(mod, 1:length(mod[[1]]))
x = sort(unique(lats))


sumLat <- function(lat, rs) {
    if (is.null(dim(rs))) out = mean(rs[lats == lat])
        else out = apply(rs[lats == lat,], 2, mean, na.rm = TRUE)
    return(out)
}
sumLats <- function(r) 
    t(sapply(x, sumLat, r[]))

polygon0 <- function(x, y, ...) {
    test = !is.na(y)
    x = x[test]; y = y[test]
    polygon(c(x[1], x, tail(x, 1)), c(0, y, 0), border = NA, ...)
}
addOverlay <- function(r, poly = TRUE) {
    covs =  sumLats(r)
    
    for (i in 2:ncol(covs)) covs[,i] = covs[,i] + covs[, i-1]
    covs = apply(covs, 2, function(i) i/covs[, ncol(covs)])
    
    
    addCov <- function(i, col) {
        y = covs[,i] * 100       
        if (poly) 
           polygon0(x, y, col = col)
        else {
            lines(x, y, lwd = 2)
            lines(x, y, col = col, lwd = 2, lty = 2)
        }
     }
    mapply(addCov, 4:1, cols)
}

#addOverlay(mod)

modBA = openMod('', dir, 'burnt_area_gb', 2001:2005, 1, extent = extent,  TRUE)
obsBA = openObs("../fireMIPbenchmarking/data/benchmarkData/GFED4s_v2.nc", 1:12, 1, modEG = mod, TRUE); cat("\n")

addFireLine <- function(r, scale) {
    y = sumLats(r) * scale
    polygon0(x, y*100, col = "black", density = 20)
}
    

doPlot <- function(covs, BAs, BAscale, name = '')  {
    plot(range(x), c(0, 100), type = 'n', xlab = '', ylab = '', xaxt = 'n')
    mtext(side = 3, adj = 0.1, name)
    addOverlay(covs)
    addFireLine(BAs, BAscale)
}
png("figs/ISIMIP_transect.png", height = 10, width = 7, units = 'in', res = 300)
par(mfrow = c(2, 1), mar = rep(0.5, 4), oma = rep(3.5, 4))
doPlot(mod, modBA, 60*60*24*365) 
doPlot(obs, obsBA, 12) 
axis(1)
mtext(outer = TRUE, side = 2, 'Cover (%)', line  = 2, 'HadGEM2-ES')
mtext(outer = TRUE, side = 1, 'Latitude', line = 2, 'CCI/GFED4s')
dev.off()
