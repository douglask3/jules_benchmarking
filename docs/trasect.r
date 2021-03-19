source("cfg.r")
graphics.off()
dir = "../jules_outputs/RUNC20C_u-by276_isimip_0p5deg_origsoil_dailytrif_fire/hadgem2-es/"

extent = c(-20, 55, -35, 33)
extent = c(19.75, 20.75, -35, 33)

ModLevels = list(list(1:5, 12:13, c(6:11), 16), list(c(1:5, 12:13), c(6:11), 16))
obsLayers = list(list(1:2, 5, c(3:4), 8), 1:60)
years_in = list(2005:2006, 2001:2006)

dir_dat = "../fireMIPbenchmarking/data/benchmarkData/"
obss = list(CCI = "vegfrac_refLC_refCW.nc",
            VCF = c("treecover2000-2014.nc", "nontree2000-2014.nc"))
obss = lapply(obss, function(i) paste0(dir_dat, i))

cols = list(c(Bare = '#a6611a', Grass = '#dfc27d', Shrub = '#80cdc1', Tree = '#018571'),
            c(Bare = '#d8b365', Grass = '#ffffbf', Wood = '#5ab4ac'))



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
addOverlay <- function(r, cols, poly = TRUE) {
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
    mapply(addCov, nlayers(r):1, cols)
}

addFireLine <- function(r, scale) {
    y = sumLats(r) * scale
    polygon0(x, y*100, col = "black", density = 20)
}
    

doPlot <- function(covs, BAs, BAscale, name = '', ...)  {
    plot(range(x), c(0, 100), type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs = 'i')
    mtext(side = 3, adj = 0.1, name)
    addOverlay(covs, ...)
    addFireLine(BAs, BAscale)
}

plotModObs <- function(obs, obsName, obsLayers, ModLevels, years, axis, cols, ...) {
    years <<- years
    mod = openMod('', dir, 'frac', years, 1, extent = extent,
                  levels = ModLevels, cover = TRUE)
    
    obs = openObs(obs, obsLayers, 1, modEG = mod, cover = TRUE); cat("\n")
    if (nlayers(obs) == 2) {
        obs = obs / (100*length(obsLayers))
        obs = addLayer(obs, 1 - sum(obs))
    }
    modBA = openMod('', dir, 'burnt_area_gb', 2001:2005, 1, extent = extent,  TRUE)
    obsBA = openObs("../fireMIPbenchmarking/data/benchmarkData/GFED4s_v2.nc", 1:12, 1,
                    modEG = mod, TRUE); cat("\n")

    doPlot(mod, modBA, 60*60*24*365, cols = cols, ...) 
    if (axis == 2) mtext("Model + fire", side = 3, line = -2, adj = 0.1)
    axis(axis)
    
    doPlot(obs, obsBA, 12, , cols = cols, ...) 
    mtext(paste("Observations:", obsName), side = 3, line = -2, adj = 0.1)
    axis(1)
    axis(axis)
    plot.new()
    legend("center", col = cols, legend = names(cols), pch = 19, pt.cex = 3, bty = 'n', horiz = TRUE)

}
png("figs/ISIMIP_transect.png", height = 6.5, width = 7.2, units = 'in', res = 300)#
    layout(cbind(1:3, 4:6), heights = c(1,1,  0.3))
    par( mar = rep(0.5, 4), oma = rep(3.5, 4))

    mapply(plotModObs, obss, names(obss),  obsLayers, ModLevels, years_in, cols = cols,
          axis = c(2, 4))
    mtext(outer = TRUE, side = 2, 'Cover (%)', line  = 2)
    mtext(outer = TRUE, side = 1, 'Latitude', line = -3)
dev.off()
