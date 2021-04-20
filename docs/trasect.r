source("cfg.r")
modFracVar = 'frac'
graphics.off()
dirs = paste0("../jules_outputs/", c("ALL_u-bk886_isimip_0p5deg_origsoil_dailytrif",
                                     "RUNC20C_u-by276_isimip_0p5deg_origsoil_dailytrif_fire/"))

extent = c(-20, 55, -35, 33)
extent = c(19.75, 20.75, -33, 33)

xat = seq(-40, 40, by = 10)
yat = seq(0, 100, by = 20)

ModLevels = list(list(1:5, 12:13, c(6:11), 16), list(1:5, 12:13, c(6:11), 16),
                 list(c(1:5, 12:13), c(6:11), 16))
obsLayers = list(list(1:2, 5, c(3:4), 8), list(1:2, 5, c(3:4), 8), 1:60)
years_in = list(2005:2006, 2005:2006, 2001:2006)

dir_dat = "../fireMIPbenchmarking/data/benchmarkData/"
obss = list(CCI = "vegfrac_refLC_refCW.nc",
            IGBP = "vegfrac_igbp.nc",
            VCF = c("treecover2000-2014.nc", "nontree2000-2014.nc"))
obss = lapply(obss, function(i) paste0(dir_dat, i))

plotModT = c(T, F, T)
plotLegT = c(F, T, T)
axisS = c(2, 2, 4)

cols = list(c(Bare = '#a6611a', Grass = '#dfc27d', Shrub = '#80cdc1', Tree = '#018571'),
            c(Bare = '#a6611a', Grass = '#dfc27d', Shrub = '#80cdc1', Tree = '#018571'),
            c(Bare = '#d8b365', Grass = '#ffffbf', Wood = '#5ab4ac'))


sumLat <- function(lat, rs, lats) {
    if (is.null(dim(rs))) out = mean(rs[lats == lat])
        else out = apply(rs[lats == lat,], 2, mean, na.rm = TRUE)
    return(out)
}
sumLats <- function(x, r, ...) {
    #browser()
    t(sapply(x, sumLat, r[], ...))
}

polygon0 <- function(x, y, yb = NULL, ...) {
    if (is.null(yb)) yb = rep(0, length(y))
    test = !is.na(y + yb)
    x = x[test]; y = y[test]; yb = yb[test]
    x = c(x, rev(x))
    y = c(y, rev(yb))
    polygon(x, y, border = NA, ...)
}
addOverlay <- function(r, x, lats, cols, poly = TRUE) {
    
    covs =  sumLats(x, r, lats)
    
    for (i in 2:ncol(covs)) covs[,i] = covs[,i] + covs[, i-1]
    covs = apply(covs, 2, function(i) i/covs[, ncol(covs)]) * 100
    
    addCov <- function(i, col) {
        yt = covs[,i] 
        if (i> 1) yb = covs[,i-1] else yb = rep(0, length(yt))
        if (poly) {
            polygon0(x, yt, yb, col = col)
        } else {
            lines(x, yt, lwd = 2)
            lines(x, yt, col = col, lwd = 2, lty = 2)
        }
     }
    mapply(addCov, nlayers(r):1, cols)
}

addFireLine <- function(r, x, lats, scale, ...) {
    y = sumLats(x, r, lats) * scale
    polygon0(x, y*100, ...)
}
    

doPlot <- function(x, lats, covs, BAs, BAscale, name = '', add = FALSE,...)  {
    if (!add) {
        plot(range(x), c(0, 100), type = 'n', xlab = '', ylab = '',
             xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs = 'i')
        lblank = rep('', length(xat))
        lapply(c(1,3), axis, at = xat, labels = lblank)
        lblank = rep('', length(yat))
        lapply(c(2,4), axis, at = yat, labels = lblank)
    }  
    mtext(side = 3, adj = 0.1, name)
    
    for (nn in 1:4) {
        lapply(covs, addOverlay, x, lats, ...)
        lapply(rev(covs), addOverlay, x, lats, ...)
    }
    if (length(BAs)==1) {
        col = "black"
        density = 30
    }  else {
        col = make.transparent("black", 0.5)
        density = 15
    }
    if (!is.null(BAs)) {
        mapply(addFireLine, BAs, BAscale, density = density, col = col,
               angle =  c(45, 90, 135, 180)[1:length(BAs)], MoreArgs = list(x = x, lats = lats))
    }
    
}


plotModObs <- function(obs, obsName, obsLayers, ModLevels, years, axis, cols,
                       plotModT, plotLegT = plotModT, ...) {
    
    if (plotLegT) {
        plot.new()
        legend("center", col = cols, legend = names(cols), pch = 19, pt.cex = 3, bty = 'n',                     horiz = TRUE)
    }

    years <<- years
    c(mods, mod):= openModsFromDir(dirs, years, ModLevels, extent)

    lats = yFromCell(mod[[1]][[1]], 1:length(mod[[1]][[1]][[1]]))
    x = sort(unique(lats))
    
    obs = openObs(obs, obsLayers, 1, modEG = mod[[1]][[1]], cover = TRUE); cat("\n")
    if (nlayers(obs) == 2) {
        obs = obs / (100*length(obsLayers))
        obs = addLayer(obs, 1 - sum(obs))
    }
    modBA = lapply(mods[[2]], openMod, dirs[2],
                   'burnt_area_gb', 2001:2005, 1, extent = extent,  TRUE)
    obsBA = openObs("../fireMIPbenchmarking/data/benchmarkData/GFED4s_v2.nc", 1:12, 1,
                    modEG = modBA[[1]], TRUE); cat("\n")
   
    if (plotModT) {
        doPlot(x, lats, mod[[1]], NULL, 60*60*24*365, cols = cols, ...)        
        doPlot(x, lats, mod[[1]], NULL, 60*60*24*365, cols = make.transparent(cols, 0.8), 
               add = TRUE, ...)
        mtext("Model without fire", side = 3, line = -2, adj = 0.1)
        axis(axis, at = yat)
        
        axis(3, at = xat)        
        doPlot(x, lats, mod[[2]], modBA, 60*60*24*365, cols = cols, ...)
        doPlot(x, lats, mod[[2]], modBA, 60*60*24*365, cols = make.transparent(cols, 0.8), 
               add = TRUE,  ...)
        mtext("Model with fire", side = 3, line = -2, adj = 0.1)
        axis(axis)
        
    }
    
    doPlot(x, lats, list(obs), list(obsBA), 12, , cols = cols, ...) 
    mtext(paste("Observations:", obsName), side = 3, line = -2, adj = 0.1)
    axis(1, at = xat)
    axis(axis, at = yat)
    
    mtext(side = 1, 'Latitude', line = 2)
}
png("figs/ISIMIP_transect.png", height = 9, width = 7.2, units = 'in', res = 300)#
    layout(cbind(1:5, 6:10), heights = c(0.3,1, 1, 1, 1))
    par( mar = rep(0.25, 4), oma = rep(3.5, 4))

    mapply(plotModObs, obss, names(obss),  obsLayers, ModLevels, years_in, cols = cols,
          axis = axisS, plotModT = plotModT)
    mtext(outer = TRUE, side = 2, 'Cover (%)', line  = 2)
    #mtext(outer = TRUE, side = 1, 'Latitude', line = -3)
dev.off()


axisS
