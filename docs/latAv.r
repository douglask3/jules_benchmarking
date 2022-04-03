source("cfg.r")
modFracVar = 'frac'
graphics.off()
dirs = c("/hpc/data/d01/hadcam/jules_output/ALL_u-bk886_isimip_0p5deg_origsoil_dailytrif",
             "/hpc/data/d05/cburton/jules_output/u-cf137")

extent = c(-20, 55, -35, 33)
extent = c(19.75, 20.75, -33, 33)
extent = c(-180, 180, -60, 90)

xat = seq(-60, 90, by = 15)
yat = seq(0, 100, by = 20)

ModLevels = list(list(1:5, 12:13, c(6:11), 16))
obsLayers = list(list(1:2, 5, c(3:4), 8))
years_in = list(2000:2005)

burntAreaFile = "/data/users/dkelley/fireMIPbenchmarking/data/benchmarkData/MODIS250_q_BA_regridded0.5.nc" 
#dir_dat = "../fireMIPbenchmarking/data/benchmarkData/"
#obss = list(CCI = "vegfrac_refLC_refCW.nc")
#obss = lapply(obss, function(i) paste0(dir_dat, i))
obs = list.files("data/0.25deg_PFT_product/", full.names = TRUE)

type = c("TREE", "SHRUB", "GRASS", "BARE")

openType <- function(type) 
    sum(stack(obs[grepl(type, obs)]))


obs = layer.apply(type, openType) 
obs[[4]] = obs[[4]] + 1 - sum(obs)
obss = list("CCI observations" = obs)

plotLegT = c(F, T, T)
axisS = c(2)

modFacVar = 'frac'
cols = list(c(Bare = '#a6611a', Grass = '#dfc27d', Shrub = '#80cdc1', Tree = '#018571'))


sumLat <- function(lat, rs, lats) {
    if (is.null(dim(rs))) out = mean(rs[lats == lat], na.rm = TRUE)
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
    y = sumLats(x, r, lats) * scale * 5
    polygon0(x, y*100, ...)
}
    

doPlot.latAv <- function(x, lats, covs, BAs, BAscale, name = '', add = FALSE,...)  {
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
                       plotLegT = TRUE, ...) {
    
    if (plotLegT) {
        plot.new()
        legend("center", col = cols, legend = names(cols), pch = 19, pt.cex = 3, bty = 'n',                     horiz = TRUE)
    }

    years <<- years
    c(mods, mod):= openModsFromDir(dirs, years, ModLevels, extent)

    lats = yFromCell(mod[[1]][[1]], 1:length(mod[[1]][[1]][[1]]))
    x = sort(unique(lats))
    
    obs = raster::resample(obs, mod[[1]][[1]])
    modBA = lapply(mods, openMod, dirs[2],
                   'burnt_area_gb', 2001:2005, 1, extent = extent,  TRUE)
    obsBA = openObs(burntAreaFile, 1:12, 1,
                    modEG = modBA[[1]], TRUE); cat("\n")
    
    sec2yr = 60*60*24*365; colsT = make.transparent(cols, 0.8)
    doPlot.latAv(x, lats, mod[[1]], NULL, sec2yr, cols = cols, ...)        
    doPlot.latAv(x, lats, mod[[1]], NULL, sec2yr, cols = colsT, add = TRUE, ...)
    mtext("Model without fire", side = 3, line = -2, adj = 0.1)
    axis(axis, at = yat)
        
    axis(3, at = xat)  
    axis(4, at = seq(0, 100, 20), labels =seq(0, 100, 20)/5)
    doPlot.latAv(x, lats, mod[[2]], modBA, sec2yr, cols = cols, ...)
    doPlot.latAv(x, lats, mod[[2]], modBA, sec2yr, cols = colsT,  add = TRUE,  ...)
    mtext("Model with fire", side = 3, line = -2, adj = 0.1)
    axis(axis)
    axis(4, at = seq(0, 100, 20), labels =seq(0, 100, 20)/5)
        
    doPlot.latAv(x, lats, list(obs), list(obsBA), 12, cols = cols, ...) 
    mtext(obsName, side = 3, line = -2, adj = 0.1)
    axis(1, at = xat)
    axis(axis, at = yat)
    axis(4, at = seq(0, 100, 20), labels =seq(0, 100, 20)/5)
    mtext(side = 1, 'Latitude', line = 2)
}
png("figs/ISIMIP_latAV.png", height = 9, width = 5, units = 'in', res = 300)#
    layout((1:4), heights = c(0.3,1, 1, 1))
    par( mar = rep(0.75, 4), oma = rep(3.5, 4))
    
    mapply(plotModObs, obss, names(obss),  obsLayers, ModLevels, years_in, cols = cols,
          axis = axisS)
    mtext(outer = TRUE, side = 2, 'Cover (%)', line  = 2)
    mtext(outer = TRUE, side = 4, 'Burnt area (%)', line  = 2)
    #mtext(outer = TRUE, side = 1, 'Latitude', line = -3)
dev.off()


axisS
