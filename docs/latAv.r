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
    polygon(y, x, border = NA, ...)
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
            lines(yt, x, lwd = 2)
            lines(yt, x, col = col, lwd = 2, lty = 2)
        }
     }
    mapply(addCov, nlayers(r):1, cols)
}

addFireLine <- function(r, x, lats, scale, ...) {
    y = sumLats(x, r, lats) * scale * 5
    polygon0(x, y*100,...)
}
    
doPlot.latAv <- function(x, lats, covs, BAs, BAscale, name = '', add = FALSE,...)  {
    if (!add) {
        #plot(range(x), c(0, 100), type = 'n', xlab = '', ylab = '',
        #     xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs = 'i')
        
        plot(c(0, 100), range(x), type = 'n', xlab = '', ylab = '',
             xaxt = 'n', yaxt = 'n', xaxs = 'i', yaxs = 'i')
        #lblank = rep('', length(xat))
        #lapply(c(1,3), axis, at = xat, labels = lblank)
        #lblank = rep('', length(yat))
        #lapply(c(2,4), axis, at = yat, labels = lblank)
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

plotModObs <- function(obs, obsName, obsLayers, ModLevels, years, cols, ...,
                       plotFUN = plot.latAv) {
    years <<- years
    FUN <- function(...) {
        c(mods, mod):= openModsFromDir(..., stream = 'ilamb')
        mod = lapply(mod, function(i) i[[1]])
            
        return(list(mods, mod))
    }
    c(mods, modNF):= FUN(dirs[1], years, ModLevels, extent)
    c(mods, modF):= FUN(dirs[2], years, ModLevels, extent) 
    mod = list(modNF, modF)
    
    obs = raster::resample(obs, mod[[1]][[1]])
    
    modBA = lapply(mods[[1]], openMod, dirs[2],'burnt_area_gb', 2001:2020, 1, extent = extent, 
                   stream = 'ilamb', TRUE)
   
    obsBA = openObs(burntAreaFile, 1:12, 1, modEG = modBA[[1]], TRUE); cat("\n")
     
    plotFUN(mod, modBA, mods, obs, obsBA, obsName, cols, ...)

    forMod <- function(rs) {
        rs = rs / sum(rs)
        MMing <- function(i) {
            y = rs[[i]]
            x = obs[[i]]
            itmeize <- function(r) addLayer(r, 1-r)
            rarea = area(y)
            Yarea = sum((y * rarea)[], na.rm = TRUE)
            c(Yarea, score(MM(itmeize(y), itmeize(x), w = rarea)))
        }
        out = sapply(1:4, MMing)
    }
    out = lapply(mod, lapply, forMod)
    
}

plot.latAv <- function(mod, modBA, mods, obs, obsBA, obsName, cols, axis, ...) {
    plot.new()
    legend("center", col = rev(cols), legend = rev(names(cols)), pch = 19, pt.cex = 3, 
           bty = 'n', horiz = TRUE)
    
    lats = yFromCell(mod[[1]][[1]], 1:length(mod[[1]][[1]][[1]]))
    x = sort(unique(lats))
    sec2yr = 60*60*24*365; colsT = make.transparent(cols, 0.8)
    
    doPlot.latAv(x, lats, mod[[1]], NULL, sec2yr, cols = cols, ...)        
    doPlot.latAv(x, lats, mod[[1]], NULL, sec2yr, cols = colsT, add = TRUE, ...)
    mtext("Model without fire", side = 3, line = -2, adj = 0.9)
    
    axis(3, at = yat)
    axis(2, at = xat)  
    axis(4, at = xat, labels = rep('', length(xat))) 
    mtext(side = 2, 'Latitude', line = 2) 
    #axis(3, at = seq(0, 100, 20), labels =seq(0, 100, 20)/5)
    
    doPlot.latAv(x, lats, mod[[2]], modBA, sec2yr, cols = cols, ...)
    doPlot.latAv(x, lats, mod[[2]], modBA, sec2yr, cols = colsT,  add = TRUE,  ...)
    mtext("Model with fire", side = 3, line = -2, adj = 0.9)
    
    axis(3, at = yat)
    axis(1, at = seq(0, 100, 20), labels =seq(0, 100, 20)/5)
    axis(2, at = xat, labels = rep('', length(xat))) 
    axis(4, at = xat, labels = rep('', length(xat))) 
    mtext(side = 3, 'Vegetation cover (%)', line  = 2.5)
    mtext(side = 1, 'Burnt area (%)', line  = 2.3)
        
    doPlot.latAv(x, lats, list(obs), list(obsBA), 12, cols = cols, ...) 
    mtext(obsName, side = 3, line = -2, adj = 0.9)
    axis(3, at = yat)
    axis(1, at = seq(0, 100, 20), labels =seq(0, 100, 20)/5)
    axis(2, at = xat, labels = rep('', length(xat))) 
    axis(4, at = xat) 
    mtext(side = 4, 'Latitude', line = 2)
}
if (T) {
png("figs/ISIMIP_latAV.png", height = 5, width = 9, units = 'in', res = 300)#
    layout(rbind(2:4, 1), heights = c(1, 0.2))
    par( mar = c(3, 0.75, 1.5, 0.75), oma = c(0.2, 3.5, 3.5, 3.5))
    
    outs = mapply(plotModObs, obss, names(obss),  obsLayers, ModLevels, years_in, cols = cols,
          axis = axisS, MoreArgs = list(plotFUN = plot.latAv))
   
    #mtext(outer = TRUE, side = 1, 'Latitude', line = -3)
dev.off()
}
browser()
regionMask = raster('data/biomAssigned.nc')
regionMask[is.na(regionMask)] = 0
plot.vegMix <- function(mod, modBA, mods, obs, obsBA, obsName, cols, axis, ...) {
    browser()
    regionMask = raster::resample(regionMask, mod[[1]][[1]])
    xbin = seq(0, 100, )/100
    xb = xbin[-1] - diff(xbin)
        dev.new()
    forRegion <- function (id) {
        mask = regionMask == id
        
        fracQuant <- function(r) {
            tot = sum(r)
           
            r = r/tot
            vals = r[[1:3]][mask]
            vals = vals[sample(1:nrow(vals), size = min(c(1000, nrow(vals))), replace = FALSE),]
             
            #tot = apply(vals, 1, sum)

            for (i in 2:ncol(vals)) 
                vals[, i] = (vals[, i] + vals[, i-1])
            
            binMakeupPlot <- function(i) {
                b1 = xbin[i]; b2 = xbin[i-1]
                if (b2 == 0) b2 = -1
                if (b1 == 1) b1 = 2
                index = which(vals[,3] >= b2 & vals[,3] <= b1)
                  
                addLine <- function(j) {
                    forLine <- function(v, col) {
                        x = id + c(-1, 1) * wds[i-1] * v/vals[j, 3]
                        if (vals[j, 3] == 0) vals[j,3] =  vals[j, 3] + runif(1, -0.01, 0.01)
                        #if (id == 8) browser()
                        lines(x, rep(vals[j, 3], 2)^4, col = col, lwd = 0.1)
                    }
                    mapply(forLine, rev(vals[j,]), cols[-1])
                }
                sapply(index, addLine)
            }

            wds = hist(vals[,3], breaks = xbin, plot = FALSE)$density
            wds = 0.5*wds / max(wds)            

            vals = sapply(2:length(xbin), binMakeupPlot)
           
            return()

            vals = apply(vals, 2, function(x)
                            hist(x, breaks = xbin/100, plot = FALSE)[['counts']])

            vals = matrix2list(vals)
            
            mv = max(unlist(vals))
            addPoly <- function(vs, col) {
                browser()
            }
            mapply(addPoly, rev(vals), cols[-1])
            browser()
        }
        extractFrac <- function(rs) 
            layer.apply(rs, fracQuant)
        
        #vmod = lapply(mod, extractFrac)
        extractFrac(mod[[1]])
        #browser()
    }
    ids = unique(regionMask)[-1]
    plot(range(ids) + c(-0.5, 0.5), range(xbin))
    lapply(ids, forRegion)
    browser()
}


png("figs/ISIMIP_vegFracs.png", height = 9, width = 5, units = 'in', res = 300)#
    layout((1:4), heights = c(0.3,1, 1, 1))
    par( mar = rep(0.75, 4), oma = rep(3.5, 4))
    
    mapply(plotModObs, obss, names(obss),  obsLayers, ModLevels, years_in, cols = cols,
          axis = axisS, MoreArgs = list(plotFUN = plot.vegMix))
    mtext(outer = TRUE, side = 2, 'Cover (%)', line  = 2)
    mtext(outer = TRUE, side = 4, 'Burnt area (%)', line  = 2)
    #mtext(outer = TRUE, side = 1, 'Latitude', line = -3)
dev.off()
