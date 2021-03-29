source("cfg.r")

dirs = paste0("../jules_outputs/", c("ALL_u-bk886_isimip_0p5deg_origsoil_dailytrif",
                                     "RUNC20C_u-by276_isimip_0p5deg_origsoil_dailytrif_fire/"))
extent = NULL
ModLevels = list(Tree = 1:5, Wood = c(1:5, 12:13), Shrub = 12:13, Grass = c(6:11),
                 "Total Vegetation" = 16)

dir_dat = "../fireMIPbenchmarking/data/benchmarkData/"
obss = list(CCI = "vegfrac_refLC_refCW.nc",
            VCF = c("treecover2000-2014.nc", "nontree2000-2014.nc"))
obss = lapply(obss, function(i) paste0(dir_dat, i))
obsLayers = list(list(1:2, c(1:2, 5), 5, c(3:4), 8), 1:60)

cols = c('#a6611a','#dfc27d','#80cdc1','#018571')

limits_cover = seq(0, 0.9, 0.1)*100
    cols_cover = c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')

c(mods, mod):= openModsFromDir(dirs, 2001:2006, ModLevels, extent)
mod = lapply(mod, function(i) lapply(i, function(k) 100*k/sum(k[[c(1, 3, 4, 5)]])))

obs = mapply(openObs,obss, obsLayers, 1, MoreArgs = list(modEG = mod[[1]][[1]], cover = TRUE)); cat("\n")
obs[[1]] = 100*obs[[1]]/sum(obs[[1]][[c(1, 3, 4, 5)]])
obs[[2]] = obs[[2]]/length(obsLayers[[2]])
obs[[2]] = addLayer(obs[[2]], 100 - sum(obs[[2]]))
blank = obs[[2]][[1]]
blank[] = NaN
obs[[2]] = addLayer(blank, obs[[2]][[1]], blank, obs[[2]][[2:3]])


summeriseModSet <- function(mod) {
    op <- function(FUN, i) 
        FUN(layer.apply(mod, function(md) md[[i]]))

    lapply(1:nlayers(mod[[1]]), function(i) addLayer(op(mean, i), op(sd.raster, i)))
}
mod = lapply(mod, summeriseModSet)
list2layer <- function(r) layer.apply(r, function(i) i)

png("figs/coverMaps.png", height = 9, width = 12, res = 300, units = 'in')
    par(mfcol = c(5, 4), mar = rep(0, 4), oma = rep(2, 4))
    layout(rbind(cbind(1:5, 6:10, 11:15, 16:20), 21), heights = c(rep(1, 5), 0.4))
    
    
    lapply(mod, function(i) lapply(i, plotStandardMap, txt = '',
                                   limits = limits_cover, cols = cols_cover))
    
    plotObs <- function(i, ...) {
        if (sum.raster(!is.na(i)) == 0) plot.new()
        else plotStandardMap(i, ...)
    }
    lapply(obs, function(i) layer.apply(i, plotObs,
                                        txt = '', limits = limits_cover, cols = cols_cover))
    StandardLegend(dat = mod[[1]], cols = cols_cover, limits = limits_cover, add = FALSE)
    

    mapply(mtext, c("Tree", "Wood", "Shrub", "Grass", "Bare"),
            adj = 1-(seq(0.1, 1, 0.2)*6)/6.4, side = 2, outer = TRUE)

    mapply(mtext, c("Model no fire", "Model with fire", "CCI", "VCF"),
           adj = seq(0.125, 1, 0.25), side = 3, outer = TRUE)
dev.off()


c(mods, modBA):= openModsFromDir(dirs[2], 2001:2005, 1,
                                 varname = 'burnt_area_gb')

modBA = list2layer(modBA[[1]])
modBA = addLayer(mean(modBA)*60*60*24*365, sd.raster(modBA))

obs = paste0(dir_dat,c("GFED4s_v2.nc", "../GFED4.fBA.r0d5.1995.2013.nc",
                        "MCD45.nc", "meris_v2.nc",
                        "MODIS250_q_BA_regridded0.5.nc"))
names(obs) = c("GFED4s", "GFED4", "MCD45", "CCI Meris", "CCI MODIS")
obsLayers = list(1:60, 8:67, 1:60, 1:60, 1:60)
obsBA = mapply(openObs,obs, obsLayers, MoreArgs = list(100*12, modEG = modBA[[1]], TRUE))

limits_BA = c(0, 0.1, 1, 2, 5, 10, 20, 50)
cols_BA = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
png("figs/fireMaps2.png", height = 9, width = 12, res = 300, units = 'in')
    par(mfcol = c(5, 4), mar = rep(0, 4), oma = rep(2, 4))
    layout(rbind(matrix(1:6, ncol = 2), 7), heights = c(rep(1, 3), 0.4))
    plotStandardMapTitle <- function(r, txt, ...) {
        plotStandardMap(r, ...)
        mtext(side = 3, txt, line = -1)
    }
    mapply(plotStandardMapTitle, c(modBA, obsBA), c("Model", names(obs)),
           MoreArgs = list(limits = limits_BA, cols = cols_BA))

    StandardLegend(dat = modBA, cols = cols_BA, limits = limits_BA, add = FALSE)
dev.off()
browser()
    mapply(mtext, c("Tree", "Shrub", "Grass", "Bare"),
           adj = seq(0.125, 1, 0.25), side = 3, outer = TRUE)

mapply(mtext, c("HadGEM-ES-2", "CCI"), adj = seq(0.5, 1-0.33/2), side = 2, outer = TRUE)

