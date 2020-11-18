library(raster)
library(gitBasedProjects)
library(rasterExtras)
sourceAllLibs("libs/")
sourceAllLibs('../benchmarkmetrics_github/benchmarkMetrics/R/')
sourceAllLibs('../rasterextrafuns/rasterPlotFunctions/R/')
graphics.off()

dir = '../jules_outputs/'

mods = c('u-by849', 'u-by851')
varName = 'burnt_area_gb'
modScale = 60*60*24*360*100

obss = paste0("../fireMIPbenchmarking/data/benchmarkData/",c(
              "GFED4s_v2.nc", "GFED4.nc", "MCD45.nc", "meris_v2.nc",
              "MODIS250_q_BA_regridded0.5.nc"))
names(obss) = c("GFED4s", "GFED4", "MCD45", "Meris", "CCI")
years = list(1997:2014, 1997:2014, 2001:2014, 2006:2013, 2001:2014)
obsLayers = list(1:12, 1:12, 1:12, 1:12, 1:12)
obsScale = 12*100

limits_aa = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)                
cols_aa = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
            '#fc4e2a','#e31a1c','#bd0026','#800026')

dlimits_aa = c(-20, -10, -5, -1, -0.5, -0.1, 0.1, 0.5, 1, 5, 10, 20)                
dcols_aa = rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

openMod <- function(mod, dir, varName, years, annualAver = FALSE, modScale) {
    tempFile = paste(c('temp/', mod, varName, range(years), annualAver, modScale, '.nc'),
                     collapse = '-')
    if (file.exists(tempFile)) return(brick(tempFile))
    files = list.files(paste0(dir, mod, '/'), full.names = TRUE)
    fyr = substr(files, nchar(files)-6, nchar(files)-3)
   
    files = files[apply( sapply(fyr, '==', years), 2, any)]
    dat = layer.apply(files, process.jules.file, 1, varName)
    if (annualAver) dat = mean(dat)
    dat = convert_pacific_centric_2_regular(dat)
    dat = dat * modScale
    writeRaster(dat, file = tempFile, overwrite = TRUE)
    return(dat)
}

openObs <- function(obs, Layers, scale, annualAver = FALSE, modEG) { 
    
    tempFile = paste(c('temp/', filename.noPath(obs, TRUE),
                     range(Layers), annualAver, scale, '.nc'),
                     collapse = '-')
    if (file.exists(tempFile)) return(brick(tempFile))
    obs = brick(obs)
    obs = obs[[Layers]]
    if (annualAver) obs = mean(obs)
    obs = obs * scale
    modEG = modEG[[1]]
    obs = raster::resample(obs, modEG)
    obs[is.na(modEG)] = NaN
    return(obs)
}

annual_average_NME <- function(mods, obss, fname) {   
    png(fname, height = 4 * 15 *(length(obss)+2)/36, width = 4*(length(mods)+1),
        units = 'in', res = 300)
    par(mfcol = c(length(obss)+2, length(mods)+1), mar = rep(0, 4))
    plot.new()

    openYears <- function(year) 
        lapply(mods, openMod, dir, varName, year, TRUE, modScale)

    #mods = openMod(dir, mod, varName, year, TRUE) * modScale
    mods_names = mods
    mods = lapply(years, openYears)
    
    obs_names = names(obss) #sapply(obss, filename.noPath, TRUE)
    obss = mapply(openObs, obss, obsLayers, obsScale,
                  MoreArgs = list(annualAver = TRUE, modEG = mods[[1]][[1]]))
    mapply(plotStandardMap, obss, obs_names,
           MoreArgs = list(cols = cols_aa, limits = limits_aa))

    StandardLegend(cols_aa, limits_aa, dat = obss[[1]], oneSideLabels = FALSE)    
    
    plotMod <- function(mod, name) {
        plotStandardMap(mod[[1]], '', cols = cols_aa, limits = limits_aa)
        mtext(name, side = 3, adj = 0.1)
        plotDiff <- function(obs, md) 
            plotStandardMap(md[[1]] - obs, '', cols = dcols_aa, limits = dlimits_aa)
        mapply(plotDiff, obss, mod)
        
        StandardLegend(dcols_aa, dlimits_aa, dat = mod[[1]][[1]]-obss[[1]], oneSideLabels = FALSE)  
    }
    mods_switch = lapply(1:2, function(i) lapply(mods, function(m) m[[i]]))
    mapply(plotMod, mods_switch, mods_names)
    dev.off()
    NMEout <- function(mod, obs) {
        scores = NME( obs, mod[[1]])
        null = null.NME(obs, n = 5)
        null = c(summary(null)[1:2], summary(null)[3] + c(-1, 1) * summary(null)[4])
        out = t(sapply(score(scores), function(i) c(null, i)))
    }
    score = mapply(function(obs, mod) lapply(mod, NMEout, obs), obss, mods)
    scoreI = score[1,]
    if (nrow(score) == 2)
        scoreO = mapply(function(i, j) cbind(i, j[,ncol(j)]), scoreI, score[2,], SIMPLIFY = FALSE)
    scoreO = do.call(rbind, scoreO)
    rownames(scoreO) = paste(rep(obs_names, each = 3), ' - step', 1:3)
    colnames(scoreO) = c("Median null", "Mean null", "RR null lower", "RR null upper",
                         mods_names)
    return(scoreO)
}

score_aa = annual_average_NME(mods, obss, 'figs/burnt_area_aa.png')




