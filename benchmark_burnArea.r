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

obss = paste0("../fireMIPbenchmarking/data/benchmarkData/",
              c("GFED4s_v2.nc", "GFED4.nc", "MCD45.nc", "meris_v2.nc",
                "MODIS250_q_BA_regridded0.5.nc"))
names(obss) = c("GFED4s", "GFED4", "MCD45", "Meris", "CCI")
years = list(1997:2013, 1997:2013, 2001:2013, 2006:2012, 2001:2013)
obsLayers = list(48:204, 48:204, 1:156, 1:84, 1:156)
obsScale = 12*100

limits_aa = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)                
cols_aa = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
            '#fc4e2a','#e31a1c','#bd0026','#800026')

dlimits_aa = c(-20, -10, -5, -1, -0.5, -0.1, 0.1, 0.5, 1, 5, 10, 20)                  
dcols_aa = rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

limits_tr = c(-2, -1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1, 2)/10         
cols_tr = rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

dlimits_tr = limits_tr            
dcols_tr = cols_tr

cols_phase = c('#313695', '#a50026', '#ffff00','#313695')
limits_phase = 0.5:11.5 

dcols_phase = c("#f7f7f7", "#b35806", "#003300", "#542788", "#f7f7f7")
dlimits_phase = (-5.5:5.5)


cols_conc = c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4',
                                       '#1d91c0','#225ea8','#253494','#081d58')
limits_conc =  seq(0, 0.9, 0.1)


dcols_conc = c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30')
dlimits_conc =  c(-0.8, -0.6, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.6, 0.8)

cols_modal = c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1','#dd3497','#ae017e','#7a0177','#49006a')#c('#34eeba','#d95f02','#7570b3')
limits_modal = c(1, 1.1, 1.2, 1.5, 2)


dcols_modal = c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')
dlimits_modal = c(-1, -0.5, -0.2, -0.1, 0.1, 0.2, 0.5, 1)


logit <- function(x) {
    x[x < 0.00000000001] = 0.00000000001
    x[x > 0.99999999999] = 0.99999999999
    log(x/(1-x))
}

logistic <- function(x) 
    1/(1+exp(x*(-1)))

TrendFun <- function(x) {
    
    fit = lm(x ~ t, data = data.frame(x = x, t = 1:length(x)))
    
    #res = try(summary(fit)[[4]][2, 3:4])
    res = try(coefficients(fit)[2], summary(fit)[[4]][2, 4])
    if (is.na(res[1])) return(c(0, 1))  
    if (class(res) == "try-error") return(c(-999, 0.0))
    
    #if (abs(res[1]) > 8) browser()
    return(res)
}

#x = yay
#xDT = x
#for (i in 1:length(x)) {
#    yr = ceiling(i/12)-1
#    mn = i - yr*12
#    xDT[i] = x[i] + wow[mn] * yr
#}

makeTrendCoe <- function(dat) {
    
    dat = dat/1200
    dat = logit(dat)
    yay <<- dat
    datO = dat[[1:2]]
    mask = !is.na(dat[[1]])
    climtol.pnt <- function(m, x) {
        m = seq(m, nlayers(dat), by = 12)
        mdat = x[m]
        return(TrendFun(mdat))
    }

    climtol <- function(m) {
        m = seq(m, nlayers(dat), by = 12)
        mdat = dat[[m]][mask]
        mtr = apply(mdat, 1, TrendFun)
        
        datO[mask] = t(mtr)
        return(datO)
    }
    if (nlayers(dat) > 24) tr = lapply(1:12, climtol)
    else browser()
    #browser()
    
    #trend = layer.apply(tr, function(i) i[[1]])
    #browser()
    #mdat = dat[mask]
    #mtr = apply(mdat, 1, TrendFun)
    #datO[mask] = t(mtr)
    #return(datO)
    #tr = datO
    #browser()
    #trend = mean(layer.apply(tr, function(i) i[[1]]*(1-i[[2]])))
    #pTot = sum(layer.apply(tr, function(i) 1-i[[2]]))
    #trend = trend / pTot
    pval = sum(layer.apply(tr, function(i) -2 * log(i[[2]])))
    pval[mask] = 1-pchisq(pval[mask], 24)
    
    trend = layer.apply(tr, function(i) i[[1]])#*(1-i[[2]])))
    
    datDT = dat
    browser()
    for (i in 1:nlayers(dat)) {
        yr = ceiling(i/12)-1
        mn = i - yr*12
        print(yr)
        print(mn)
        print('---')
        datDT[[i]] = dat[[i]] - trend[[mn]] * yr
    }
    ldat = logistic(dat)
    ldatDT = logistic(datDT)
    
    mdat = mean(logistic(dat))
    mdatDT =  mean(logistic(datDT))
    diff = mdat - mdatDT
    testn = diff < 0
    testp = !testn  
    
    diff[testn] = diff[testn] / mdatDT[testn]
    diff[testp] = diff[testp] / (1-mdatDT[testp])
    
    return(addLayer(diff, pval))
}

PhaseConcMod <- function(dat, modEG = NULL) {
    clim = layer.apply(1:12, function(mn) mean(dat[[seq(mn, nlayers(dat), by = 12)]]))
    if (!is.null(modEG)) clim = raster::resample(clim, modEG)
    PC = PolarConcentrationAndPhase(clim, 'months')
    Modal = Modalise(clim)
    return(addLayer(PC, Modal))
}


datConvert <- function(dat, scale, annualAver = FALSE, trend = FALSE, Seasonal = FALSE, 
                       modality = FALSE,
                       modEG = NULL) {
    RS <- function(dat) {
        if (!is.null(modEG)) dat = raster::resample(dat, modEG)
        return(dat)
    }
    if (annualAver) {
        dat = mean(dat) * scale
        dat = RS(dat)
    } else if (trend) {        
        dat = RS(dat)
        dat = makeTrendCoe(dat * scale)
    } else if (Seasonal) 
        dat = PhaseConcMod(dat, modEG)
    return(dat)
}

openMod <- function(mod, dir, varName, years, modScale, ..., layer = NULL) {
    tempFile = paste(c('temp/', mod, varName, range(years), ..., modScale, '.nc'),
                     collapse = '-')
    cat("\nopening:", mod)
    cat("\n\tinto:", tempFile, "\n")
    if (file.exists(tempFile)) dat = brick(tempFile)
    else {
        files = list.files(paste0(dir, mod, '/'), full.names = TRUE)
        fyr = substr(files, nchar(files)-6, nchar(files)-3)
   
        files = files[apply( sapply(fyr, '==', years), 2, any)]
        dat = layer.apply(files, process.jules.file, 1, varName)
        dat = datConvert(dat, modScale, ...)
        dat = convert_pacific_centric_2_regular(dat)
    
        writeRaster(dat, file = tempFile, overwrite = TRUE)
    }
    if (!is.null(layer)) dat = dat[[layer]]
    return(dat)
}

openObs <- function(obs, Layers, scale, modEG, ..., layer = NULL) { 
    
    tempFile = paste(c('temp/xx', filename.noPath(obs, TRUE),
                     range(Layers), ..., scale, '.nc'),
                     collapse = '-')
    cat("\nopening:", obs)
    cat("\n\tinto:", tempFile)
    
    if (file.exists(tempFile)) obs = brick(tempFile)
    else {
        obs = brick(obs)
        obs = obs[[Layers]]
        obs = datConvert(obs, scale, modEG = modEG[[1]], ...)

        #obs[is.na(modEG)] = NaN
        #browser()
        writeRaster(obs, file = tempFile, overwrite = TRUE)
    }
    if (!is.null(layer)) obs = obs[[layer]]

    return(obs)
}

minusStandard <- function(r1, r2) r1 - r2
annual_average_NME <- function(mods, obss, fname, cols, limits, dcols, dlimits, 
                                ..., FUN = NME, nullFUN = null.NME, diffFUN = minusStandard,
                                    plotFun = plotStandardMap, legendFun = StandardLegend) {   
    png(fname, height = 4 * 15 *(length(obss)+2)/36, width = 4*(length(mods)+1),
        units = 'in', res = 300)
    par(mfcol = c(length(obss)+2, length(mods)+1), mar = rep(0, 4))
    plot.new()

    openYears <- function(year) 
        lapply(mods, openMod, dir, varName, year, modScale,...)

    #mods = openMod(dir, mod, varName, year, TRUE) * modScale
    mods_names = mods
    mods = lapply(years, openYears)
    
    obs_names = names(obss) #sapply(obss, filename.noPath, TRUE)
    
    obss = mapply(openObs, obss, obsLayers, obsScale,
                  MoreArgs = list(modEG = mods[[1]][[1]][[1]], ...))
    
    mapply(plotStandardMap, obss, obs_names,
           MoreArgs = list(cols = cols, limits = limits))

    legendFun(cols, limits, dat = obss[[1]], oneSideLabels = FALSE)    
    
    plotMod <- function(mod, name) {
        
        plotFun(mod[[1]], '', cols = cols, limits = limits)
        mtext(name, side = 3, adj = 0.1, line = -1)
        plotDiff <- function(obs, md) 
            plotFun(diffFUN(md[[1]], obs), '', cols = dcols, limits = dlimits)
        mapply(plotDiff, obss, mod)
        
        legendFun(dcols, dlimits, dat = mod[[1]][[1]]-obss[[1]], oneSideLabels = FALSE)  
    }
    mods_switch = lapply(1:2, function(i) lapply(mods, function(m) m[[i]]))
    mapply(plotMod, mods_switch, mods_names)
    dev.off()
    NMEout <- function(mod, obs) {
        
        scores = FUN( obs, mod[[1]])
        if (class(scores) != "numeric") scores = score(scores)
        null = nullFUN(obs, n = 5)
        
        null = c(summary(null)[1:2], summary(null)[3] + c(-1, 1) * summary(null)[4])
        if (length(scores) == 1) out = t(c(null, scores))
            else out = t(sapply(scores, function(i) c(null, i)))
        return(out)
    }
    
    score = mapply(function(obs, mod) lapply(mod, NMEout, obs), obss, mods)
    scoreI = score[1,]
    if (nrow(score) == 2)
        scoreO = mapply(function(i, j) cbind(i, j[,ncol(j)]), scoreI, score[2,],
                        SIMPLIFY = FALSE)
    else browser()
    scoreO = do.call(rbind, scoreO)
    if (nrow(scoreO) == (3* length(obs_names)))
            rownames(scoreO) = paste(rep(obs_names, each = 3), ' - step', 1:3)
    else rownames(scoreO) = obs_names
            
    colnames(scoreO) = c("Median null", "Mean null", "RR null lower", "RR null upper",
                         mods_names)
    return(scoreO)
}

#score_aa = annual_average_NME(mods, obss, 'figs/burnt_area_aa.png', 
#                               cols_aa, limits_aa, dcols_aa, dlimits_aa, TRUE)

#score_tr = annual_average_NME(mods, obss, 'figs/burnt_area_tr.png',
#                              cols_tr, limits_tr, dcols_tr, dlimits_tr, 
#                              FALSE, TRUE)

MPD_only <- function(x, y, w = NULL) 
    out = MPDonly(x, y)[[1]]
    
null.MPD_only <- function(x, n = 5) 
    null.FUN(x, FUN = MPDonly, n = 5)


score_ph = annual_average_NME(mods, obss, 'figs/burnt_area_ph.png', 
                              cols_phase, limits_phase, dcols_phase, dlimits_phase,
                              FALSE, FALSE, TRUE, FUN = MPD_only, nullFUN = null.MPD_only,
                              diffFUN = phaseDiff,
                              legendFun = SeasonLegend, layer = 1)

score_cn = annual_average_NME(mods, obss, 'figs/burnt_area_cn.png', 
                              cols_conc, limits_conc, dcols_conc, dlimits_conc,
                              FALSE, FALSE, TRUE, layer = 2)

score_md = annual_average_NME(mods, obss, 'figs/burnt_area_md.png', 
                              cols_modal, limits_modal, dcols_modal, dlimits_modal,
                              FALSE, FALSE, TRUE, layer = 3)



