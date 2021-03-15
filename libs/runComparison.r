

runComparison <- function(mods, obss, fname, cols, limits, dcols, dlimits, 
                                ..., FUN = NME, nullFUN = null.NME, diffFUN = minusStandard,
                                plotFun = plotStandardMap, legendFun = StandardLegend,
                                extend_max = TRUE, extend_min = FALSE,
                                dextend_max = TRUE, dextend_min = TRUE) {   
    png(fname, height = 4 * 15 *(length(obss)+2)/36, width = 4*(length(mods)+1),
        units = 'in', res = 300)
    par(mfcol = c(length(obss)+2, length(mods)+1), mar = rep(0, 4))
    plot.new()

    openYears <- function(year) 
        lapply(mods, openMod, dir, varName, year, modScale,...)

    #mods = openMod(dir, mod, varName, year, TRUE) * modScale
    mods_names = mods
    mods = lapply(years, openYears)
    modsEG = mods[[1]][[1]][[1]]
    mods = lapply(mods, function(mod) lapply(mod, function(m) raster::resample(m, modsEG)))
    obs_names = names(obss) #sapply(obss, filename.noPath, TRUE)
    
    obss = mapply(openObs, obss, obsLayers, obsScale,
                  MoreArgs = list(modEG = modsEG, ...))

    mapply(plotFun, obss, obs_names,
           MoreArgs = list(cols = cols, limits = limits))

    legendFun(cols, limits, dat = obss[[1]], oneSideLabels = FALSE,
              extend_max = extend_max, extend_min = extend_min)    
    
    plotMod <- function(mod, name) {
        
        plotFun(mod[[1]], '', cols = cols, limits = limits)
        mtext(name, side = 3, adj = 0.1, line = -1)
        plotDiff <- function(obs, md) 
            plotFun(diffFUN(md[[1]], obs), '', cols = dcols, limits = dlimits)
        mapply(plotDiff, obss, mod)
        
        legendFun(dcols, dlimits, dat = diffFUN(mod[[1]][[1]],obss[[1]]), oneSideLabels = FALSE,
                  extend_max = dextend_max, extend_min = dextend_min)  
    }
    mods_switch = lapply(1:length(mods[[1]]), function(i) lapply(mods, function(m) m[[i]]))
    
    mapply(plotMod, mods_switch, mods_names)
    dev.off.gitWatermark()
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
    scoreO = score[1,]
    if (nrow(score) == 2)
        scoreO = mapply(function(i, j) cbind(i, j[,ncol(j)]), scoreO, score[2,],
                        SIMPLIFY = FALSE)
    else {
        for (id in 2:nrow(score)) {
            scoreO = mapply(function(i, j) cbind(i, j[,ncol(j)]), scoreO, score[id,],
                        SIMPLIFY = FALSE)
        }
    }
        
    scoreO = do.call(rbind, scoreO)
    if (nrow(scoreO) == (3* length(obs_names)))
            rownames(scoreO) = paste(rep(obs_names, each = 3), ' - step', 1:3)
    else rownames(scoreO) = obs_names
            
    colnames(scoreO) = c("Median null", "Mean null", "RR null lower", "RR null upper",
                         mods_names)
    return(scoreO)
}
