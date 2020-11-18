source("libs/PolarConcentrationAndPhase.r")
library(plotrix)

phaseDiff <- function(x, mn = NULL) {
    if (is.null(mn)) {
        if (is.na(x)) return(NA)
    } else {
        if (is.na(x) || is.na(mn)) return(NA)
        x = x - mn
    }
    if (abs(x) <= 6) return(x)
    if (x <(-6)) return(x+12)
    if (x >  6 ) return(x-12)
}
plotConPhase <- function(r, let = c('', ''), addLegend = FALSE,
                         phase_colsi = phase_cols,
                         conc_colsi = conc_cols,
                         phase_limsi = phase_lims,   conc_limsi =conc_lims,
                         regions = NULL) {
    
    if (is.list(r) && length(r) == 1) {
        pc = PolarConcentrationAndPhase.RasterBrick(r[[1]], phase_units = "months")
        
        e = list(NULL, NULL)
    } else {
        pcs = lapply(r, PolarConcentrationAndPhase.RasterBrick, phase_units = "months")
        phases = layer.apply(pcs, function(i) i[[1]])
        phase = phases[[1]]
        phase[] = apply(phases[], 1, median)
        #browser()
        dphase = phases-phase
        dphase = layer.apply(dphase, function(i) {i[] = sapply(i[], phaseDiff); i})
        dphase = max(abs(dphase))

        concs = layer.apply(pcs, function(i) i[[2]])
        conc = mean(concs)
        dconc = eFromRange(concs)

        pc = list(phase, conc) 
        e = list(dphase, dconc)
    }
    
    plotStandardMap(pc[[1]], e = e[[1]], limits = phase_limsi, cols = phase_colsi, limit_error = c(1, 3))
    addLetLab(let[1])
    if (addLegend) SeasonLegend(phase_limsi, cols = phase_colsi, add = FALSE)
    plotStandardMap(pc[[2]], e = e[[2]], limits = conc_limsi, cols = conc_colsi)
    addLetLab(let[2])
    if (addLegend) conLegend(pc, conc_limsi, conc_colsi)

    if (!is.null(regions)) {
        browser()
    }
}

conLegend <- function(pc, limits = conc_lims, cols = conc_cols) {
    if (is.list(pc) && length(pc) ==2) pc = pc[[2]]
    else if (nlayers(pc) == 2) pc = pc[[2]]
    StandardLegend(limits = limits, cols = cols, extend_max = FALSE,
                        maxLab = 1, dat = pc, add = TRUE, oneSideLabels = FALSE)
}
