
PhaseConcMod <- function(dat, modEG = NULL) {
    clim = layer.apply(1:12, function(mn) mean(dat[[seq(mn, nlayers(dat), by = 12)]]))
    if (!is.null(modEG)) clim = raster::resample(clim, modEG)
    PC = PolarConcentrationAndPhase(clim, 'months')
    
    Modal = Modalise(clim)
    return(addLayer(PC, Modal))
}

