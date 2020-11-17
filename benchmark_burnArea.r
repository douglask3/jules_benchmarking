library(raster)
library(gitBasedProjects)
library(rasterExtras)
sourceAllLibs("libs/")
sourceAllLibs('../benchmarkmetrics_github/benchmarkMetrics/R/')

dir = '../jules_outputs/'

mods = c('u-by849', 'u-by851')
varName = 'burnt_area_gb'
modScale = 60*60*24*360

obs = "../fireMIPbenchmarking/data/benchmarkData/GFED4s_v2.nc"
years = 1997:2014
obsLayers = 1:12
obsScale = 12

openMod <- function(dir, mod, varName, years, annualAver = FALSE) {
    files = list.files(paste0(dir, mod, '/'), full.names = TRUE)
    fyr = substr(files, nchar(files)-6, nchar(files)-3)
   
    files = files[apply( sapply(fyr, '==', years), 2, any)]
    dat = layer.apply(files, process.jules.file, 1, varName)
    if (annualAver) dat = mean(dat)
    dat = convert_pacific_centric_2_regular(dat)
    return(dat)
}


annual_average_NME <- function(mod, obs) {
    mod = openMod(dir, mod, varName, years, TRUE) * modScale
    
    obs = brick(obs)
    obs = obs[[obsLayers]]
    obs = mean(obs) * obsScale
    obs = raster::resample(obs, mod[[1]])
    obs[is.na(mod)] = NaN
    scores = NME(mod, obs)
    null = null.NME(obs, n = 5)
    null = c(summary(null)[1:2], summary(null)[3] + c(-1, 1) * summary(null)[4])
    out = t(sapply(score(scores), function(i) c(null, i)))
    return(out)
}

score = annual_average_NME(mods[1], obs)
