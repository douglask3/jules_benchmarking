source("cfg.r")





regions = raster('../ConFIRE_attribute/data/GFEDregions.nc')

region_names = c("BONA", "TENA", "CEAM", "NHSA", "SHSA", "EURO", "MIDE", "NHAF", 
                 "SHAF", "BOAS", "CEAS", "SEAS", "EQAS", "AUST")


obs_files = paste0("/data/dynamic/dkelley/fireMIPbenchmarking/data/benchmarkData/",
              c("GFED4s_v2.nc", "GFED4.nc", "MCD45.nc", "meris_v2.nc",#../GFED4.fBA.r0d5.1995.2013.nc"
                "MODIS250_q_BA_regridded0.5.nc"))



annualise <- function(dat) { 
    FUN <- function(yr) 
        sum(dat[[((yr-1)*12+1):(yr*12)]])
    
    out = layer.apply(1:(nlayers(dat)/12), FUN)
    out = raster::resample(out, regions)
    
}

regionalise <- function(dat) {
    regionSumms <- function(region) {
        if (region >0) ys = datr[regions == region]
        else ys = datr[!is.na(regions)]
        
        out = apply(ys, 2 , sum, na.rm = TRUE)
    }
    datr = dat * raster::area(dat)
    lapply(0:max(regions[], na.rm = TRUE), regionSumms)
}

openObs <- function(file) {
    dat = annualise(brick(file)) * obs_scale
    if (grepl('MCD45', file)) dat = dat[[-nlayers(dat)]]
    if(max(dat[[1]][dat[[1]]<9E9], na.rm = TRUE) > 50) dat = dat/100
    dat[dat>9E9] = NaN
    reg_dat = regionalise(dat)
    return(list(dat, reg_dat))
}

openMod <- function(model) {
    openPeriod <- function(dir, layer) {
        file = list.files(paste0(dir, '/', model, '/'), full.names = TRUE)
        dat = brick( file[grepl(varname, file)])
        dat = dat[[layer]]
        
        dat = annualise(dat)*mod_scale
        return(dat)
    }
    dat = mapply(openPeriod, dirs, layers)
    dat = addLayer(dat[[1]], dat[[2]])
    
    reg_dat = regionalise(dat)
    return(list(dat, reg_dat))
}



output4region <- function(region, name, mods, obss) {
    
    modsR = sapply(mods, function(mod) mod[[2]][[region+1]])
    modsR = cbind(modsR, apply(modsR, 1, mean))
    obssR = lapply(obss, function(obs) obs[[2]][[region+1]])
    
    years = seq(mod_start, length.out = dim(modsR)[1])
    
    fillBlanks <- function(obs, start) {
        out = rep(NaN, length(years))
        start = which(years == start)
        out[seq(start, length.out = length(obs))] = obs
        return(out)
    }
    obssR_fill = mapply(fillBlanks, obssR, obs_start)


    out = cbind(years, obssR_fill, modsR) 
    colnames(out) = c('years', paste0('Obs-', names(obs_files)),
                               paste0('Sim-', c(models, 'enemble')))
    dir = paste0("outputs/", varname, dir_lab)
    makeDir(dir)
    write.csv(out, file = paste0(dir,  name, '-', region, ".csv"))

}

regionArea = sapply(1:14, function(i)
                    sum(raster::area(regions)[regions == i]))
totArea = sum(regionArea)

MM_single <- function(x, y, ...)
    MM(cbind(x, 1-x), cbind(y, 1-y), ...)

null.MM_single <- function(x, ...)
    null.MM(cbind(x, 1-x), ...)


plot4Region <- function(region) {
    modsR = lapply(mods, function(mod) mod[[2]][[region]])
    obssR = lapply(obss, function(obs) obs[[2]][[region]])
    ylim = range(0, c(unlist(modsR), unlist(obssR)))
    plot(xlim, ylim, type = 'n')

    linesDat <- function(mod, lty, start = xlim[1], col = 'red') 
        lines(start:(start+length(mod)-1), mod, col = col, lwd = 2, lty = lty)
    
    mapply(linesDat, modsR, 1:4)
    mapply(linesDat, obssR, c(1, 1, 1:4), obs_start, col = 'blue')
}

RUN <- function() {
    mods = lapply(models, openMod)
    obss = lapply(obs_files, openObs)
    
    mapply(output4region, 0:14, c('GLOB', region_names), MoreArgs = list(mods, obss))
    #plot4Region(2)
    return(list(obss, mods, obs_start, mod_start))
}

dirs = c("../ConFIRE_ISIMIP/INFERNOonOff/Global/historic_on/", "../ConFIRE_ISIMIP/INFERNOonOff/Global/RCP6.0_on/")

models = list.files(dirs[1])

layers = list(1621:1728, 1:(15*12))

obs_files =  c(paste0("/home/h01/cburton/GitHub/ISIMIP3a/Observations/", 
                     c("GFED4.1s_Burned_Fraction.nc", "FireCCI5.1_Burned_Fraction.nc", "GFED500m_Burned_Percentage.nc")),
            paste0("/data/dynamic/dkelley/fireMIPbenchmarking/data/benchmarkData/",
                  c("GFED4.nc", "MCD45.nc", "meris_v2.nc")))
names(obs_files) = c("GFED4s", "fire_CCI", "GFED500", "GFED4", "MCD45", "meris")
obs_start = c(2001, 2002, 1998, 1998, 2001, 2006)
mod_start = 1997

xlim = c(1997, 2020)

varname = 'burnt_area'
mod_scale = 60*60*24*30
obs_scale = 1
dir_lab = "_obs_isimip2b/"

obsLayers = NULL
if (T) {
burnt_area = RUN()

obs_files =  c(VCF = "/data/dynamic/dkelley/fireMIPbenchmarking/data/benchmarkData/treecover2000-2014.nc")
obs_start = c(2000)
mod_start = 2000

layers = list(1645:1728, 1:(10*12))
xlim = c(1999, 2015)

varname = 'trees'
mod_scale = 1/12
obs_scale = 1/(0.8*12)
dir_lab = "_obs_isimip2b_fireOn/"
trees_on = RUN()

varname = 'tallTrees'
tallTrees_on = RUN()

dirs = c("../ConFIRE_ISIMIP/INFERNOonOff/Global/historic_off/", "../ConFIRE_ISIMIP/INFERNOonOff/Global/RCP6.0_off/")

dir_lab = "_obs_isimip2b_fireOff/"

varname = 'trees'
trees_off = RUN()


tallTrees_off = RUN()
}
dirs = c("../ConFIRE_ISIMIP/INFERNOonOff/Global/historic_on/", "../ConFIRE_ISIMIP/INFERNOonOff/Global/RCP6.0_on/")
obs_files =  c(GFAS = "/data/dynamic/dkelley/fireMIPbenchmarking/data/benchmarkData/GFAS.nc")
varname = 'fire_emissions_veg'
dir_lab = "_obs_isimip2b/"

obs_start = 2000
fire_emissions = RUN()
obs_scale = 60*60*24*30
save(burnt_area, trees_on, tallTrees_on, trees_off, tallTrees_off, fire_emissions, regionArea, file = "outputs/regional_process_ts.Rd")
