source("cfg.r")

source("../ConFIRE_attribute/libs/plotStandardMap.r")
source("../LPX_equil/libs/legendColBar.r")

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
        if (class(ys) == 'matrix') out = apply(ys, 2 , sum, na.rm = TRUE)
            else out = sum(ys, 2, na.rm = TRUE)
    }
    datr = dat * raster::area(dat)
    lapply(0:max(regions[], na.rm = TRUE), regionSumms)
}

openObs <- function(file) {
    print(file)
    dat = brick(file)
    if (nlayers(dat) > 12) dat = annualise(dat)
    dat = dat  * obs_scale
    if (grepl('fraccci', file)) dat = dat * 12 * 0.8
    if (grepl('MCD45', file)) dat = dat[[-nlayers(dat)]]
    if (grepl('avitabile', file)) dat = dat * 45 
    if (grepl("ESACCI-BIOMASS-", file)) dat = mean(dat[[-1]])
    if (grepl('GFAS', file)) dat[is.na(mods[[1]][[1]][[1]])] = NaN
    if(max(dat[[1]][dat[[1]]<9E9], na.rm = TRUE) > 50) dat = dat/100
    dat[dat>9E9] = NaN
    dat = raster::resample(dat, regions)
    reg_dat = regionalise(dat)
    return(list(dat, reg_dat))
}

openMod <- function(model) {
    print(model)
    openPeriod <- function(dir, layer) {
        if (is.null(layer)) return(NULL)
        file = list.files(paste0(dir, '/', model, '/'), full.names = TRUE)
        dat = brick( file[grepl(varname, file)])
        dat = dat[[layer]]
        
        dat = annualise(dat)*mod_scale
        return(dat)
    }
    dat = mapply(openPeriod, dirs, layers)
    testNll = sapply(dat, is.null)
    if (any(testNll) ) dat = dat[!testNll][[1]] else dat = addLayer(dat[[1]], dat[[2]])
    
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
    dir = paste0("outputs/", varname, dir_lab, '/')
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

plotMap <- function(r, txt = '', cols, limits, side = 3, ...) {
    
    letterN <<- letterN + 1
    plotStandardMap(r, cols = cols, limits = limits,y_range = c(-60, 90), #projection = 'robinson',
                     interior = FALSE)
    mtext(side = 1, line = -2, adj = 0.2, letters[letterN])
    mtext(txt, side = side, ...)
}

plotMaps <- function(mods, obss, obsNames, obsStart, scale = 1, units = '%', regional = FALSE,
                     readyStiched = FALSE) {
    
    if (readyStiched) lmods = length(mods[[1]]) else lmods = length(mods)
    lmat = length(obss) + matrix(1:(lmods*(length(obss)+1)), ncol = length(obss)+1)
    
    lmat = rbind(c(0, 1:length(obss)), lmat, max(lmat) +1)
    lmat = rbind(0, lmat, 0); lmat = cbind(0, lmat)
    heights =c(0.1, rep(1, lmods +1), 0.8, 0.1)
    widths =c(0.1, rep(1, ncol(lmat)-1))    
    
    png(paste0("figs/", varname, '_', dir_lab, "regioanl_is_", regional, "_isimip2b_maps.png"), height = sum(heights)*1.2,
        width = sum(widths)*2.5, res = 300, units = 'in')
    
    layout(lmat, heights = heights, widths = widths)    
    par(mar = rep(0, 4))

    if (regional) {
        meanFUN <- function(r, layers = 1:nlayers(r[[1]]), ...) {
            r[[1]] = r[[1]][[1]]
            for (reg in 1:14) 
                r[[1]][regions == reg] = mean(r[[2]][[reg+1]])/regionArea[reg]  

            return(r[[1]])
        }
    } else {
        meanFUN <- function(r, layers = 1:nlayers(r[[1]]), ...) mean(r[[1]][[layers]], ...)
    }
    if (readyStiched) {
        absPlot_fun <- function(i) scale*meanFUN(i[[1]], na.rm = TRUE)
        diff_from_obs <- function(obs, mod) 
            lapply(mod, function(i) meanFUN(i) - meanFUN(obs[[1]]))
                   
        mods_diff = unlist(mapply(diff_from_obs, obss, mods, SIMPLIFY = FALSE))
    } else {
        diff_from_obs <- function(obs, stYr) {
            lys  =  1:length(obs[[2]][[1]]) + stYr - mod_start
            lapply(mods, function(i)
                            scale*(meanFUN(i, lys, na.rm = TRUE) - meanFUN(obs) ))
    
        }
        absPlot_fun <- function(i) scale * meanFUN(i, na.rm = TRUE)
        mods_diff = unlist(mapply(diff_from_obs, obss, obsStart, SIMPLIFY = FALSE))
    }
    mods_plot = lapply(mods, absPlot_fun)         
    obss_plot = lapply(obss, absPlot_fun)
    browser()
    mapply(plotMap, obss_plot, obsNames, MoreArgs = list(cols, limits, side = 3, line = -0.5))
    mapply(plotMap, mods_plot, models, 
           MoreArgs = list(cols, limits, side = 2, line = -0.5)) 
    lapply(mods_diff, plotMap, '', dcols, dlimits) 
    plot.new()
    legendFUN <- function(cols, limits, extend_min, xx) {
        if (extend_min) { maxLab = ''; minLab = ''}
        legendColBar(cols = cols, limits = limits, maxLab = maxLab, minLab = minLab,
                     transpose = TRUE, oneSideLabels = TRUE, 
                     add = TRUE, extend_min = extend_min,  extend_max = maxLab == '',
                     xx = xx, yy = c(0, 1))
    }
    legendFUN( cols,  limits, FALSE, c(0.66, 0.95))
    legendFUN(dcols, dlimits, TRUE, c(0.17, 0.46))  
    mtext(side = 1, line = -1, '%')
    dev.off() 
}


RUN <- function() {
    mods <<- lapply(models, openMod)
    obss <<- lapply(obs_files, openObs)
    
    mapply(output4region, 0:14, c('GLOB', region_names), MoreArgs = list(mods, obss))

    letterN <<- 0
    for (reg in c(T, F)) 
    plotMaps(mods, obss[plotObs],
             names(obs_files)[plotObs], obs_start[plotObs], plotScale, regional = reg)
    return(list(obss, mods, obs_start, mod_start))
}

#########################
## Burnt Area          ##
#########################
dirs = c("../ConFIRE_ISIMIP/INFERNOonOff/Global/historic_on/", 
         "../ConFIRE_ISIMIP/INFERNOonOff/Global/RCP6.0_on/")


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

plotObs = 2:3
plotScale = 100
cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
dcols = c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')

limits = c(0, 1, 2, 5, 10, 20, 40, 60)
dlimits = c(-20, -10,-5, -2, -1,  1, 2, 5, 10, 20)
maxLab = minLab = ''

burnt_area = RUN()

#########################
##Tree Cover          ##
#########################

obs_files =  c(VCF = "/data/dynamic/dkelley/fireMIPbenchmarking/data/benchmarkData/treecover2000-2014.nc", CCI = "data/fraccci/TREE.nc")
obs_start = c(2000, 2010)
mod_start = 2000

layers = list(1645:1728, 1:(10*12))
xlim = c(1999, 2015)

varname = 'trees'
mod_scale = 1/12
obs_scale = 1/(0.8*12)
dir_lab = "_obs_isimip2b_fireOn"


plotObs = 1:2
plotScale = 100
cols = c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476','#41ab5d','#238b45','#006d2c','#00441b')
dcols = c('#762a83','#9970ab','#c2a5cf','#e7d4e8','#f7f7f7','#d9f0d3','#a6dba0','#5aae61','#1b7837')

limits = seq(10, 90, by = 10)
dlimits = seq(-40, 40, by = 10)
maxLab = 100
minLab = 0

trees_on = RUN()

varname = 'tallTrees'
#tallTrees_on = RUN()

dirs = c("../ConFIRE_ISIMIP/INFERNOonOff/Global/historic_off/", "../ConFIRE_ISIMIP/INFERNOonOff/Global/RCP6.0_off/")

dir_lab = "_obs_isimip2b_fireOff"

varname = 'trees'
#trees_off = RUN()


varname = 'tallTrees'
#tallTrees_off = RUN()

dirs = c("../ConFIRE_ISIMIP/INFERNOonOff/Global/historic_on/", "../ConFIRE_ISIMIP/INFERNOonOff/Global/RCP6.0_on/")
obs_files =  c(GFAS = "/data/dynamic/dkelley/fireMIPbenchmarking/data/benchmarkData/GFAS.nc")
varname = 'fire_emissions_veg'
dir_lab = "_obs_isimip2b"

obs_start = 2000
obs_scale = 60*60*24*30

plotObs = 1
plotScale = 12
cols = c('#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')
dcols = rev(c('#b35806','#e08214','#fdb863','#fee0b6','#f7f7f7','#d8daeb','#b2abd2','#8073ac','#542788'))

limits = c(0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 4)
dlimits = c(-1.2, -0.8, -0.4, -0.2, -0.1, 0.1, 0.2, 0.4, 0.8, 1.2)
maxLab = ''
minLab = ''

#fire_emissions = RUN()
#save(burnt_area, trees_on, tallTrees_on, trees_off, tallTrees_off, fire_emissions, regionArea, file = "outputs/regional_process_ts.Rd")


obs_files = c(paste0("/data/dynamic/dkelley/fireMIPbenchmarking/data/benchmarkData/", 
                     c("avitabile_carbon_veg_05.nc", "Carvalhais.cVeg_50.360.720.1-1.nc")),
              paste0("outputs/", 
                     c("ESACCI-BIOMASS-L4-AGB-MERGED-0.5Degree-2010-fv3.0-ensembles-10.nc",
                       "ESACCI-BIOMASS-L4-AGB-MERGED-0.5Degree-2018-fv3.0-ensembles-10.nc")))

varname = 'cveg'
names(obs_files) = c('Avitabile', 'Carvalhais', 'CCI Carbon 2010', 'CCI Carbon 2018')
years = c(1997, 1997, 2010, 2018)
obs_scale = 1
mod_scale = 1
layerss = list(list(1621:1728, NULL), list(1621:1728, NULL), 
               list(NULL, 49:60), list(NULL, 289:300))

plotObs = c(1:2, 4)

plotScale = 1
cols = c('#fff7fb','#ece2f0','#d0d1e6','#a6bddb','#67a9cf','#3690c0','#02818a','#016c59','#014636')
dcols = c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e')

limits = c(0, 1,5, 10, 20, 50, 100, 150, 200)
dlimits = c(-100, -50, -20, -10, -5, -1, 1, 5, 10, 20, 50, 100)
maxLab = ''
minLab = ''


FUN <- function(obs_file, yr, layer) {
    obss = lapply(obs_file, openObs)

    obs_start <<- yr
    mod_start <<- yr
    layers <<- layer
    mods = lapply(models, openMod)
    obss = lapply(obs_file, openObs)
    letterN <<- 0
    
    return(list(mods, obss))
}

#vcarbon = mapply(FUN, obs_files, years, layerss)

#save(burnt_area, trees_on, tallTrees_on, trees_off, tallTrees_off, 
#     fire_emissions, vcarbon, regionArea, file = "outputs/regional_process_ts.Rd")

#plotMaps(vcarbon[1,plotObs], vcarbon[2,plotObs],
#             names(obs_files)[plotObs], obs_start[plotObs], plotScale, regional = FALSE,
#        readyStiched = TRUE)


dvcarbon = vcarbon[,3:4]
dobss = list(mean(brick('outputs/ESACCI-BIOMASS-L4-AGB-MERGED-0.5Degree-CHANGE-2018-2010-fv3.0-ensembles-10.nc')[[-1]][[1:10]]))
dmods = list(lapply(1:4, function(i) 
                    list(mean(dvcarbon[[1,2]][[i]][[1]]) - mean(dvcarbon[[1,1]][[i]][[1]]))))

varname = 'dcveg'
cols = c('#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221')
limits = c(-20, -10, -5, -2, -1, 1, 2, 5, 10, 20)
dlimits = limits
lettersN <<- 0
plotMaps(dmods, dobss, 'AGB change 2018-2010', 2010, plotScale, 
         readyStiched = TRUE, regional = FALSE)

browser()
frac_cci = brick("../fireMIPbenchmarking/data/benchmarkData/vegfrac_refLC_refCW.nc")
frac_cci = convert_pacific_centric_2_regular(frac_cci)
frac_cci = sum(raster::resample(frac_cci, regions)[[1:2]])

dirss = list(c("../ConFIRE_ISIMIP/INFERNOonOff/Global/historic_off/", "../ConFIRE_ISIMIP/INFERNOonOff/Global/RCP6.0_off/"),
             c("../ConFIRE_ISIMIP/INFERNOonOff/Global/historic_on/" , "../ConFIRE_ISIMIP/INFERNOonOff/Global/RCP6.0_on/" ))

mod_scale = 1/12
mod_start = 2000
layers = list(NULL, 49:60)
mods = lapply(models, openMod)
varname = 'tallTrees'
#obsLayers = list(1:2, 5, 3:4, 8)
#cover_cci = layer.apply(obsLayers, function(i) sum(frac_cci[[i]]))


obs = unlist(regionalise(frac_cci)[-1])/regionArea
compare <- function(dir) {
    dirs <<- dir
    mods = mods0 = lapply(models, openMod)
    mods = lapply(mods, function(mod) unlist(mod[[2]][-1])/regionArea)

    out = sapply(mods, function(i) score(NME(obs, i, w = regionArea)))
    return(out)    
}
outs = lapply(dirss, compare)
outs = do.call(rbind, outs)
#mods = lapply(models, openMod)
nulls = summary(null.NME(obs, w = regionArea))
nulls = nulls[c(1, 2, 3, 3)] + c(0, 0, -1, 1) * nulls[4]
outs = cbind(t(matrix(nulls, nrow = 4, ncol = nrow(outs))), outs)
write.csv('outs', file = 'outputs/cci_veg_frac_region_isimip2b_nme.csv')
