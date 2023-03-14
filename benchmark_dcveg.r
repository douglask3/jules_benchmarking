source("cfg.r")

mod_dirs = c("/hpc/data/d01/hadcam/jules_output/", "/hpc/data/d05/cburton/jules_output/")

mods = c('ALL_u-bk886_isimip_0p5deg_origsoil_dailytrif/', 'u-cf137/')
isimip_mods = c("GFDL-ESM2M", "HADGEM2-ES", "IPSL-CM5A-LR", "MIROC5")

mods = paste0(rep(mods, each = 4), isimip_mods, '/')
mod_dirs = rep(mod_dirs, each = 4)

stream = '.ilamb.'

varName = 'cv'
modScale = 12#60*60*24*360

extent = c(-180, 180, -90, 90)

years = c(2010, 2018)


obs_files = paste0('../forest_emissions_attribution/data/raw/biomass_cci/ESACCI-BIOMASS-L4-AGB-MERGED-100m-', c('2010', '2018'), '-fv2.0_regridded_JULES_01deg.nc')

openRegrid <- function(file) {
    dat = raster(file)
    aggregate(dat, fact = 5)
}

mod2obsGrid <- function(mod) {
    mod = addLayer(mean(mod[[1:12]]), mean(mod[[13:24]]))
    mod = raster::resample(mod, dobs)
}
diffMod <- function(mod) 
    mod[[2]] - mod[[1]]

if (T) {
mods = mapply(openMod, mods, mod_dirs, 
               MoreArgs = list(varName, years, modScale, levels = 1, stream = stream))

obs = layer.apply(obs_files, openRegrid)/12
dobs = obs[[1]] - obs[[2]]
mods = lapply(mods, mod2obsGrid)
dmods = lapply(mods, diffMod)
dobs[is.na(dmods[[1]])] = NaN#obs = layer.apply(obs_files, openRegrid)

Dnme = lapply(dmods, function(i) NME(dobs, i, w = raster::area(dobs)))
nme  = lapply(1:2, function(i) 
                lapply(mods, function(mod) NME(obs[[i]], mod[[i]], w = raster::area(obs))))

Dnull = null.NME(dobs, n = 5, w = raster::area(dobs))
null = lapply(1:2, function(i) null.NME(obs[[i]], n = 5, w = raster::area(obs[[i]])))
}
browser()
makeRow <- function(nme, null, nr = 1) {
    nulls = standard.round(summary(null)[1:4], 2)
    nulls[3] = paste(nulls[3], '+/-', nulls[4])
    nulls = nulls[1:3]
    scrs = standard.round(sapply(nme, function(i) score(i)[nr]))
    c(nulls, scrs)
}

#out = rbind(makeRow(Dnme, Dnull), 
#        makeRow(nme[[1]], null[[1]]), 
#        makeRow(nme[[1]], null[[1]], 2), 
#        makeRow(nme[[1]], null[[1]], 3), 
#        makeRow(nme[[2]], null[[2]]), 
#        makeRow(nme[[2]], null[[2]], 2), 
 #       makeRow(nme[[2]], null[[2]], 3))
#outr = do.call(rbind, out)

regions = raster("../ConFIRE_attribute/data/GFEDregions.nc")
regions = round(raster::resample(regions, dobs))

forRegion <- function(id) {
    mask = regions!=id
    dobs[mask] = NaN
    area = raster::area(dobs, na.rm = TRUE)
    sarea = sum(area[], na.rm = TRUE)
    
    change <- function(r) {
        r[mask] = NaN
        sum((r*area)[], na.rm = TRUE)/sarea
    }
    cobs = change(dobs)
    cmod = sapply(dmods, change)
    return(c(sarea, cobs, cmod))
}

#out = sapply(unique(regions[!is.na(regions)]), forRegion)

out[2,] = -out[2,]
#out[2,] = out[2,] /2
xrange = range(out[-1,])
plot(xrange, xrange, type = 'n', xlab = '', ylab = '')
grid()

addRegion <- function(xys) {
    points(rep(xys[2], 4), xys[3:6], pch = 19, col = "blue")
    points(rep(xys[2], 4), xys[7:10], pch = 19, col = "red")
}

apply(out, 2, addRegion)

bench <- function(x) 
    score(NME(out[2,], x, w = out[1,]))[1]

nme_regions = apply(out[-c(1:2),], 1, bench)
