source("cfg.r")
modFracVar = 'frac'
extent = c(-180, 180, -90, 90)
years = 2001:2002
ModLevels = list(c(1:5, 12:13), c(6:11), 16)

sims = paste0("../jules_outputs/", c("ALL_u-bk886_isimip_0p5deg_origsoil_dailytrif/",
                                     "RUNC20C_u-by276_isimip_0p5deg_origsoil_dailytrif_fire/"))
names(sims) = c("Fire off", "Fire on")
dir_dat = "../fireMIPbenchmarking/data/benchmarkData/"
obss = list(CCI = "vegfrac_refLC_refCW.nc",
            IGBP = "vegfrac_igbp.nc",
            VCF = c("treecover2000-2014.nc", "nontree2000-2014.nc"))
obss = lapply(obss, function(i) paste0(dir_dat, i))

vars = c("trees", "herb", "bare") 

realmNames = c("Global", "Australia", "Southern\nAfrica", "Northern\nAfrica", "SE\nAsia",
               "Southern\nAmerica", "Temperate\nAmerica", "Boreal\nAmerica",
               "Temperate\nEurasia", "Boreal\nEurasia")

cols = c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d',
         '#238443','#006837','#004529')

rmask = raster('../UKESM-land-evaluation/outputs/realms.nc')
rmask = raster::resample(rmask, raster('data/seamask.nc'))
rmask = convert_regular_2_pacific_centric(rmask)
meanBrick <-function(...) mean(brick(...))
#files = list.files(dir, full.names=TRUE)

forDataset <- function(id, idn, layers = NULL, modid = 1) {   
    forRegion <- function(rid, name) {        
        forModel <- function(dat, modname, add = FALSE) {      
             modname = gsub("_noFire", "", modname)
            if (modname != '') idn = paste0(modname, '\n', idn)
            tree = dat[[1]]
            herb = dat[[2]]
            if (grepl("VCF", idn)) {
                tree = tree
                bare = 1-tree-herb
            } else bare = dat[[3]]
            
    #open <- function(f) {
    #    if (length(f) > 1) out = mean(layer.apply(f, meanBrick))
    #    else out = raster(f)
    #    return(out)
    #}
    #tree = open(tree)
    #herb = open(herb)
    #bare = open(bare)   
           
            tot = tree + herb + bare
        
            tree = tree/tot
            herb = herb/tot
            bare = bare/tot
            rmask = round(raster::resample(rmask,tree))
            
            mask = !is.na(tot + rmask)
        
            if (!is.na(rid)) mask = mask & ( rmask== rid)
            hs = herb[mask]
            ts = tree[mask]
            bs = 1 - ts - hs
            trianglePlot(bs, ts, cols, name, add = add)

            if (name == realmNames[1]) mtext(side = 3, adj = 0.9, idn, line = -0.67)
            if (tail(realmNames, 1) == name)
                text(x = 1, y = -0.19, 'Bare (%)', xpd = TRUE, adj = 1, xpd = NA)
            if (id == obss[1]) {
                mtext(side = 2, adj = 0.9, name, line = -2)
                text(x = 0.06, y = 0.38, 'Herb (%)', xpd = TRUE, srt = 60, adj = 1)
            }
            if (id == tail(obss, 1)) 
                text(x = 0.78, y = 0.67, 'Wood (%)', xpd = TRUE, srt = -60, adj = 1)
        
            return(cbind(ts, hs, bs))
        }
        if (is.null(layers)) {
            c(mods, mod):= openModsFromDir(id, years, ModLevels, extent)
            #mapply(forModel, mod[[1]], mods[[1]], c(F, rep(F, length(mod) -1)))
            forModel(mod[[1]][[modid]], mods[[1]][[modid]])
        } else {
            if (is.list(layers)) {                
                dat = brick(id)
                dat = layer.apply(layers, function(i) sum(dat[[i]]))
            } else {
                dat = layer.apply(id, function(i) mean(brick(i)[[layers]]))
            }
            if (max.raster(dat, na.rm = TRUE) > 1) dat = dat/100
            forModel(dat, '')    
        }
    }
    mapply(forRegion, c(NaN, 1:9), realmNames)
}
obsLayers = list(list(c(1:2, 5), c(3:4), 8), list(c(1:2, 5), c(3:4), 8), 1:60)
png("figs/VegDistTriangle.png", width = 7.2*11/4, height = 1.1*7.2 * sqrt(3) * 0.5 * 10/4,
    res = 300, units = 'in')        
    par(mfcol = c(10, 11), mar = rep(1, 4), oma = c(0.35, 0, 1, 0.3))
    mapply(forDataset, obss, names(obss), layers = obsLayers, SIMPLIFY = FALSE)
    lapply(1:4,  function(i) 
           mapply( forDataset, sims, names(sims),modid = i, SIMPLIFY = FALSE))
dev.off()
