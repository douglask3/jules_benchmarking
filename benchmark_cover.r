graphics.off()
source("cfg.r")

mod_dirs = c("/hpc/data/d01/hadcam/jules_output/", "/hpc/data/d05/cburton/jules_output/")

mods = c('ALL_u-bk886_isimip_0p5deg_origsoil_dailytrif/', 'u-cf137_TEST-LowMort/')
isimip_mods = c("GFDL-ESM2M", "HADGEM2-ES", "IPSL-CM5A-LR", "MIROC5")

mods = paste0(rep(mods, each = 4), isimip_mods, '/')
mod_dirs = rep(mod_dirs, each = 4)

varName = 'frac'
modScale = 1#60*60*24*360
modFracVar = 'frac'

stream = '.ilamb.'

dir_dat = "../fireMIPbenchmarking/data/benchmarkData/"
obss = list("vegfrac_refLC_refCW.nc", c("treecover2000-2014.nc", "nontree2000-2014.nc"))
obss = lapply(obss, function(i) paste0(dir_dat, i))
names(obss) = c("cci_lifeForms", "VCF_lifeForms")
years = list(2010:2010, 2001:2014)
obsLayers = list(list(1:2, 3:4, 5, 8), 1)
obsScale = c(1, 1)

extent = c(-90, -30, -60, 15)
extent = c(-180, 180, -90, 90)
limits_aa = c(0, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)*10                
cols_aa =c('#fffbf7','#ecf0e2','#d0e6d1','#a6dbbd','#67cfa9','#36c090',
            '#028a81','#016c59','#014636')

dlimits_aa = c(-20, -10, -5, -1, -0.5, -0.1, 0.1, 0.5, 1, 5, 10, 20)*10                  
dcols_aa = c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30')


obss = list(paste0(dir_dat, "vegfrac_refLC_refCW.nc"))
names(obss) = 'CCI'
obsLayers = list(list(1:2, 5, 3:4, 8))
years = years[1]
obsScale = 1
#score_cciLF = runComparison(mods, obss, 'figs/vegCover_aa.png', 
#                         ModLevels = list(1:5, 12:13, 6:11, 16),
#                         mod_dirs = mod_dirs, stream = stream,
#                         FUN = MM, nullFUN = null.MM,
#                         cols_aa, limits_aa, dcols_aa, dlimits_aa, TRUE)

obsLayers = list(list(1:2, c(3:5, 8)))
score_cciTR = runComparison(mods, obss, 'figs/vegCover_aa.png', 
                         ModLevels = list(1:5, c(12:13, 6:11, 16)),
                         mod_dirs = mod_dirs, stream = stream,
                         FUN = MM, nullFUN = null.MM,
                         cols_aa, limits_aa, dcols_aa, dlimits_aa, cover = TRUE)

obsLayers = list(list(c(1:2, 5), c(3:4, 8)))
score_cciWD = runComparison(mods, obss, 'figs/vegCover_aa.png', 
                         ModLevels = list(c(1:5,12:13), c(6:11, 16)),
                         mod_dirs = mod_dirs, stream = stream,
                         FUN = MM, nullFUN = null.MM,
                         cols_aa, limits_aa, dcols_aa, dlimits_aa, cover = TRUE)

obsLayers = list(list(c(3:4), c(1:2, 5, 8)))
score_cciGS = runComparison(mods, obss, 'figs/vegCover_aa.png', 
                         ModLevels = list(c(6:11), c(1:5,12:13, 16)),
                         mod_dirs = mod_dirs, stream = stream,
                         FUN = MM, nullFUN = null.MM,
                         cols_aa, limits_aa, dcols_aa, dlimits_aa, cover = TRUE)

obsLayers = list(list(c(8), c(1:5)))
score_cciBG = runComparison(mods, obss, 'figs/vegCover_aa.png', 
                         ModLevels = list(c(16), c(1:5,12:13, 6:11)),
                         mod_dirs = mod_dirs, stream = stream,
                         FUN = MM, nullFUN = null.MM,
                         cols_aa, limits_aa, dcols_aa, dlimits_aa, cover = TRUE)




dir_dat = "../fireMIPbenchmarking/data/benchmarkData/"
obss = list( VCF = paste0(dir_dat, c("treecover2000-2014.nc", "nontree2000-2014.nc")))
years = list(2001:2014)
obsLayers = list(1:168)
obsScale = 1/100
score_vcfLF = runComparison(mods, obss, 'figs/vegCover_aa.png', 
                         ModLevels = list(c(1:5, 12:13), 6:11, 16),
                         mod_dirs = mod_dirs, stream = stream,
                         FUN = MM, nullFUN = null.MM,
                         cols_aa, limits_aa, dcols_aa, dlimits_aa, TRUE)

obss = list(VCF = paste0(dir_dat, c("treecover2000-2014.nc")))
obsScale = 1/80
score_vcfWD = runComparison(mods, obss, 'figs/vegCover_aa.png', 
                         ModLevels = list(c(1:5, 12:13), c(6:11, 16)),
                         mod_dirs = mod_dirs, stream = stream,
                         FUN = MM, nullFUN = null.MM,
                         cols_aa, limits_aa, dcols_aa, dlimits_aa, TRUE)


out = mapply(function(i, j) cbind(j, rownames(i), round(i,2)),
             list(score_cciLF, score_cciTR, score_cciWD, score_cciGS, score_cciBG, 
                  score_vcfLF, score_vcfWD),
             c("CCI Life Form", "CCI Tree", "CCI Wood", "CCI Grass", "CCI Bare",
               "VCF Life", "VCF Wood"), SIMPLIFY = FALSE)
out = do.call(rbind, out)
write.csv(out, "outputs/cover_comparison_global.csv")
