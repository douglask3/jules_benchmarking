#########
## cfg ##
#########          
graphics.off()
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("libs/")
source("../gitProjectExtras/gitBasedProjects/R/makeDir.r")
library(pracma)
################
## param list ##
################
modFacVar = 'frac'
jules_out_dir = "/hpc/data/d05/cburton/jules_output/"
outs = cbind("u-cb020_EXP4.61" = c(4.61, prior = FALSE),
             "u-cb020_EXP0.26" = c(0.26, prior = FALSE),
             "u-cb020_EXP" = c(0, prior = FALSE),
             "u-cb020_CTRL" = c(1, prior = FALSE),
             "u-cb020_EXP1.29" = c(1.29, prior = FALSE),
             "u-cb020_EXP1.20" = c(1.20, prior = FALSE),
             "u-cb020_EXP0.58" = c(0.58, prior = FALSE),
             "u-cb020_EXP0.14" = c(0.14, prior = FALSE),
             "xxp0.0" = c(0.0, prior = TRUE),
             "xxp0.0" = c(9E9, prior = TRUE))

extent = c(-180, 180, -90, 90)
pc_sample = 5

param_trans = list("BL_av_BA" = c("fun" = modParamTrans.zeroInf, "ParamsGuess" = c(-1),
                                  "funInverse" = modParamTrans.zeroInf.inverse),
                                  prior = list(prior.LogNormal, 0.4, 1)) #list(prior.Uniform, 0, 1)

variables = list("Burnt Area" = list(ObsOpenArgs   = list(obs = 'data/GFED4s_burnt_area.nc',
                                                           Layers = 1:24,
                                                           annualAver = TRUE),
                                     JulesOpenArgs = list(years = 2001:2002,
                                                           annualAver = TRUE,
                                                           fileID = 'ilamb',
                                                           varName = "burnt_area_gb",
                                                           modScale = 60*60*24*30),
                                     OptFun = logitNormal),

                 "Tree Cover" = list(ObsOpenArgs   = list(obs = 'data/TreeCover.nc',
                                                           Layers = 1,
                                                           annualAver = TRUE),
                                     JulesOpenArgs = list(years = 2001:2002,
                                                           annualAver = TRUE,
                                                           fileID = 'ilamb',
                                                           varName = "frac",
                                                           levels = list(c(1:5, 12:13)),
                                                           modScale = 1),
                                     OptFun = logitNormal))
temp_dir_move = "/data/dynamic/dkelley/jules_ES-UKESM1.1_runs"

paramDetail = 0.01
nNewParams = 30

################
## Open JULES ##
################
openJulesDat <- function(mod) {   
    openVar <- function(var, modin,...)
        dat = do.call(openMod, c(modin, jules_out_dir, var$JulesOpenArgs,
                                 temp_dir_move = temp_dir_move,...))  
    if (substr(mod, 1, 2) =='xx') {
        
        if (substr(mod, 3, 3) == "p") return(as.numeric(substr(mod, 4, nchar(mod))))
        if (substr(mod, 3, 3) == "s") {
            fill= as.numeric(substr(mod, 4, nchar(mod)))
            dat = lapply(variables, openVar, colnames(outs)[1], 
                         fill = fill, fileID2 = paste0('filled', fill))
            #browser()
        }
    } else dat = lapply(variables, openVar, mod)
    #files = list.files(dir)  
    return(dat) 
}

simss = lapply(colnames(outs), openJulesDat)

## Bootstrapo random samples
if (pc_sample < 100) {
    tfile = paste0('temp/mask-', pc_sample, filename.noPath(simss[[1]][[1]]))
    if (file.exists(tfile)) masks = brick(tfile)
    else {
        masks = !is.na(simss[[1]][[1]])
        layerMask <- function(mask) {
            ids = which(mask[])
            ids = sample(ids, size = round( length(ids) * (100-pc_sample)/100), replace = FALSE)
            mask[ids] = 0
            mask
        }
        masks = layer.apply(masks, layerMask)
        masks = writeRaster(masks, file = tfile)
    }
}

##############
## Open obs ##
##############       
openObsi <- function(var) 
    dat = do.call(openObs, c(var$ObsOpenArgs, modEG = simss[[1]][[1]]))

obss = lapply(variables, openObsi)
obss = lapply(variables, openObsi)

ps = sapply(simss, calProbs)

##############
## Open obs ##
##############  
png("figs/parameter_opt.png", height = 7, width =7, res =300, units = 'in')
if (nrow(outs) ==2) {
    optiPlot1D()
} else if (nrow(outs)==1) {
    browser()
} else {
    browser()
}
dev.off()
