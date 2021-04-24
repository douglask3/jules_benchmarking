library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("libs/")

jules_out_dir = "/hpc/data/d05/cburton/jules_output/"
outs = cbind("u-cb020_CTRL" = c(1),
             "u-cb020_EXP" = c(0))


extent = c(-180, 180, -90, 90)
pc_sample = 5

param_trans = list("BL_av_BA" = c("fun" = modParamTrans.zeroInf, "ParamsGuess" = c(-1)))

variables = list("Burnt Area" = list(ObsOpenArgs   = list(obs = 'data/GFED4s_burnt_area.nc',
                                                           Layers = 1:24,
                                                           annualAver = FALSE),
                                     JulesOpenArgs = list(years = 2001:2002,
                                                           annualAver = FALSE,
                                                           fileID = 'ilamb',
                                                           varName = "burnt_area_gb",
                                                           modScale = 60*60*24*30),
                                     OptFun = logitNormal))

openJulesDat <- function(mod) {    
    #files = list.files(dir)
    openVar <- function(var)
        dat = do.call(openMod, c(mod, jules_out_dir, var$JulesOpenArgs))    
    dat = lapply(variables, openVar)
}

simss = lapply(colnames(outs), openJulesDat)
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
        
openObsi <- function(var) 
    dat = do.call(openObs, c(var$ObsOpenArgs, modEG = simss[[1]][[1]]))

obss = lapply(variables, openObsi)
obss = lapply(variables, openObsi)


calProbs <- function(sims) {
    calProb <- function(obs, sim, var) {
        nl = min(nlayers(obs), nlayers(sim))
        grabV <- function(i, r) {
            mask = !is.na(r[[i]])& masks[[i]]
            out = r[[i]][mask]
            return(out)
        }
        grabVs <- function(r) {
            vName = paste0('temp/', filename.noPath(r, TRUE), maskName, '.csv')
            
            if (file.exists(vName)) vs = read.csv(vName)[,1]
            else {
                vs = lapply(1:nl, grabV, r)
                vs = unlist(vs)
                write.csv(vs, file = vName, col.names = FALSE, row.names = FALSE)               
            }
            return(vs)
        }
        maskName = filename.noPath(masks, TRUE)
        vObs = grabVs(obs)        
        vSim = grabVs(sim)
        p = var$OptFun(vObs, vSim)
    }
    ps = mapply(calProb, obss, sims, variables)   
}

ps = sapply(simss, calProbs)

