library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("libs/")

jules_out_dir = "/hpc/data/d05/cburton/jules_output/"
outs = cbind("u-cb020_CTRL" = c(1, prior = FALSE),
             "u-cb020_EXP" = c(0, prior = FALSE),
             "xxs0.99999" = c(9E9, prior = TRUE))


extent = c(-180, 180, -90, 90)
pc_sample = 5

param_trans = list("BL_av_BA" = c("fun" = modParamTrans.zeroInf, "ParamsGuess" = c(-1)),
                                  prior = list(prior.Uniform, 0, 1))

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
    openVar <- function(var, modin,...)
        dat = do.call(openMod, c(modin, jules_out_dir, var$JulesOpenArgs, ...))  
    if (substr(mod, 1, 2) =='xx') {
        if (substr(mod, 3, 3) == "p") return(as.numeric(substr(mod, 3, nchar(mod))))
        if (substr(mod, 3, 3) == "s") {
            fill= as.numeric(substr(mod, 4, nchar(mod)))
            dat = lapply(variables, openVar, colnames(outs)[1], 
                         fill = fill, fileID2 = paste0('filled', fill))
            #browser()
        }
    } else dat = lapply(variables, openVar, mod)
    #files = list.files(dir)
      
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
        if (is.numeric(sim)) {
            if (sim == 0) return(-1000000)
            else return(log(sim))
        }
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

if (nrow(outs) ==2) {
    x = seq(0, 1, 0.01)
    params = param_trans[[1]]$fun(outs[1,])
    

    test = outs[2,]==1
    nps = length(ps)
    if (any(test)) {
        fit = polyfix(params[test], ps[test], nps , params[!test], ps[!test])
    } else fit = polyfit(params, ps, nps)
    y = polyval(fit, x) * do.call(param_trans$prior[[1]], list(x, param_trans$prior[-1]))
    xrange = range(c(params, x)); yrange = range(c(ps, y))
    plot(xrange, yrange, type = 'n')
    points(params[test], ps[test], pch = 19, col = 'blue')
    points(params[!test], ps[!test], pch = 19, col = 'red')
    lines(x, y)
    browser()
} else if (nrow(out)==1) {
    browser()
} else {
    browser()
}

