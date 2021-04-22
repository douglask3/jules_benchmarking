library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("libs/")

jules_out_dir = "/hpc/data/d05/cburton/jules_output/"
outs = cbind("u-cb020_CTRL" = c(1),
             "u-cb020_EXP" = c(0))
extent = c(-180, 180, -90, 90)


zeroInf <- function(x, a = -1) 1-exp(a*x)

logitNormal <- function(obs, sim) {
    obs0 = obs == 0
    obs1 = !obs0
    obs = logit(obs)
    simt = logit(sim)
    sim2 = (1-sim)^2
    FUN <- function(params) {
        P0 = sim2 * (1 - params[1])
        P1 =  (1 - P0[obs1]) * (dnorm(obs[obs1], simt[obs1], params[2]))
        
        -sum(log(P0[obs0])) - sum(log(P1))
    }
    -optim(c(0.1, 1), FUN)$value
}

param_trans = list("BL_av_BA" = c("fun" = zeroInf, "ParamsGuess" = c(-1)))

variables = list("Burnt Area" = list(ObsOpenArgs   = list(obs = 'data/GFED4s_burnt_area.nc',
                                                           Layers = 1:12,
                                                           annualAver = TRUE),
                                      JulesOpenArgs = list(years = 2001:2002,
                                                           annualAver = TRUE,
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

openObsi <- function(var)
    dat = do.call(openObs, c(var$ObsOpenArgs, modEG = sim[[1]][[1]]))

obss = lapply(variables, openObsi)

calProbs <- function(sims) {
    calProb <- function(obs, sim, var) {
        mask = !is.na(obs + sim)
        p = var$OptFun(obs[mask], sim[mask])
    }
    ps = mapply(calProb, obss, sims, variables)   
}

ps = sapply(simss, calProbs)

