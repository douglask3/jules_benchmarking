graphics.off()
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
sourceAllLibs("libs/")

jules_out_dir = "/hpc/data/d05/cburton/jules_output/"
outs = cbind("u-cb020_CTRL" = c(1, prior = FALSE),
             "u-cb020_EXP" = c(0, prior = FALSE),
             "u-cb020_EXP0.26" = c(0.26, prior = FALSE),
             "u-cb020_EXP1.29" = c(1.29, prior = FALSE),
             "xxs0.99999" = c(9E9, prior = TRUE))


extent = c(-180, 180, -90, 90)
pc_sample = 5

param_trans = list("BL_av_BA" = c("fun" = modParamTrans.zeroInf, "ParamsGuess" = c(-1),
                                  "funInverse" = modParamTrans.zeroInf.inverse),
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
    xf = x#logit(x)
    params0 = params
    params = param_trans[[1]]$fun(outs[1,])
    
    test = is.na(outs[2,])
    if (ncol(outs) <=3) outs[2,test] = TRUE else {
        outs = outs[, !test]
        ps = ps[!test]
        params = params[!test]
    }
    #browser()
    test = outs[2,]==1
    nps = length(ps)
    
    paramsF = params#logit(params)
    if (any(test)) {
        fit = polyfix(paramsF[test], ps[test], nps , paramsF[!test], ps[!test])
    } else fit = polyfit(paramsF, ps, nps)
    y = polyval(fit, xf) * do.call(param_trans$prior[[1]], list(x, param_trans$prior[-1]))
    xrange = range(c(params, x)); yrange = range(c(ps, y))
    par(mfrow = c(2, 1))
    plot(xrange, yrange, type = 'n', xaxt = 'n', xlab = '', ylab = '')
    mtext.units(side = 2, line = 2, 'log(P(~beta~|Y)')
    mtext.units(side = 1, line = 2, 'log(P(~beta~)')
    
    labs = signif(param_trans[[1]]$funInverse(seq(0, 1, length.out = 7)),1)
    at = param_trans[[1]]$fun(labs)
    axis(side = 1, at = at, label = labs)

    points(params[test], ps[test], pch = 19, col = 'blue')
    points(params[!test], ps[!test], pch = 19, col = 'red')
    lines(x, y)

    dy = diff(y)
    turning = which(dy[-1] <0 & head(dy, -1) >0) + 1
    yturn = y[turning]; xturn = x[turning]
    turning = sort.int(yturn, decreasing = 2, index.return = TRUE)[[2]]
    yturn = yturn[turning] ;  xturn = xturn[turning]

    params = params[!test]; ps = ps[!test]; target = max(ps)
    checkLikelyMissing <- function(addinX, addinY) {
        if (!is.na(addinX)) {
            params = c(params, addinX);
            ps = c(ps, addinY)
        }
        index = sort.int(params, index.return = TRUE)[[2]]
        params = params[index]
        ps = ps[index]
        
        params = c(0, params, 1)
        ps = c(ps[1], ps, tail(ps, 1))
        dif = target - ps

        midway <- function(i) {
            #if (i == 4) browser()
            (params[i-1] * dif[i] + params[i] * dif[i-1])/(dif[i] + dif[i-1])
        }
        mids = sapply(2:length(ps), midway)
        dif = mapply(mean, dif[-1], head(dif, -1))
        test = dif >0
        dif = dif[test]
        liki = which.min(dif/diff(params)[test])
        
        return(mids[liki])   
    }
    likliMis = c()
    for (i in 1:length(xturn)) {
         
        likliMis = c(likliMis, checkLikelyMissing(c(xturn[1:i], likiMis), c(yturn[1:i], rep(target, length(likiMis)))))
       
    }
    
    likliMis = unique(likliMis)
    sapply(xturn, function(i) lines(c(i, i), c(-9E9, 9E9), lty = 2))
    sapply(likliMis, function(i) lines(c(i, i), c(-9E9, 9E9), lty = 3))
    xturn = round(param_trans[[1]]$funInverse(xturn), 2)
    
    likliMis = round(param_trans[[1]]$funInverse(likliMis), 2)
    
    plot(c(0, 1), c(0, 1), type = 'n', axes = FALSE, xlab = '', ylab = '')
    text.units(x = 0.01, y = 0.9, adj = 0,
         paste0("Best performance: ", outs[1,which.max(ps)],
                "; log(P(~beta~|Y)) = ", round(max(ps))))
    text(x = 0.01, y = 0.7, adj = 0, "likley max(P) in order of preference:")
    
    tab1 = cbind(xturn, round(yturn))
    
    text(x = 0.01, y = 0.6, adj = 0, "Param. value")
    text.units(x = 0.26, y = 0.6, adj = 0, "log(P(~beta~|Y))")
    
    textTfun <- function(i, y, ...)
        text(i, adj = 0, y = y - seq(0, by = 0.1, length.out = length(i)),...)
    textTfun(tab1[,1], x = 0.01, y = 0.5)
    textTfun(tab1[,2], x = 0.26, y = 0.5)
    
    text(x = 0.61, y = 0.7, adj = 0, "likley false optimization points")
    text(x = 0.61, y = 0.6, adj = 0, "in order of preference:")
    textTfun(likliMis, x = 0.61, y = 0.5)
} else if (nrow(out)==1) {
    browser()
} else {
    browser()
}

