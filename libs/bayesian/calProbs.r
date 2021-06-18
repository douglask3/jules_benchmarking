calProbs <- function(sims) {
    calProb <- function(obs, sim, var) {
        if (is.numeric(sim)) {
            if (sim == 0) return(-1000000)
            else return(log(sim))
        }
        nl = min(nlayers(obs), nlayers(sim))
        grabV <- function(i, r) {
            mask = masks[[i]]
            out = r[[i]][mask]
            
            return(out)
        }
        grabVs <- function(r) {
            vName = paste0('temp/', filename.noPath(r, TRUE), maskName, '.csv')
            
            if (file.exists(vName) & F) vs = read.csv(vName)[,1]
            else {
                vs = lapply(1:nl, grabV, r)
                vs = unlist(vs)
                write.csv(vs, file = vName, col.names = FALSE, row.names = FALSE)               
            }
            return(vs)
        }
        maskName = filename.noPath(masks, TRUE)
        vObs = vObs0 =  grabVs(obs)        
        vSim = vSim0 =  grabVs(sim)
        
        tst = !(is.na(vObs) | is.na(vSim))
        vObs = vObs[tst]; vSim = vSim[tst]
        p = var$OptFun(vObs, vSim)
    }
    ps = mapply(calProb, obss, sims, variables)  
    #ps[1] = 0 
    return(mean(ps))
}

