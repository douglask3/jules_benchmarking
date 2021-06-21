source("cfg.r")

jules_dir = "/hpc/data/d05/cburton/jules_output/"
outstream = "S2.Annual"

Hist_yrs  = 1860:2014
Futr_yrs  = 2015:2100

varname = "cveg"
extent = c(-180, 180, -90, 90)

tVTas = read.csv("data/time_vs_temp.csv")

wilcoxP <- function(...) {
    out = wilcox.test(...)[[3]]
    if (is.na(out)) out = 1
    out
}

findPntsInterest <- function(y, x, bins, xbin) {
    nn <<- nn+ 1
    print(nn)
    out = rep(NaN, 8)
    #if (sum(y) < 1) return(out)
    if (length(unique(y))<5) return(out)
    
    binned = mapply(function(l, u) y[x>=l & x<=u], bins[,1], bins[,2],  SIMPLIFY = FALSE)
    nl = length(binned)
    mbin = sapply(binned, mean)

    lastBin = binned[[nl]]
    mxpnt = which.max(mbin)
    
    pmx = wilcoxP(binned[[mxpnt]], lastBin)
    out[1:4] = c(xbin[mxpnt], pmx, xbin[mxpnt], xbin[mxpnt])
    if (mxpnt == nl) 
        return(out)

    
    for (maxMin in mxpnt:1) 
        if (wilcoxP(binned[[mxpnt]], binned[[maxMin]]) < 0.1) break

    for (maxMax in mxpnt:nl)
         if (wilcoxP(binned[[mxpnt]], binned[[maxMax]]) < 0.1) break
    xbin0 = xbin; mbin0 = mbin
    maxMin = xbin[maxMin]; maxMax = xbin[maxMax]
    mbin = mbin[mxpnt:nl]; binned = binned[mxpnt:nl]; xbin = xbin[mxpnt:nl]
    
    mnpnt = which.min(diff(mbin))
    pmn = wilcox.test(unlist(binned[mnpnt - c(1, 0)]), unlist(tail(binned, 2)))[[3]]#
    pmn = 1-(1-pmx) * (1-pmn)
    out[3:8] = c(maxMin, maxMax,  xbin[mnpnt], pmn, xbin[mnpnt], xbin[mnpnt])
    if (mnpnt == length(binned))
        return(out)

    for (minMin in mnpnt:1) 
        if (wilcoxP(binned[[mnpnt]], binned[[minMin]]) < 0.1) break

    for (minMax in mnpnt:length(xbin))
         if (wilcoxP(binned[[mnpnt]], binned[[minMax]]) < 0.1) break

    minMin = xbin[minMin]; minMax = xbin[minMax]
    out[7:8] = c(minMin, minMax)
    
    return(out)
}

findMax_andMaxLoss <- function(dat, nm, id, ts, sbin = 20, neg = FALSE) {   
    tfile = paste0(c("temp/", id, nm, range(ts), neg, ".nc"), collapse = '-')
    print(tfile)
    if (file.exists(tfile)) return(brick(tfile))
    mask = !is.na(dat[[1]])
    vdat = dat[mask]
    
    if (nlayers(dat) == length(ts)) {
        xs = ts
    } else {
        if (nlayers(dat) == 12*length(ts))
            xs = head(seq(ts[1], tail(ts, 1)+1, length.out = nlayers(dat)+1),-1)
        else browser()
    }
    #browser()
    if (grepl("tempBased", id)) {
        bins = seq(floor(min(ts)/sbin[1])*sbin[1], ceiling(max(ts)/sbin[1])*sbin[1], sbin[1])
        #browser()
        bins = cbind(head(bins, -sbin[2]/sbin[1]), tail(bins, -sbin[2]/sbin[1]))#
    } else 
        bins = cbind(xs[1:(length(xs)-round(sbin/2)-1)], xs[round(sbin/2):length(xs)])
    
    xbin = apply(bins, 1, mean)
    # = seq(min(xs), max(xs))
    nn <<- 0
    pnts = apply(vdat, 1, findPntsInterest, xs, bins, xbin)

    out = dat[[1:8]]
    out[mask] = t(pnts)
    names(out) = c("MaxPnt", "MaxPval", "MaxRangeMin", "MaxRangeMax",
                   "LossPnt", "LossPval", "LossRangeMin", "LossRangeMax")
    out = writeRaster(out, file = tfile, overwrite = TRUE)
    return(out)
}

plotFUN <- function(dat, v1, v2) {
    #ncol = length(dat)
    #nrow = 4
    #lmat = lmat0 = matrix(1:(nrow*ncol), ncol = ncol)
    #for (i in 1:(nrow-1)) 
    #    lmat = rbind(lmat[1:(2*i-1),], nrow*ncol+i, lmat[((2*i)):nrow,])
    lmat = rbind(1:3, 4:6, 7, 8:10, 11:13, 14)
    
    layout(lmat)
    par(mar = rep(0, 4), oma = c(0, 2, 2, 0))
    cols = c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695')
    #cols = c('#a50026','#f46d43','yellow','#74add1','#313695')
    limits = 2015:2099
    labels = limits
    labels[] = ''
    labels[seq(1, length(limits), by = 5)] = limits[seq(1, length(limits), by = 5)]
    plotStandardMap.maskMax <- function(x, ..., mv = NULL, Ttxt = '', Stxt = '') {
        if (is.null(mv)) mv =  max.raster(x, na.rm = TRUE)
        overtime = x == mv
        x[overtime] = NaN
        plotStandardMap(x, ...)
        try(plotStandardMap(overtime, limits = c(0.5), col = c("transparent", "black"),
                             add = TRUE))
        mtext(side = 3, line = 1, Ttxt)
        mtext(side = 2, line = -1, Stxt)
        
        #xy = xyFromCell(overtime, which(overtime[]))
        
    }
    pNorm <- function(j, expN = '') {
        if (j == 1) nms = names(Futr_mods) else nms = ''
        exp = c(expN, "", "")
        mapply(function(i, txt1, txt2) plotStandardMap.maskMax(i[[j]], '', limits, cols = cols,
                                    e = 1-i[[j+1]], limits_error = c(0.05, 0.1), 
                                    Ttxt = txt1, Stxt = txt2),
            dat, nms, exp)
    }#lNorm <- function(j) 
    pNorm(1, "Year max.")    
    pNorm(5, "Year Loss")    
    StandardLegend(cols, limits, dat[[1]][[5]], oneSideLabels = T, extend_min = TRUE, labelss = labels)
    limits = c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8)
    cols = rev(cols)
    
    pNorm <- function(j, expN) {
        #if (j == 1) nms = names(Futr_mods) else nms = ''
        exp = c(expN, "", "")
        FUN <- function(i, v2i, txt2) {
            r = i[[j]]
            r[] = v2i[round(i[[j]])[]-1860]
            
            plotStandardMap.maskMax(r, '', limits, cols = cols,
                                    e = 1-i[[j+1]], limits_error = c(0.05, 0.1),
                                    mv = v2i[round(max.raster(i[[j]], na.rm = TRUE)-1860)],
                                    Ttxt = '', Stxt = txt2)
        }
        mapply(FUN, dat, v2, exp)
    }

    pNorm(1, "Temp max.")
    pNorm(5, "Temp Loss")
    StandardLegend(cols, limits, dat[[1]][[5]], oneSideLabels = T)
   
}

openRunExp <- function(hist_mod, Futr_mods, varname, fname) {
    hist = openMod(hist_mod, jules_dir, varname, Hist_yrs, 1, stream = outstream)
    futr = lapply(Futr_mods, openMod, jules_dir, varname, Futr_yrs, 1, stream = outstream)

    dat = lapply(futr, function(i) addLayer(hist, i))
    datID = sapply(Futr_mods, paste, Hist_mod, sep = '--')

    years = c(Hist_yrs, Futr_yrs)

    tVTas = tVTas[, c(1, 1+which(grepl('NF', colnames(tVTas)[-1]) == grepl("NF", hist_mod)))]
    tass = sapply(years, function(yr) tVTas[yr == tVTas[,1],])
    tass = lapply(as.list(data.frame(t(tass))), unlist)

    #datTas  = mapply(findMax_andMaxLoss, dat, datID, paste(varname, '-tempBased'), tass[-1], 
    #             MoreArgs = list(c(0.1, 1)))
    
    datTime = mapply(findMax_andMaxLoss, dat, datID, paste(varname, '-timeBased'),
                 MoreArgs = list(years))

    png(paste0("figs/max_loss", fname, ".png"), height = 7.2,width = 8, units = 'in', res = 300)
    plotFUN(datTime, years, tass[-1])
    dev.off()
}


Hist_mod  = "u-bi607_HIST"
Futr_mods = c(SSP1 = "u-cd136_SSP1", SSP3 = "u-cd136_SSP3", SSP5 = "u-cd136_SSP5")
dats_WF = openRunExp(Hist_mod, Futr_mods, varname, "withFire")

Hist_mod  = "u-bi607_HIST_NF"
Futr_mods = c(SSP1 = "u-bi607_SSP1_NF", SSP3 = "u-bi607_SSP3_NF", SSP5 = "u-bi607_SSP5_NF")
dats_NF = openRunExp(Hist_mod, Futr_mods, varname, "noFire")
