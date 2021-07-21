source("cfg.r")
graphics.off()
jules_dir = "/hpc/data/d05/cburton/jules_output/"
outstream = "S2.Annual"

Hist_yrs  = 1860:2014
Futr_yrs  = 2015:2100

varnames = c("landCoverFrac", "burnt_area_gb", "cveg")
levels   = list(1:5, 1, 1)
negs = c(FALSE, TRUE, FALSE)
cols = list(c('#f7fcf5','#e5f5e0','#c7e9c0','#a1d99b','#74c476',
              '#41ab5d','#238b45','#006d2c','#00441b'), 
            c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
              '#fc4e2a','#e31a1c','#bd0026','#800026'),
            c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4',
              '#1d91c0','#225ea8','#253494','#081d58'))

limits = list(c(5, 10, 20, 40, 60, 80), c(0.1, 1, 2, 5, 10, 20, 30), c(1, 2, 4, 6, 8, 10, 15, 20))
scale = c(100, 60*60*24*365*100, 1)
extent = c(-180, 180, -90, 90)

YrCols = c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9',
         '#74add1','#4575b4','#313695')

TasCols = c('#fff7ec','#fee8c8','#fdd49e','#fdbb84','#fc8d59','#ef6548','#d7301f','#b30000','#7f0000')

tVTas = read.csv("data/time_vs_temp.csv")

pntsOfInterest = list(Cerrado = c(-52.5, -20), Catinga = c(-40, -7.5), Amozonia = c(-62.5, 0),
                      Congo = c(22.5, 0), Miombo = c(22.5, -10), 
                      "Cordilleran\nforest" = c(-125, 53.5),
                      "America\nBoreal" = c(-100, 50),
                      "Kakadu" = c(132.5, -13),
                      "Daintree" = c(145.5, -17.5),
                      "Borneo"  = c(115, 0),
                      "SE Asia\nrainforest" = c(107, 15),
                      "Siberia\nforest" = c(100, 60))

wilcoxP <- function(a, b, ...) {
    if (is.null(a) || is.null(b)) return(1)
    out = wilcox.test(a, b, ...)[[3]]
    if (is.na(out)) out = 1
    out
}

findPntsInterest <- function(y, x, bins, xbin, halfMax = FALSE, 
                             pltPoints = FALSE, col = "red", spi = 1, neg = FALSE) {
    nn <<- nn+ 1
    print(nn)
    out = rep(c(9E9, rep(NaN, 3)), 2)
    #if (sum(y) < 1) return(out)
    if (length(unique(y))<5) return(out)
    
    binned = mapply(function(l, u) y[x>=l & x<=u], bins[,1], bins[,2],  SIMPLIFY = FALSE)
    tbin = sapply(binned, length) > 1
    binned = binned[tbin]; bins = bins[tbin]; xbin = xbin[tbin]
    
    nl = length(binned)
    mbin = sapply(binned, mean)
    

    lastBin = binned[[nl]]

    if (neg) mxpnt = which.min(mbin)
    else mxpnt = which.max(mbin)
    
    pmx = wilcoxP(binned[[mxpnt]], lastBin)
    out[1:4] = c(xbin[mxpnt], pmx, xbin[mxpnt], xbin[mxpnt])
    
    if (pltPoints) {        
        prec = sapply(binned, quantile, c(0.1, 0.25, 0.75, 0.9))

        pointPol <- function(index, bindex,  col) {
            xh = x[index]; yh = y[index]; xbinh = xbin[bindex]; ph = prec[,bindex]
            points(xh, yh, pch = 4, col = col)
            for (i in 1:2) polygon(c(xbinh, rev(xbinh)), c(ph[i,], rev(ph[5-i,])),
                                   col = make.transparent(col, 0.8), border = NA)
            
        }
        pointPol(1:length(Hist_yrs), xbin<x[length(Hist_yrs)], make.transparent('black', 0.67))
       
        pointPol(length(Hist_yrs):length(x), xbin>=x[length(Hist_yrs)], col)
        #pointPol(1:length(x), 1:length(xbin), col)            
    }
   
    if (mxpnt == nl) 
        return(out)
    
    for (maxMin in mxpnt:1) 
        if (wilcoxP(binned[[mxpnt]], binned[[maxMin]]) < 0.1) break

    for (maxMax in mxpnt:nl)
         if (wilcoxP(binned[[mxpnt]], binned[[maxMax]]) < 0.1) break
    xbin0 = xbin; mbin0 = mbin
    maxMin = xbin[maxMin]; maxMax = xbin[maxMax]
    mbin = mbin[mxpnt:nl]; binned = binned[mxpnt:nl]; xbin = xbin[mxpnt:nl]
    
    #xbin<x[length(Hist_yrs)]
    if(neg) FUN = which.max else FUN = which.min
    if (halfMax) {
        index = FUN(abs(xbin-x[length(Hist_yrs)]))
        if (neg)  mnpnt= max(which(mbin > (2*mbin[index])))
        else mnpnt= max(which(mbin > (mbin[index]/2)))
    } else mnpnt = FUN(diff(mbin))
    
    pmn = wilcoxP(unlist(binned[mnpnt - c(1, 0)]), unlist(tail(binned, 2)))#[[3]]#
    pmn = 1-(1-pmx) * (1-pmn)
   # print(xbin[mnpnt])
    if (is.na(xbin[mnpnt])) {
        if (mnpnt < 0) mnpnt = 1
        else browser()
    }
    out[3:8] = c(maxMin, maxMax,  xbin[mnpnt], pmn, xbin[mnpnt], xbin[mnpnt])

    if (pltPoints) {
        xp = c(maxMin, maxMax); yp = rep(par("usr")[4] * spi * 0.02, 2)#rep(max(y) + 0.02*diff(range(y)), 2)
        lines(xp, yp, col = col)
        points(xp, yp, pch = 19, col = col, cex = 0.5)
    }
    
    if (mnpnt == length(binned))
        return(out)

    for (minMin in mnpnt:1) 
        if (wilcoxP(binned[[mnpnt]], binned[[minMin]]) < 0.1) break

    for (minMax in mnpnt:length(xbin))
         if (wilcoxP(binned[[mnpnt]], binned[[minMax]]) < 0.1) break

    minMin = xbin[minMin]; minMax = xbin[minMax]
    out[7:8] = c(minMin, minMax)
    
    if (pltPoints) {
        xp = c(minMin, minMax); yp = par("usr")[4] - yp#; yp = mbin[mnpnt] + diff(mbin)[mnpnt] *(xp - mean(xp)) + 0.1*diff(range(y))
        lines(xp, yp, col = col)
        points(xp, yp, pch = 19, col = col, cex = 0.5)
    }
    
    return(out)
}
nn <<- 0
findMax_andMaxLoss <- function(dat, nm, id, ts, ..., halfMax = FALSE, neg = FALSE) {   
    tfile = paste0(c("temp/", id, nm, range(ts), c('', 'neg')[neg+1],
                   halfMax, ".nc"), collapse = '-')
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
    #if (grepl("tempBased", id)) {
    #    bins = seq(floor(min(ts)/sbin[1])*sbin[1], ceiling(max(ts)/sbin[1])*sbin[1], sbin[1])
    #    #browser()
    #    bins = cbind(head(bins, -sbin[2]/sbin[1]), tail(bins, -sbin[2]/sbin[1]))#
    #} else 
    #    bins = cbind(xs[1:(length(xs)-round(sbin/2)-1)], xs[round(sbin/2):length(xs)])
    
    
    # = seq(min(xs), max(xs))
    nn <<- 0
    pnts = apply(vdat, 1, findPntsInterest, xs, halfMax = halfMax, neg = neg,...)

    out = dat[[1:8]]
    out[mask] = t(pnts)
    names(out) = c("MaxPnt", "MaxPval", "MaxRangeMin", "MaxRangeMax",
                   "LossPnt", "LossPval", "LossRangeMin", "LossRangeMax")
    out = writeRaster(out, file = tfile, overwrite = TRUE)
    
    return(out)
}

plotFUN <- function(dat1, dat2, v1, v2) {
    #ncol = length(dat)
    #nrow = 4
    #lmat = lmat0 = matrix(1:(nrow*ncol), ncol = ncol)
    #for (i in 1:(nrow-1)) 
    #    lmat = rbind(lmat[1:(2*i-1),], nrow*ncol+i, lmat[((2*i)):nrow,])
    lmat = rbind(1:3, 4:6, 7:9, 10, 11:13, 14:16, 17:19, 20)
    
    layout(lmat)
    par(mar = rep(0, 4), oma = c(0, 3, 2, 0))
    
    #cols = c('#a50026','#f46d43','yellow','#74add1','#313695')
    limits = 2015:2099
    labels = limits
    labels[] = ''

    labels[seq(1, length(limits), by = 5)] = limits[seq(1, length(limits), by = 5)]
    plotStandardMap.maskMax <- function(x, ..., mv = NULL, Ttxt = '', Stxt = '') {
        if (is.null(mv)) mv =  max.raster(x, na.rm = TRUE)
        overtime = x == mv
        x[overtime] = NaN
        plotStandardMap(x, ..., ePatternRes = 30,  ePatternThick = 0.6)
        try(plotStandardMap(overtime, limits = c(0.5), col = c("transparent", "black"),
                             add = TRUE))
        mtext(side = 3, line = 1, Ttxt)
        mtext(side = 2, line = -1, Stxt)
        grid() 
        #xy = xyFromCell(overtime, which(overtime[]))
        
    }
    pNorm <- function(dat, j, expN = '') {
        if (j == 1) nms = names(Futr_mods) else nms = ''
        exp = c(expN, "", "")
        mapply(function(i, txt1, txt2) plotStandardMap.maskMax(i[[j]], '', limits, 
                                    cols = YrCols,
                                    e = 1-i[[j+1]], limits_error = c(0.05, 0.1), 
                                    Ttxt = txt1, Stxt = txt2),
               dat, nms, exp)
    }#lNorm <- function(j) 
    pNorm(dat1, 1 ,"Year of\nmax.")    
    pNorm(dat1, 5, "Year of\nmax. loss")    
        
    pNorm(dat2, 5, "Year of\n50% loss")
    StandardLegend(YrCols, limits, dat1[[1]][[5]], oneSideLabels = T, extend_min = TRUE, 
                   labelss = labels)
    limits = c(0, 0.5, 1, 1.5, 2, 3, 4, 6, 8)
    #cols = rev(cols)
    
    pNorm <- function(dat, j, expN) {
        #if (j == 1) nms = names(Futr_mods) else nms = ''
        exp = c(expN, "", "")
        FUN <- function(i, v2i, txt2) {
            r = i[[j]]
            r[] = v2i[round(i[[j]])[]-1860]
            
            plotStandardMap.maskMax(r, '', limits, cols = TasCols,
                                    e = 1-i[[j+1]], limits_error = c(0.05, 0.1),
                                    mv = v2i[round(max.raster(i[[j]], na.rm = TRUE)-1860)],
                                    Ttxt = '', Stxt = txt2)
        }
        mapply(FUN, dat, v2, exp)
    }
    
    pNorm(dat1, 1, "Temp of\nmax.")
    pNorm(dat2, 5, "Temp of\nmax. loss")
    pNorm(dat2, 5, "Temp of\n50% loss")
    StandardLegend(TasCols, limits, dat1[[1]][[5]], oneSideLabels = T)  
}

openRunExp <- function(hist_mod, Futr_mods, varname, levels = 1,
                       cols = NULL, limits = NULL, scale, neg, fname) {
    #if (varname == varnames[1]) return()
    hist = openMod(hist_mod, jules_dir, varname, Hist_yrs, 1, 
                   levels = levels, stream = outstream)
    futr = lapply(Futr_mods, openMod, jules_dir, varname, Futr_yrs, 1, 
                  levels = levels, stream = outstream)
    
    Hist_mod <<- hist_mod
    Futr_mods <<- Futr_mods
     
    if (!is.null(cols) && !is.null(limits)) {
        png(paste0("figs/historic_map-", varname, '-', fname, ".png"), 
            height = 5, width = 7, res = 300, units = "in")
            plotStandardMap(mean(hist[[nlayers(hist) + -10:0]])*scale, cols = cols,
                            limits = limits)
            lapply(pntsOfInterest, function(i) points(i[1], i[2], pch = 4, cex = 1.5, lwd = 1.5))
            StandardLegend(cols, limits, hist[[1]], add = TRUE, oneSide = FALSE)
        dev.off()
    }
    #browser()
    dat = lapply(futr, function(i) addLayer(hist, i))
    
    plotPoints <- function(pntsOI, name, x = years, bins, xbin, halfMax) {
        pnt = cellFromXY(dat[[1]], pntsOI)
        tfile = paste0("temp/", varname, '-', paste0(datID, collapse = '_'), '-', 
                       pnt, 'halfMax', '.Rd')
        if (file.exists(tfile)) load(tfile)
        else {
            dats = lapply(dat, function(i) i[pnt])
            save(dats, file = tfile)
        }
        
        plot(range(x), c(0, 1.1*max(unlist(dats))), type = 'n')
        polygon(par("usr")[c(1, 2, 2, 1, 1)], par("usr")[c(3, 3, 4, 4, 3)], col = "white")
        grid()
        mtext(side = 2, adj = 0.1, padj = 1, line = -1.5, name)
        if (!is.list(x)) x = list(x, x, x)
        #browser()
        mapply(findPntsInterest, dats, x, col = c('#7570b3', '#1b9e77','#d95f02'), spi = 1:3, 
               MoreArgs = list(bins, xbin, halfMax = halfMax, neg = neg, pltPoints = TRUE))
        #apply(dats, 2, findPntsInterest, years, bins, xbin, TRUE)
    }    
    
    
    datID = sapply(Futr_mods, paste, hist_mod, sep = '--')

    years = c(Hist_yrs, Futr_yrs)
   

    tVTas = tVTas[, c(1, 1+which(grepl('NF', colnames(tVTas)[-1]) == grepl("NF", hist_mod)))]
    tass = sapply(years, function(yr) tVTas[yr == tVTas[,1],])
    tass = lapply(as.list(data.frame(t(tass))), unlist)

    nrow = ceiling(sqrt(length(pntsOfInterest))); ncol = length(pntsOfInterest)/nrow
    TSplot <- function(halfMax, name, xlab, xs, bins, xbin) {
        png(paste0("figs/site_", name, '-', varname, '-', fname, '-', halfMax, ".png"),
            height = 9*2, width = 8*1.5, units = 'in', res = 300)
            par(mfrow = c(nrow, ncol), mar = c(0.1, 2.5, 0.1, 0.1), oma = c(4, 2, 0, 0)) 
            mapply(plotPoints, pntsOfInterest, names(pntsOfInterest),
                   MoreArgs = list(xs, bins, xbin, halfMax = halfMax))
            mtext.units(side = 1, xlab, line = 2, outer = TRUE)
            mtext.units(side = 2, 'kg C~m-2~', line = -1, outer = TRUE)
        dev.off()
    }
    sbin = c(0.1, 1)
    bins = seq(floor(min(unlist(tass[-1]))/sbin[1])*sbin[1], 
               ceiling(max(unlist(tass[-1]))/sbin[1])*sbin[1], sbin[1])
    bins = cbind(head(bins, -sbin[2]/sbin[1]), tail(bins, -sbin[2]/sbin[1]))
    #browser()
    xbin = apply(bins, 1, mean)
    lapply(c(T, F), TSplot,"tasSeries", '~DEG~C', tass[-1], bins, xbin)
    #datTas  = mapply(findMax_andMaxLoss, dat, datID, paste(varname, '-tempBased'), tass[-1], 
    #             MoreArgs = list(c(0.1, 1), neg = neg))

        
    xs = years; sbin = 21
    bins = cbind(xs[1:(length(xs)-round(sbin)+1)], xs[round(sbin):length(xs)]) 
    xbin = apply(bins, 1, mean)
    lapply(c(T, F), TSplot,"timeSeries", 'year', years, bins, xbin)

    datTime1 = mapply(findMax_andMaxLoss, dat, datID, paste(varname, '-timeBased'),
                 MoreArgs = list(years, bins, xbin, neg = neg))
    datTime2 = mapply(findMax_andMaxLoss, dat, datID, paste(varname, '-timeBased'),
                 MoreArgs = list(years, bins, xbin, halfMax = TRUE, neg = neg))
    

    png(paste0("figs/max_loss", '-', varname, '-', fname, ".png"), 
        height = 10,width = 8, units = 'in', res = 300)
    plotFUN(datTime1, datTime2, years, tass[-1])
    dev.off()
}


runVarname <- function(varname, levels, cols, limits, scale, neg) {
    Hist_mod  = "u-bi607_HIST"
    Futr_mods = c(SSP1 = "u-cd136_SSP1", SSP3 = "u-cd136_SSP3", SSP5 = "u-cd136_SSP5")
    dats_WF = openRunExp(Hist_mod, Futr_mods, varname, levels, cols, limits, scale, neg, 
                         "withFire")

    if (grepl("burnt_area", varname)) return()
    Hist_mod  = "u-bi607_HIST_NF"
    Futr_mods = c(SSP1 = "u-bi607_SSP1_NF", SSP3 = "u-bi607_SSP3_NF", SSP5 = "u-bi607_SSP5_NF")
    dats_NF = openRunExp(Hist_mod, Futr_mods, varname, levels, cols, limits, scale, neg, 
                         "noFire")
}

mapply(runVarname, varnames, levels, cols, limits, scale, negs)
