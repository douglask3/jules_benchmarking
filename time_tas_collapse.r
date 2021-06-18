source("cfg.r")

jules_dir = "/hpc/data/d05/cburton/jules_output/"
Hist_mod  = "u-bi607_HIST_NF"
Futr_mods = c(SSP1 = "u-bi607_SSP1_NF", SSP3 = "u-bi607_SSP3_NF", SSP5 = "u-bi607_SSP5_NF")
outstream = "S2.Annual"

Hist_yrs  = 1860:2014
Futr_yrs  = 2015:2100

varname = "cveg"
extent = c(-180, 180, -90, 90)

hist = openMod(hist_mod, jules_dir, varname, Hist_yrs, 1, stream = outstream)
futr = lapply(Futr_mods, openMod, jules_dir, varname, Futr_yrs, 1, stream = outstream)

dat = lapply(futr, function(i) addLayer(hist, i))
datID = sapply(Futr_mods, paste, Hist_mod, sep = '--')
years = c(Hist_yrs, Futr_yrs)

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

findMax_andMaxLoss <- function(dat, nm, id, ts, sbin = 20) {   
    tfile = paste0(c("temp/", id, nm, range(ts), ".nc"), collapse = '-')
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


dat = mapply(findMax_andMaxLoss, dat, datID, 'timeBased', MoreArgs = list(years))
