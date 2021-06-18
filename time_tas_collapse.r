source("cfg.r")

jules_dir = "../jules_outputs/"
hist_mod = "u-bi607_Hist_NF_new"

Hist_yrs = 1990:2014

varname = "cv"
extent = c(-180, 180, -90, 90)

dat = openMod(hist_mod, jules_dir, varname, Hist_yrs, 1)

nn= 0
findPntsInterest <- function(y, x, bins, xbin) {
    nn<<-nn+1
    print(nn)
    out = rep(NaN, 8)
    if (sum(y) < 1) return(out)
    if (length(unique(y))<5) return(out)
    #ys = filter(y, rep(1/sbin, sbin), sides = 2)
    
    binned = mapply(function(l, u) y[x>=l & x<=u], bins[,1], bins[,2],  SIMPLIFY = FALSE)#head(bins, -1), bins[-1], SIMPLIFY = FALSE)
    #pv = mapply(function(b1, b2) wilcox.test(b1, b2, paired = F)[[3]],
    #                        head(binned, -1), binned[-1])
    
    nl = length(binned)
    mbin = sapply(binned, mean)

    lastBin = binned[[nl]]
    mxpnt = which.max(mbin)
    
    pmx = wilcox.test(binned[[mxpnt]], lastBin)[[3]]
    out[1:4] = c(xbin[mxpnt], pmx, xbin[mxpnt], xbin[mxpnt])
    if (mxpnt == nl) 
        return(out)

    for (maxMin in mxpnt:1) 
        if (wilcox.test(binned[[mxpnt]], binned[[maxMin]])[[3]] < 0.1) break

    for (maxMax in mxpnt:nl)
         if (wilcox.test(binned[[mxpnt]], binned[[maxMax]])[[3]] < 0.1) break
    
    maxMin = xbin[maxMin]; maxMax = xbin[maxMax]
    mbin = mbin[mxpnt:nl]; binned = binned[mxpnt:nl]; xbin = xbin[mxpnt:nl]
    
    mnpnt = which.min(diff(mbin))
    pmn = wilcox.test(unlist(binned[mnpnt - c(1, 0)]), unlist(tail(binned, 2)))[[3]]
    out[3:8] = c(maxMin, maxMax,  xbin[mnpnt], pmn, xbin[mnpnt], xbin[mnpnt])
    if (mnpnt == length(binned))
        return(out)

    for (minMin in mnpnt:1) 
        if (wilcox.test(binned[[mnpnt]], binned[[minMin]])[[3]] < 0.1) break

    for (minMax in mnpnt:length(xbin))
         if (wilcox.test(binned[[mnpnt]], binned[[minMax]])[[3]] < 0.1) break

    minMin = xbin[minMin]; minMax = xbin[minMax]
    out[7:8] = c(minMin, minMax)
    return(out)
}

mask = !is.na(dat[[1]])
vdat = dat[mask]

if (nlayers(dat) == length(Hist_yrs)) {
    xs = Hist_yrs
} else {
    if (nlayers(dat) == 12*length(Hist_yrs))
        xs = head(seq(Hist_yrs[1], tail(Hist_yrs, 1)+1, length.out = nlayers(dat)+1),-1)
    else browser()
}


bins = cbind(xs[1:(length(xs)-23)], xs[24:length(xs)])
xbin = apply(bins, 1, mean)
# = seq(min(xs), max(xs))
#pnts = apply(vdat, 1, findPntsInterest, xs, bins, xbin)

out = dat[[1:8]]
out[mask] = t(pnts)


