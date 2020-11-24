TrendFun <- function(x) {
    
    fit = lm(x ~ t, data = data.frame(x = x, t = 1:length(x)))
    
    #res = try(summary(fit)[[4]][2, 3:4])
    res = try(c(coefficients(fit)[2], summary(fit)[[4]][2, 4]))
    
    if (is.na(res[1])) return(c(0, 1))  
    if (class(res) == "try-error") return(c(-999, 0.0))
    
    return(res)
}


makeTrendCoe <- function(dat, diff = FALSE) {
    dat0 = dat
    dat = dat/1200
    if (nlayers(dat) > 24) 
        dat = layer.apply(12:nlayers(dat), function(i) mean(dat[[(i-11):i]]))
    
    
    dat = logit(dat)
    mask = !is.na(dat[[1]])

    mdat = dat[mask]
    mtr = apply(mdat, 1, TrendFun)


    datO = dat[[1:2]]
    datO[[1]][mask] =  mtr[1,]
    datO[[2]][mask] =  mtr[2,]
    
    if (!diff) return(datO)
    
    datDT = layer.apply(1:nlayers(dat), function(t) dat[[t]] - t  * datO[[1]])

    mdat = mean(logistic(dat))
    mdatDT =  mean(logistic(datDT))
    diff =  mdat

    diff = mdat - mdatDT
    testn = mdatDT > mdat
    testp = !testn  

    diff[testn] =(mdat[testn] - mdatDT[testn])*mdatDT[testn]
    diff[testp] = (mdat[testp] - mdatDT[testp])*(1 - mdatDT[testp])


    return(2*(logit(mdat)- logit(mdatDT))/((nlayers(dat)/12)^2))
}

