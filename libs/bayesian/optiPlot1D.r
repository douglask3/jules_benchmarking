optiPlot1D <- function() {
    x = seq(0, 1, paramDetail)
    xf = x#logit(x)
    params = param_trans[[1]]$fun(outs[1,])
    
    PriorTest = is.na(outs[2,])
    if (ncol(outs) <=3) outs[2,PriorTest] = TRUE else {
        outs = outs[, !PriorTest]
        ps = ps[!PriorTest]
        params = params[!PriorTest]
    }
    #browser()
    PriorTest = outs[2,]==1
    nps = length(ps)
    
    paramsF = params#logit(params)
    if (any(PriorTest)) {
        fit = polyfix(paramsF[PriorTest], ps[PriorTest], nps ,
                      paramsF[!PriorTest], ps[!PriorTest])
    } else fit = polyfit(paramsF, ps, nps)
    y = polyval(fit, xf) * do.call(param_trans$prior[[1]], list(x, param_trans$prior[-1]))
    xrange = range(c(params, x)); yrange = range(c(ps, y))
    par(mfrow = c(2, 1))
    plot(xrange, yrange, type = 'n', xaxt = 'n', xlab = '', ylab = '')

    mtext.units(side = 2, line = 2, 'log(P(~beta~|Y) P(Y))')
    mtext.units(side = 1, line = 2, '~beta~')
    
    labs = signif(param_trans[[1]]$funInverse(seq(0, 1, length.out = 7)),1)
    at = param_trans[[1]]$fun(labs)
    axis(side = 1, at = at, label = labs)

    points(params[PriorTest], ps[PriorTest], pch = 19, col = 'blue')
    points(params[!PriorTest], ps[!PriorTest], pch = 19, col = 'red')
    lines(x, y)

    dy = diff(y)
    turning = which(dy[-1] <0 & head(dy, -1) >0) + 1
    yturn = y[turning]; xturn = x[turning]
    turning = sort.int(yturn, decreasing = 2, index.return = TRUE)[[2]]
    yturn = yturn[turning] ;  xturn = xturn[turning]

    params = params[!PriorTest]; ps = ps[!PriorTest]; target = max(ps)
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
        PriorTest = dif >0
        dif = dif[PriorTest]
        liki = which.min(dif/diff(params)[PriorTest])
        
        return(mids[liki])   
    }
    likliMis = c()
    for (i in 1:length(xturn)) {
         
        likliMis = c(likliMis, checkLikelyMissing(c(xturn[1:i], likliMis), c(yturn[1:i], rep(target, length(likliMis)))))
       
    }
    
    likliMis = unique(likliMis)
    sapply(xturn, function(i) lines(c(i, i), c(-9E9, 9E9), lty = 2))
    sapply(likliMis, function(i) lines(c(i, i), c(-9E9, 9E9), lty = 3))
    xturn = round(param_trans[[1]]$funInverse(xturn), 2)
    
    likliMis = round(param_trans[[1]]$funInverse(likliMis), 2)
    
    plot(c(0, 1), c(0, 1), type = 'n', axes = FALSE, xlab = '', ylab = '')
    text.units(x = 0.01, y = 0.9, adj = 0,
         paste0("Best performance: ", outs[1,which.max(ps)],
                "; log(P(~beta~|Y) P(Y)) = ", round(max(ps))))
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
}
