optiPlot1D <- function() {    
    
    x = seq(0, 0.9999999999999, paramDetail)
    xf = logit(x)
    xp = param_trans[[1]]$funInverse(x)
    params = outs[1,]
    
    PriorTest = is.na(outs[2,]) | outs[2,]
    
    #if (ncol(outs) <=3) outs[2,PriorTest] = TRUE else {
    if (sum(!PriorTest) >=3) {
        outs = outs[, !PriorTest]
        ps = ps[!PriorTest]
        params = params[!PriorTest]
    }
    
    nps = length(ps)  
    paramsP = params
    params = param_trans[[1]]$fun(params)  
    paramsF = logit(params)

    weightedApprox <- function(x) {
        w = paramsF-x
        findP <- function(test, FUN) {
            if (!any(test)) return(NULL)
            p = FUN(w[test])
            c(ps[test][p], abs(w[test][p]), paramsF[test][p])
        }
        
        p1 = findP(w < 0, which.max); p2 = findP(w >=0, which.min)
        #browser()
        if(is.null(p1)) return(p2[1]) else if(is.null(p2)) return(p1[1])
        #return(((1/p1[2]) + (1/p2[2])))
        ((p1[1]/p1[2]) + (p2[1]/p2[2]))/((1/p1[2]) + (1/p2[2]))
    }
    
    By = sapply(xf, weightedApprox)
    Py = log(do.call(param_trans$prior[[1]], list(xp, param_trans$prior[-1])))
    y = By + Py
    
    ps = ps + log(do.call(param_trans$prior[[1]], list(paramsP, param_trans$prior[-1])))
    ps[is.infinite(ps)] = min(y[!is.infinite(y)]) -1
    ## new samples
    P0 = exp(y - min(y[!is.infinite(y)])) 
    
    Y = X = rep(NaN, nNewParams)
    for (i in 1:nNewParams) {
        
        P = P0 * sapply(x, function(j) min(abs(c(X, params) - j), na.rm = TRUE))
        
        id = sample(1:length(x), size = 1, prob = P)
        X[i] = x[id]
        Y[i] = y[id]
    }
        
    #X = sample(1:length(x), size = nNewParams, replace = FALSE, prob = P)
    #Y = y[X]; X = x[X]
    
    
    xrange = range(c(params, x)); yrange = range(c(ps, y[!is.infinite(y)]))
    par(mfrow = c(2, 1))
    plot(xrange, yrange, type = 'n', xaxt = 'n', xlab = '', ylab = '')

    mtext.units(side = 2, line = 2, 'log(P(~beta~|Y) P(Y))')
    mtext.units(side = 1, line = 2, '~beta~')
    
    labs = signif(param_trans[[1]]$funInverse(seq(0, 1, length.out = 7)),1)
    at = param_trans[[1]]$fun(labs)
    axis(side = 1, at = at, label = labs)

    points(params[PriorTest], ps[PriorTest], pch = 19, col = 'blue')
    points(params[!PriorTest], ps[!PriorTest], pch = 19, col = 'red')

    yt = y[!is.infinite(y)]
    rescale2By <- function(x) {
        test = !is.infinite(x)
        minX = min(x[test])
        x = x - minX
        x = x / diff(range(x[test]))
        x = x * diff(range(yt)) + min(yt)
        x[!test] = min(yt) - 1
        x
    }
    lines(x, rescale2By(y))
    lines(x, rescale2By(Py), lty = 2, col = 'blue')
    lines(x, rescale2By(By), lty = 2, col = 'red')

    mapply(function(x, wd) lines(c(x, x), c(-9E9, 9E9), lty = 3, lwd = wd, col = '#00000099'),
           X, seq(2, 0.1, length.out = length(X)))
    
    par(mar = rep(0, 4))
    plot(c(0, 1), c(0, 1), type = 'n', axes = FALSE, xlab = '', ylab = '')

    points(c(0.1, 0.1), c(0.98, 0.93), pch = 19, col = c('blue', 'red'))
    text.units(c(0.12), c(0.98), c('Seed'), adj = 0)
    text.units(c(0.12), c(0.93), c('JULES run'), adj = 0)

    addLegLine <- function(x, y = 0.98, txt, ...) {
        lines(c(x, x + 0.07), c(y, y), ...)
        text.units(x + 0.1, y, txt, adj = 0)
    }
    addLegLine(0.4, txt = 'P(~beta~)', col = 'blue', lty = 2)
    addLegLine(0.65, txt = 'P(Y|~beta~)', col = 'red', lty = 2)
    addLegLine(0.4, y = 0.93, 'P(~beta~| Y)')

    lines(c(0.9, 0.9), c(0.99, 0.8), lty = 3, col = '#00000099')
    text(0.93, 0.9, 'New samples', srt = 90, xpd = NA)
    
    
    text.units(x = 0.01, y = 0.83, adj = 0,
         paste0("Best performance: ", outs[1,which.max(ps)],
                "; log(P(~beta~|Y) P(Y)) = ", round(max(ps))))
    
    
    text(x = 0.01, y = 0.75, adj = 0, "New Parameters to test:")   
    X =  round(param_trans[[1]]$funInverse(X), 2)
    if (length(X) != length(unique(X))) browser()
    tab1 = cbind(1:length(X), round(X, 2), round(Y, 2))
    
    AddNews <- function(x, i1, i2) {
        i = i1:i2
        text(x = x + 0.05, y = 0.65, adj = 0, "Param.")
        text.units(x = x + 0.15, y = 0.65, adj = 0, "log(P(~beta~|Y))")
    
        textTfun <- function(i, y, ...)
            text(i, adj = 0, y = y - seq(0, by = 0.05, length.out = length(i)),...)
        textTfun(paste0(tab1[i,1], '.'), , x = x, y = 0.5)
        textTfun(tab1[i,2], x = 0.05 + x, y = 0.5)
        textTfun(tab1[i,3], x = 0.15 + x, y = 0.5)
    
    }
    i1 = seq(1, nNewParams, by = 10)
    i2 = c(tail(i1, -1)-1, nNewParams)
    xS = seq(0, by = 0.33, length.out = length(i1))
    mapply(AddNews, xS, i1, i2)
    browser()
}
