checkPriorPoints <- function(outs, ps, nthresh) {
    PriorTest = is.na(outs[nrow(outs),]) | outs[nrow(outs),]
    
    if (sum(!PriorTest) >= nthresh) {
        outs = outs[, !PriorTest]
        ps = ps[!PriorTest]
        #params = params[!PriorTest]
        PriorTest = PriorTest[!PriorTest]
    }    
    
    paramsP = params = outs[1:(nrow(outs)-1),] 
    if ((nrow(params)) != length(param_trans)) {
        print("param info doesnt not match params from runs")
        browser()
    }
    for (i in 1:length(param_trans)) 
        params[i,] = param_trans[[i]][[1]]$fun(params[i,]) 
    
    paramsF = logit(params)
    return(list(paramsP, params, paramsF, PriorTest))
}

optiPlot2D <- function() {    
  
    #######################################################
    ## removing proior points if we have enough new ones ##
    #######################################################    
    c(paramsP, params, paramsF, PriorTest) := checkPriorPoints(outs, ps, 4)
    
    ##############################
    ## interoplating P(Y|beta)  ##
    ##############################
    x  = seq(0, 1.0, paramDetail)
    xf = logit(x)
    xp = lapply(param_trans, function(pt) pt[[1]]$funInverse(x))

    weightedApprox <- function(x) {
        if (length(ps) == 1) return(ps)
        w = apply(paramsF, 2, function(i) sqrt(sum((i - x)^2)))
        return(ps[which.min(w)])
        ## Just neighest neightbour for now. Will work on thi later.
        index = sort.int(w, index.return=TRUE)[[2]]
        w = w[index]
        paramsF = paramsF[,index]
        ps = ps[index]
        
        angleW <- function(i) {
            A = x
            B = paramsF[,i]
            ab = B - A
            AB = sqrt(sum(ab^2))
            calAngle <- function(C) {
                bc = B - C
                BC =  sqrt(sum(bc^2))
                
                angle = acos(dot(ab, bc)/(AB*BC))
                if (is.na(angle)) browser()
                browser()
            }
            angles = apply(paramsF[,1:(i-1)], 2, calAngle)
            browser()
        }
        w[4:length(w)] = sapply(4:length(w), angleW)
        browser()
        
    }
    xfs = cbind(xf, rep(xf, each = length(xf)))
    
    By = apply(xfs, 1, weightedApprox)
    
    CalPrior <- function(pt, xp) log(do.call(pt$prior[[1]], list(xp, pt$prior[-1])))   
    Py = mapply(CalPrior, param_trans, xp)
    Py = unlist(lapply(Py[,2], function(p) p + Py[,1]))
    
    y = By + Py
    
    ##########################
    ## finding new samples  ##
    ##########################
    
    P0 = exp(y - min(y[!is.infinite(y)])) 
    P0[is.infinite(P0)] = 0
    
    xs = cbind(x, rep(x, each = length(x)))
    XY = matrix(NaN, nrow = nNewParams, ncol = 3)
    for (i in 1:nNewParams) {
        weigthDist <- function(j) {
            plist = rbind(XY[,-1], t(params))
            plist = plist[!apply(is.na(plist), 1, any),]
            dist = apply(plist, 1, function(pl) sqrt(sum((pl - j)^2)))
            
            min(dist, na.rm = TRUE)
        } 
        P = P0 * apply(xs, 1, weigthDist)    
        id = sample(1:length(P), size = 1, prob = P)
        XY[i,] = c(y[id], xs[id,])
    }
    ##########################
    ## plotting             ##
    ##########################     
    ps0 = ps
    matrix2list <- function(x) 
        lapply(apply(x, 2, list), function(i) i[[1]])
        
    
    pntPriors = mapply(CalPrior, param_trans, matrix2list(t(paramsP)))
    
    ps = ps + apply(pntPriors, 1, sum)
    
    ps[is.infinite(ps)] = min(y[!is.infinite(y)])
    
    ## setting up plot
    #xrange = range(c(params, x)); yrange = range(c(ps, y[!is.infinite(y)]))
    dev.new()
    layout(rbind(c(1, 3, 5), c(2, 4, 6), 7), height = c(0.7, 0.3, 1))
    par(mar = c(2.5, .5, 0.5, 0.5), oma = c(0, 3, 1.5, 0))
    
    imageP <- function(Pz, psZ = NULL, cols, title = '', 
                       xaxt = TRUE, xlbt = FALSE, yaxt = FALSE, newbies = FALSE) {
        Pz[is.infinite(Pz)] = NaN
        limits = quantile(Pz, seq(0, 1, by = 0.1), na.rm = TRUE)
        
        limits[1] = -9E9; limits[length(limits)] = 9E9
        cols = make_col_vector(cols, limits = limits[1:10])
        image(x, x, matrix(Pz, ncol = length(x)), breaks = limits, col = cols, 
              xaxt = 'n', yaxt = 'n')
        mtext.units(side = 3, adj = 0.5, title)
        addAxis <- function(side, pt) {
            labs = signif(pt[[1]]$funInverse(seq(0, 1, length.out = 7)),1)
            at = pt[[1]]$fun(labs)
            axis(side = side, at = at, label = labs)
        }
        if (xaxt) addAxis(1, param_trans[[1]])
        if (xlbt) mtext(names(param_trans)[1], line = 2)
        if (yaxt) { addAxis(2, param_trans[[2]]); mtext(names(param_trans)[2], line = 2)}
        if (newbies) points(XY[,2], XY[,3], pch = 17)
        if (!is.null(psZ)) {
            pch = c(19, 4)[PriorTest+1]
            points(params[1,], params[2,], pch = pch)
            points(params[1,], params[2,], pch = pch, 
                   cex = 0.8, col = cols[cut_results(psZ, limits)-1])
        }
        plot(c(0, 1), c(0, 1), axes = FALSE, xlab = '', ylab = '', type = 'n')
        add_raster_legend2(cols = cols, limits = round(head(limits[-1], -1), 2), 
                           transpose = FALSE, 
                           srt = 45, plot_loc = c(0.0, 1.0, 0.8,0.9))
    } 
 
    imageP(Py, NULL, cols = c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4',
                             '#1d91c0','#225ea8','#0c2c84'), 'log(P(~beta~)', yaxt = TRUE)    
    
    imageP(By, ps0, cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c',
                             '#fc4e2a','#e31a1c','#b10026'), 'log(P(Y|~beta~)', xlbt = TRUE)

    imageP(y, ps, cols = c('#fff7f3','#fde0dd','#fcc5c0','#fa9fb5','#f768a1',
                           '#dd3497','#ae017e','#7a0177'), 
           'log(P(~beta~|Y) P(Y))', newbies = TRUE)
    
    browser()
plot(xrange, yrange, type = 'n', xaxt = 'n', xlab = '', ylab = '')

    mtext.units(side = 2, line = 2, 'log(P(~beta~|Y) P(Y))')
    mtext.units(side = 1, line = 2, '~beta~')
    
    labs = signif(param_trans[[1]]$funInverse(seq(0, 1, length.out = 7)),1)
    at = param_trans[[1]]$fun(labs)
    axis(side = 1, at = at, label = labs)

    ## adding samples
    points(params[PriorTest], rescale2By(ps0[PriorTest]), pch = 4, col = 'red')
    points(params[!PriorTest], rescale2By(ps0[!PriorTest]), pch = 19, col = 'red')
    points(params[PriorTest], ps[PriorTest], pch = 4, col = 'black')
    points(params[!PriorTest], ps[!PriorTest], pch = 19, col = 'black')

    ## adding mapped postiriors
    
    lines(x, rescale2By(y))
    lines(x, rescale2By(Py), lty = 2, col = 'blue')
    lines(x, rescale2By(By), lty = 2, col = 'red')

    mapply(function(x, wd) lines(c(x, x), c(-9E9, 9E9), lty = 3, lwd = wd, col = '#00000099'),
           X, seq(2, 0.1, length.out = length(X)))

    ############
    ## legend ##
    ############  
    par(mar = rep(0, 4))
    plot(c(0, 1), c(0, 1), type = 'n', axes = FALSE, xlab = '', ylab = '')

    
    text.units(c(0.03), 1.05, c('P(Y|~beta~)'), adj = 0.67, xpd = NA)
    text.units(c(0.1), 1.05, c('P(~beta~|Y)'), adj = 0.33, xpd = NA)
    points(c(0.03, 0.03), c(0.98, 0.93), pch = c(4, 19), col = 'red')
    points(c(0.1, 0.1), c(0.98, 0.93), pch = c(4, 19), col = 'black')
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
    
    ############################
    ## sample/new sample info ##
    ############################  
    text.units(x = 0.01, y = 0.83, adj = 0,
         paste0("Best performance: ", outs[1,which.max(ps)],
                "; log(P(~beta~|Y) P(Y)) = ", round(max(ps), 2)))
    
    
    text(x = 0.01, y = 0.75, adj = 0, "New Parameters to test:")   
    X =  round(param_trans[[1]]$funInverse(X), 2)
    if (length(X) != length(unique(X))) browser()
    tab1 = cbind(1:length(X), round(X, 2), round(Y, 2))
    
    AddNews <- function(x, i1, i2) {
        i = i1:i2
        text(x = x + 0.05, y = 0.67, adj = 0, "Param.")
        text.units(x = x + 0.15, y = 0.67, adj = 0, "log(P(~beta~|Y) P(Y))")
    
        textTfun <- function(i, y, ...)
            text(i, adj = 0, y = y - seq(0, by = 0.05, length.out = length(i)),...)
        textTfun(paste0(tab1[i,1], '.'), , x = x, y = 0.6)
        textTfun(tab1[i,2], x = 0.05 + x, y = 0.6)
        textTfun(tab1[i,3], x = 0.15 + x, y = 0.6)
    
    }
    i1 = seq(1, nNewParams, by = 10)
    i2 = c(tail(i1, -1)-1, nNewParams)
    xS = seq(0, by = 0.33, length.out = length(i1))
    mapply(AddNews, xS, i1, i2)
}

