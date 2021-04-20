trianglePlot <- function(xs, ys, cols, name, add = FALSE) {
    x = xs + 0.5 * ys
    y = ys  
    if (add) FUN = points else FUN = plot
    if (length(x) == 0) {
        FUN(0, 0, pch = 20, cex = 2, col = cols, axes = FALSE, type = 'n',
             xlab = '', ylab = '', xlim = c(0, 1), ylim = c(0,1))
    } else {
        cols = cols[unlist(mapply(rep, 1:9, 9 + (1:9)^3))]
        cols = densCols(x,y, colramp = colorRampPalette(cols))
        FUN(y~x, pch = 20, cex = 2, col = cols, axes = FALSE,
             xlab = '', ylab = '', xlim = c(0, 1), ylim = c(0,1))
    }
    
    fitT = triangularRegressoion(xs, ys)
    fitL = triangleLinearReg(xs, ys) 

    x = seq(0, 1, 0.2)
    lines(c(0, 1, 0.5, 0), c(0, 0, 1,0))
    lapply(x, function(i) {
           lines(c(0.5*i, 1-0.5*i), c(i,i), lty = 2, col = '#00000099');
           text(x=1-0.5*i+0.05, y=i, paste(i*100, ''), xpd = TRUE, srt = -60)})

    lapply(x, function(i) {
           lines(c(i, 0.5+0.5*i), c( 0, 1-i), lty = 3, col = '#00000099');
           text(x=i, y=-0.07, paste(i*100, ''), xpd = TRUE, srt = 0)})

    lapply(x, function(i) {
           lines(c(0.5-0.5*i, 1-i), c( 1-i, 0), lty = 4, col = '#00000099');
           text(x=0.45-0.5*i, y=1-i, paste(i*100, ''), xpd = TRUE, srt = 60)})
}

triangleLinearReg <- function(xs, ys) {
    test = !is.na(xs+ys)
    xs = xs[test]; ys = ys[test]; ys2 = 1-xs -ys
    cutLin <- function(x, a, b) {
        y = a*x + b
        y[y<0] = 0
        y[y>1] = 0
        return(y)
    }
    
    xnew = xnew0 = seq(0, 1, 0.001)
    fitY <- function(yi, a = -1, b = 1) 
         nls(yi ~cutLin(xs, a, b), start = list(a = a, b = b))
    
    r2 <- function(yi, fit) 
        if (class(fit) == "try-error") return(0) else return(cor(yi, predict(fit))^2)

    fit  = try(fitY(ys )); rs  = r2(ys , fit )
    fit2 = try(fitY(ys2)); rs2 = r2(ys2, fit2)
    
    tFUN <- function(y) y
    if (r2(ys2, fit2) > r2(ys, fit)) {
        fit = fit2
        rs = rs2
        tFUN <- function(z) 1-xnew-z
    }
    #if (round(rs, 2) == 0.07) browser()
    ynew = tFUN(predict(fit, newdata = list(xs = xnew)))
    
    text.units(x = 0.9, y = 0.95, paste0("Linear ~R2~: ", round(rs, 2)), 
               adj = 0, xpd = NA, srt = -60)
    lines(xnew+0.5*ynew, ynew, srt = -60)
}

triangularRegressoion <- function(xs, ys) {
    A = (ys/xs)
    A[is.na(A)] = 1
    A[A>1] = 1/A[A>1]
    A = A+1

    logit <- function(x) {
        x[x<0.000001] = 0.000001
        x[x>0.999999] = 0.999999
        log(x/(1-x))
    }
    xf = logit(xs *A); yf = logit(ys*A)
    xf = xf[!(ys == 0 & xs == 1)]; yf = yf[!(ys == 0 & xs == 1)]
    fit = lm(yf~xf)
    xnew = xnew0 = seq(0, 1, 0.001)
    xnewt = logit(xnew)
    ynewt = predict(fit, newdata = data.frame(xf = xnewt))
    ynew = ynew0 = 1/(1+exp(-ynewt))
    
    A = (ynew/xnew)
    A[A>1] = 1/A[A>1]
    A = A+1
    xnew = xnew /A; ynew = ynew/A

    lines(xnew+0.5*ynew, ynew, lwd = 2)
    text.units(x = 0.8, y = 0.9, paste0("Triangle ~R2~: ", round(summary(fit)$r.squared,2)), adj = 0, xpd = NA, srt = -60)

    
    return(fit)
}
