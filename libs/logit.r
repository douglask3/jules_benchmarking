logit <- function(x, mn = 1E-50) {
    #browser()
    x[x < mn] = mn
    x[x > (1-mn)] = 1-mn

    out = log(x/(1-x))
    out[out > -log(mn)] = -log(mn)
    out
}
