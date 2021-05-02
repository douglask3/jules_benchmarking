prior.Uniform <- function(x, ps) {
    y = rep(1, length(x))
    y[x < ps[[1]]] = 0
    y[x > ps[[2]]] = 0
    return(y)
}
