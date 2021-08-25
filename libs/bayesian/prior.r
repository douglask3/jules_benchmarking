prior.Uniform <- function(x, ps) {
    y = rep(1, length(x))
    y[x < ps[[1]]] = 0
    y[x > ps[[2]]] = 0
    return(y)
}

prior.LogNormal <- function(x, ps) {   
    x = log(x)    
    mu = log(ps[[1]])
    y = dnorm(x, mu, ps[[2]])
    return(y)
}
