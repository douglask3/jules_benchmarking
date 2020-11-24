
minusStandard <- function(r1, r2) {
    if (nlayers(r1) == 1) r1 = r1[[1]]
    if (nlayers(r2) == 1) r2 = r2[[1]]
    return(r1 - r2)
}
