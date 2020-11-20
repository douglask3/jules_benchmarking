Modalise <- function(r) {
    modal = testModal(r)
    modal_approx = layer.apply(2:6, function(i)
                              modal[[i]] * (0.5-cos(2*pi * modal[[i+6]]/12)/2))
    modal_approx = 1 + sum(modal_approx, na.rm = TRUE)/modal[[1]]
    return(modal_approx)
}
