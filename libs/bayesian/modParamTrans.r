
modParamTrans.zeroInf <- function(x, a = -1) 1-exp(a*x)

modParamTrans.zeroInf.inverse <- function(y, a = -1) {
    y[y ==1] = 1 - 1E-10
    (1/a)* log(1-y)
}

modParamTrans.zeroOne <- function(x, ..) x

modParamTrans.zeroOne.inverse <- function(x, ..) x
