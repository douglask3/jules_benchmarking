
modParamTrans.zeroInf <- function(x, a = -1) 1-exp(a*x)

modParamTrans.zeroInf.inverse <- function(y, a = -1) (1/a)* log(1-y)


