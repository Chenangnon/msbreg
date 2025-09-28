# Retrieve MBM dimensions from a model frame
# is.msbm.frame
get.msbm.dims <- function () {
  expression ( {
    dims <- frame$dims
    nobs <- dims$nobs
    npars <- dims$npars
    q <- dims$q
    p <- dims$p
    pj <- dims$pj
    d <- dims$d
    r <- dims$r
    intercepts <- frame$intercepts
  } )
}
