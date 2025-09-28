
#' @rdname rsquared
#' @export NagelkerkeR2
NagelkerkeR2 <- function (object, ..., adjust.size = FALSE) {

  out <- rsquared (object, ..., method = 'Nagelkerke',
                   adjust.size = adjust.size)
  Call <- match.call()
  rownames(out) <- as.character(Call[-1L])[1:NROW(out)]

  return(out)
}

get.NR2VAL <- function (object, n, rank, adjust.size) {
  R2 <- (1 - exp((object$deviance - object$null.deviance)/n))/(1 - exp(-object$null.deviance/n))
  if (adjust.size)
    R2 <- 1 - (1 - R2) * (n - 1) / (n - rank)

  return(R2)
}
