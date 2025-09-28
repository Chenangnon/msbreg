
.fisher.matrix.x.y <- function (x, y,
                                weights = 1,
                                offset = 0,
                                sigma = NULL,
                                ...) {
  # Remove NAs from x, y, weights, and offset (if any)
  eval(drop.NAs.from.x.y.weights.offset())

  if (is.null(sigma)) {
    if (length(y)) {
      fitwlm <- stats::lm.wfit(x = x, y = y, w = weights, offset = offset)
      rss <- sum(weights * fitwlm$residuals^2)
      sigma <- sqrt(rss / fitwlm$df.residual)
    }
    else {
      sigma <- 1
    }
  }
  else {
    stopifnot(is.numeric(sigma))
    stopifnot(sigma > 0)
  }

  out <- t(x) %*% (x * weights)
  out <- out / (sigma^2)

  return(structure(out,
                   sigma = sigma,
                   class = c("fisher.matrix",
                             class(out))))
}
