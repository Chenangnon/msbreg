
# object of class 'msbm'
empirical.lambda <- function(object, scale.unit = TRUE,
                             min.quantile = bin.size,     # must be > eps and < 0.25
                             max.quantile = 3 * bin.size, # must be <= 0.5
                             bin.size = 0.1,
                             leftmost = TRUE,    # should we always test against the leftmost bin?
                             alpha = 0.2,         # Type I error
                             alternative = c("two.sided", "less", "greater"),
                             conf.level = 0.95,
                             eps = 1e-5, ...) {
  lambda0 <- object$lambda.values
  if (length(lambda0) > 1) {

    warning("lambda correction not implemented for vector valued object$lambda.values")

    return(structure(lambda0, lambda0 = lambda0))
  }
  theta <- object$coefficients
  mf <- object$frame %||% model.frame(object)
  X <- mf$input.matrix
  csigns <- sign(theta[mf$parindex$slopes[1:NCOL(X)]])
  csigns[csigns == 0] <- 1
  if (any(csigns > 0)) {
    X <- t(t(X) * csigns)
  }

  Linfer <- infer.lambda(x = X, y = mf$y, weights = mf$weights,
                         min.quantile = min.quantile, max.quantile = max.quantile,
                         bin.size = bin.size, leftmost = leftmost, alpha = alpha,
                         conf.level = conf.level, eps = eps)
  mu_ab <- Linfer$lambda
  pick <- Linfer$which
  npick <- length(pick)

  hvalues <- object$fitted.values[pick] / lambda0
  omega_ab <- mean(hvalues)
  omega2_ab <- var(hvalues) / npick

  Lhat <- max(min(mu_ab / (omega_ab - omega2_ab), 1-object$control$epsilon), object$control$epsilon)

  return(structure(Lhat,
                   lambda0 = lambda0,
                   mu_ab = mu_ab))
}
