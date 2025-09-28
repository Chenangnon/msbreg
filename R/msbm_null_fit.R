
msbm.null.fit <- function (frame, y = frame$y,
                           link = logit(),
                           control = list()) {
  stopifnot(length(y) == NROW(frame$input.matrix))
  control <- do.call("msbm.control", control)
  ynm <- names(y)
  frame$y <- as.numeric(y)
  names(frame$y) <- ynm

  if (length(frame$sample.weights) > 1) {
    frame$sample.weights <- frame$dims$nobs *
      frame$sample.weights / sum(frame$sample.weights)
  }

  #* Set dimensions from mbm.frame: q, pj (for j = 1:q), d, r
  eval(get.msbm.dims())
  weights <- rep(frame$weights, length.out = nobs)
  n.ok <- nobs - sum(weights * frame$sample.weights == 0)

  #* Residual deviance function
  eval (toget_dev.resids())

  #* Null deviance calculation
  theta <- numeric(dims$npars)
  eval (toget_null.msbm.fit())

  # Format output
  out <- list(null.deviance = nulldev,
              null.mu = nullmu,
              null.rank = nullrank)

  return(out)
}
