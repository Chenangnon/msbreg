# Functions from 'mirt::numerical_deriv' (21/02/2025)
# Saved here for internal use only, free of any modification of the package by the authors
# All credits to: Phil Chalmers <rphilip.chalmers@gmail.com>
mirtJacobian <- function (x, f, ..., delta = 1e-05,
                     type = "Richardson", r = 4L) {
  if (!length(x)) {
    return(numeric())
  }

  if (is.character(f))
    f <- get(f, mode = 'function', envir = parent.frame())

  if (type == "central") {
    mirtCentral_difference(par = x, f = f, delta = delta,
                              ...)
  }
  else if (type == "forward") {
    mirtForward_difference(par = x, f = f, delta = delta,
                              ...)
  }
  else if (type == "Richardson") {
    mirtRichardson(par = x, f = f, delta = delta * 10,
                      r = r, ...)
  }
}

mirtHessian <- function (x, f, ..., delta = 1e-05,
                      type = "Richardson", r = 4L) {
  if (!length(x)) {
    return(matrix(numeric()))
  }

  if (is.character(f))
    f <- get(f, mode = 'function', envir = parent.frame())

  if (type == "central") {
    mirtCentral_difference2(par = x, f = f, delta = delta,
                               ...)
  }
  else if (type == "forward") {
    mirtForward_difference2(par = x, f = f, delta = delta,
                               ...)
  }
  else if (type == "Richardson") {
    mirtRichardson2(par = x, f = f, delta = delta * 1000,
                       r = r, ...)
  }
}

mirtForward_difference <- function(par, f, delta, ...) {
  dots <- list(...)
  np <- length(par)
  g <- numeric(np)
  if (is.null(dots$ObJeCtIvE))
    fx <- f(par, ...)
  else fx <- dots$ObJeCtIvE
  for (i in seq_len(np)) {
    p <- par
    p[i] <- p[i] + delta
    g[i] <- (f(p, ...) - fx)/delta
  }
  g
}

mirtForward_difference2 <- function(par, f, delta, ...) {
  dots <- list(...)
  np <- length(par)
  hess <- matrix(0, np, np)
  if (is.null(dots$ObJeCtIvE))
    fx <- f(par, ...)
  else fx <- dots$ObJeCtIvE
  fx1 <- numeric(np)
  for (i in seq_len(np)) {
    tmp <- par
    tmp[i] <- tmp[i] + delta
    fx1[i] <- f(tmp, ...)
  }
  for (i in seq_len(np)) {
    for (j in i:np) {
      fx1x2 <- par
      fx1x2[i] <- fx1x2[i] + delta
      fx1x2[j] <- fx1x2[j] + delta
      hess[i, j] <- hess[j, i] <- (f(fx1x2, ...) -
                                     fx1[i] - fx1[j] + fx)/(delta^2)
    }
  }
  (hess + t(hess))/2
}

mirtCentral_difference <- function(par, f, delta, ...) {
  np <- length(par)
  g <- numeric(np)
  for (i in seq_len(np)) {
    p1 <- p2 <- par
    p1[i] <- p1[i] + delta
    p2[i] <- p2[i] - delta
    g[i] <- (f(p1, ...) - f(p2, ...))/(2 * delta)
  }
  g
}

mirtCentral_difference2 <- function(par, f, delta, ...) {
  np <- length(par)
  hess <- matrix(0, np, np)
  fx <- f(par, ...)
  for (i in seq_len(np)) {
    for (j in i:np) {
      if (i == j) {
        p <- par
        p[i] <- p[i] + 2 * delta
        s1 <- f(p, ...)
        p[i] <- p[i] - 4 * delta
        s3 <- f(p, ...)
        hess[i, i] <- (s1 - 2 * fx + s3)/(4 * delta^2)
      }
      else {
        p <- par
        p[i] <- p[i] + delta
        p[j] <- p[j] + delta
        s1 <- f(p, ...)
        p[j] <- p[j] - 2 * delta
        s2 <- f(p, ...)
        p[i] <- p[i] - 2 * delta
        s4 <- f(p, ...)
        p[j] <- p[j] + 2 * delta
        s3 <- f(p, ...)
        hess[i, j] <- hess[j, i] <- (s1 - s2 - s3 +
                                       s4)/(4 * delta^2)
      }
    }
  }
  (hess + t(hess))/2
}

mirtRichardson <- function(par, f, delta, r = 4L, ...) {
  R0 <- R1 <- matrix(0, length(par), r)
  R0[, 1L] <- mirtCentral_difference(par = par, f = f, delta = delta,
                                 ...)
  for (i in 1L:(r - 1L)) {
    delta <- delta/2
    R1[, 1L] <- mirtCentral_difference(par = par, f = f,
                                   delta = delta, ...)
    for (j in 1L:i) R1[, j + 1] <- (4^j * R1[, j] - R0[,
                                                       j])/(4^j - 1)
    R0 <- R1
  }
  R1[, i + 1]
}

mirtRichardson2 <- function(par, f, delta, r = 4L, ...) {
  R0 <- R1 <- matrix(0, length(par)^2, r)
  R0[, 1L] <- as.vector(mirtCentral_difference2(par = par,
                                            f = f, delta = delta, ...))
  for (i in 1L:(r - 1L)) {
    delta <- delta/2
    R1[, 1L] <- as.vector(mirtCentral_difference2(par = par,
                                              f = f, delta = delta, ...))
    for (j in 1L:i) R1[, j + 1] <- (4^j * R1[, j] - R0[,
                                                       j])/(4^j - 1)
    R0 <- R1
  }
  hess <- matrix(R1[, i + 1], length(par), length(par))
  (hess + t(hess))/2
}
