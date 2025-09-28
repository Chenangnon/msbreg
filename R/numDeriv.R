# Functions from 'numDeriv' (21/02/2025)
# Saved here for internal use only, free of any modification of the package by the authors
# All credits to: Paul Gilbert and Ravi Varadhan
numDerivhessian <- function (func, x, method = "Richardson",
                    method.args = list(), ...) {
  if (1 != length(func(x, ...)))
    stop("Richardson method for hessian assumes a scalar valued function.")
  if (method == "complex") {
    args <- list(eps = 1e-04, d = 0.1, zero.tol = sqrt(.Machine$double.eps/7e-07),
                 r = 4, v = 2)
    args[names(method.args)] <- method.args
    fn <- function(x, ...) {
      numDerivgrad(func = func, x = x, method = "complex", side = NULL,
           method.args = list(eps = .Machine$double.eps),
           ...)
    }
    return(numDerivjacobian(func = fn, x = x, method = "Richardson",
                    side = NULL, method.args = args, ...))
  }
  else if (method != "Richardson")
    stop("method not implemented.")
  args <- list(eps = 1e-04, d = 0.1, zero.tol = sqrt(.Machine$double.eps/7e-07),
               r = 4, v = 2, show.details = FALSE)
  args[names(method.args)] <- method.args
  D <- numDerivgenD (func, x, method = method, method.args = args, ...)$D
  if (1 != nrow(D))
    stop("BUG! should not get here.")
  H <- diag(NA, length(x))
  u <- length(x)
  for (i in 1:length(x)) for (j in 1:i) {
    u <- u + 1
    H[i, j] <- D[, u]
  }
  H <- H + t(H)
  diag(H) <- diag(H)/2
  H
}

numDerivgrad <- function (func, x, method = "Richardson",
                  side = NULL, method.args = list(), ...) {
  f <- func(x, ...)
  n <- length(x)
  if (is.null(side))
    side <- rep(NA, n)
  else {
    if (n != length(side))
      stop("Non-NULL argument 'side' should have the same length as x")
    if (any(1 != abs(side[!is.na(side)])))
      stop("Non-NULL argument 'side' should have values NA, +1, or -1.")
  }
  case1or3 <- n == length(f)
  if ((1 != length(f)) & !case1or3)
    stop("numDerivgrad assumes a scalar valued function.")
  if (method == "simple") {
    args <- list(eps = 1e-04)
    args[names(method.args)] <- method.args
    side[is.na(side)] <- 1
    eps <- rep(args$eps, n) * side
    if (case1or3)
      return((func(x + eps, ...) - f)/eps)
    df <- rep(NA, n)
    for (i in 1:n) {
      dx <- x
      dx[i] <- dx[i] + eps[i]
      df[i] <- (func(dx, ...) - f)/eps[i]
    }
    return(df)
  }
  else if (method == "complex") {
    if (any(!is.na(side)))
      stop("method 'complex' does not support non-NULL argument 'side'.")
    eps <- .Machine$double.eps
    v <- try(func(x + eps * (0+1i), ...))
    if (inherits(v, "try-error"))
      stop("function does not accept complex argument as required by method 'complex'.")
    if (!is.complex(v))
      stop("function does not return a complex value as required by method 'complex'.")
    if (case1or3)
      return(Im(v)/eps)
    h0 <- rep(0, n)
    g <- rep(NA, n)
    for (i in 1:n) {
      h0[i] <- eps * (0+1i)
      g[i] <- Im(func(x + h0, ...))/eps
      h0[i] <- 0
    }
    return(g)
  }
  else if (method == "Richardson") {
    args <- list(eps = 1e-04, d = 1e-04, zero.tol = sqrt(.Machine$double.eps/7e-07),
                 r = 4, v = 2, show.details = FALSE)
    args[names(method.args)] <- method.args
    d <- args$d
    r <- args$r
    v <- args$v
    show.details <- args$show.details
    a <- matrix(NA, r, n)
    h <- abs(d * x) + args$eps * (abs(x) < args$zero.tol)
    pna <- (side == 1) & !is.na(side)
    mna <- (side == -1) & !is.na(side)
    for (k in 1:r) {
      ph <- mh <- h
      ph[pna] <- 2 * ph[pna]
      ph[mna] <- 0
      mh[mna] <- 2 * mh[mna]
      mh[pna] <- 0
      if (case1or3)
        a[k, ] <- (func(x + ph, ...) - func(x - mh, ...))/(2 *
                                                             h)
      else for (i in 1:n) {
        if ((k != 1) && (abs(a[(k - 1), i]) < 1e-20))
          a[k, i] <- 0
        else a[k, i] <- (func(x + ph * (i == seq(n)),
                              ...) - func(x - mh * (i == seq(n)), ...))/(2 *
                                                                           h[i])
      }
      if (any(is.na(a[k, ])))
        stop("function returns NA at ", h, " distance from x.")
      h <- h/v
    }
    if (show.details) {
      cat("\n", "first order approximations", "\n")
      print(a, 12)
    }
    for (m in 1:(r - 1)) {
      a <- (a[2:(r + 1 - m), , drop = FALSE] * (4^m) -
              a[1:(r - m), , drop = FALSE])/(4^m - 1)
      if (show.details & m != (r - 1)) {
        cat("\n", "Richarson improvement group No. ",
            m, "\n")
        print(a[1:(r - m), , drop = FALSE], 12)
      }
    }
    return(c(a))
  }
  else stop("indicated method ", method, "not supported.")
}

numDerivjacobian <- function (func, x, method = "Richardson",
                              side = NULL, method.args = list(), ...) {
  f <- func(x, ...)
  n <- length(x)
  if (is.null(side))
    side <- rep(NA, n)
  else {
    if (n != length(side))
      stop("Non-NULL argument 'side' should have the same length as x")
    if (any(1 != abs(side[!is.na(side)])))
      stop("Non-NULL argument 'side' should have values NA, +1, or -1.")
  }
  if (method == "simple") {
    args <- list(eps = 1e-04)
    args[names(method.args)] <- method.args
    side[is.na(side)] <- 1
    eps <- rep(args$eps, n) * side
    df <- matrix(NA, length(f), n)
    for (i in 1:n) {
      dx <- x
      dx[i] <- dx[i] + eps[i]
      df[, i] <- (func(dx, ...) - f)/eps[i]
    }
    return(df)
  }
  else if (method == "complex") {
    if (any(!is.na(side)))
      stop("method 'complex' does not support non-NULL argument 'side'.")
    eps <- .Machine$double.eps
    h0 <- rep(0, n)
    h0[1] <- eps * (0+1i)
    v <- try(func(x + h0, ...))
    if (inherits(v, "try-error"))
      stop("function does not accept complex argument as required by method 'complex'.")
    if (!is.complex(v))
      stop("function does not return a complex value as required by method 'complex'.")
    h0[1] <- 0
    jac <- matrix(NA, length(v), n)
    jac[, 1] <- Im(v)/eps
    if (n == 1)
      return(jac)
    for (i in 2:n) {
      h0[i] <- eps * (0+1i)
      jac[, i] <- Im(func(x + h0, ...))/eps
      h0[i] <- 0
    }
    return(jac)
  }
  else if (method == "Richardson") {
    args <- list(eps = 1e-04, d = 1e-04, zero.tol = sqrt(.Machine$double.eps/7e-07),
                 r = 4, v = 2, show.details = FALSE)
    args[names(method.args)] <- method.args
    d <- args$d
    r <- args$r
    v <- args$v
    a <- array(NA, c(length(f), r, n))
    h <- abs(d * x) + args$eps * (abs(x) < args$zero.tol)
    pna <- (side == 1) & !is.na(side)
    mna <- (side == -1) & !is.na(side)
    for (k in 1:r) {
      ph <- mh <- h
      ph[pna] <- 2 * ph[pna]
      ph[mna] <- 0
      mh[mna] <- 2 * mh[mna]
      mh[pna] <- 0
      for (i in 1:n) {
        a[, k, i] <- (func(x + ph * (i == seq(n)), ...) -
                        func(x - mh * (i == seq(n)), ...))/(2 * h[i])
      }
      h <- h/v
    }
    for (m in 1:(r - 1)) {
      a <- (a[, 2:(r + 1 - m), , drop = FALSE] * (4^m) -
              a[, 1:(r - m), , drop = FALSE])/(4^m - 1)
    }
    return(array(a, dim(a)[c(1, 3)]))
  }
  else stop("indicated method ", method, "not supported.")
}

numDerivgenD <- function (func, x, method = "Richardson",
                          method.args = list(), ...) {
  if (method != "Richardson")
    stop("method not implemented.")
  args <- list(eps = 1e-04, d = 1e-04, zero.tol = sqrt(.Machine$double.eps/7e-07),
               r = 4, v = 2)
  args[names(method.args)] <- method.args
  d <- args$d
  r <- args$r
  v <- args$v
  if (v != 2)
    stop("The current code assumes v is 2 (the default).")
  f0 <- func(x, ...)
  n <- length(x)
  h0 <- abs(d * x) + args$eps * (abs(x) < args$zero.tol)
  D <- matrix(0, length(f0), (n * (n + 3))/2)
  Daprox <- matrix(0, length(f0), r)
  Hdiag <- matrix(0, length(f0), n)
  Haprox <- matrix(0, length(f0), r)
  for (i in 1:n) {
    h <- h0
    for (k in 1:r) {
      f1 <- func(x + (i == (1:n)) * h, ...)
      f2 <- func(x - (i == (1:n)) * h, ...)
      Daprox[, k] <- (f1 - f2)/(2 * h[i])
      Haprox[, k] <- (f1 - 2 * f0 + f2)/h[i]^2
      h <- h/v
    }
    for (m in 1:(r - 1)) for (k in 1:(r - m)) {
      Daprox[, k] <- (Daprox[, k + 1] * (4^m) - Daprox[,
                                                       k])/(4^m - 1)
      Haprox[, k] <- (Haprox[, k + 1] * (4^m) - Haprox[,
                                                       k])/(4^m - 1)
    }
    D[, i] <- Daprox[, 1]
    Hdiag[, i] <- Haprox[, 1]
  }
  u <- n
  for (i in 1:n) {
    for (j in 1:i) {
      u <- u + 1
      if (i == j)
        D[, u] <- Hdiag[, i]
      else {
        h <- h0
        for (k in 1:r) {
          f1 <- func(x + (i == (1:n)) * h + (j == (1:n)) *
                       h, ...)
          f2 <- func(x - (i == (1:n)) * h - (j == (1:n)) *
                       h, ...)
          Daprox[, k] <- (f1 - 2 * f0 + f2 - Hdiag[,
                                                   i] * h[i]^2 - Hdiag[, j] * h[j]^2)/(2 * h[i] *
                                                                                         h[j])
          h <- h/v
        }
        for (m in 1:(r - 1)) for (k in 1:(r - m)) Daprox[,
                                                         k] <- (Daprox[, k + 1] * (4^m) - Daprox[, k])/(4^m -
                                                                                                          1)
        D[, u] <- Daprox[, 1]
      }
    }
  }
  D <- list(D = D, p = length(x), f0 = f0, func = func, x = x,
            d = d, method = method, method.args = args)
  class(D) <- "Darray"
  invisible(D)
}
