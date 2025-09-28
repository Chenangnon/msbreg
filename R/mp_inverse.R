# Used by "score.test.msbm"
mp.inverse <- function (X, eps = NULL) {
  stopifnot(is.numeric(X))
  X <- as.matrix(X)

  # Perform singular value decomposition
  s <- svd(X)
  d <- s$d

  # Number of singular values (min(dim(X)))
  m <- length(d)

  # Check 'eps'
  if (missing(eps) | is.null(eps)) {
    eps <- .Machine$double.eps * max(NROW(X), NCOL(X)) * max(d)
    eps <- max(min(eps, 1e-2), .Machine$double.eps^4)
  }
  else {
    stopifnot(is.numeric(eps) | is.null(eps))
    stopifnot(eps > 0, eps < 1e-2)
  }

  # remove singular values ~ zero
  d <- d[d > eps]
  n <- length(d)

  # Inverse of positive singular values
  inv <- if (n > 1) diag(1/d) else 1/d

  # Add rows, columns and rows of zeros if X has null singular values
  if (n != m) {
    inv <- cbind(inv, matrix(0, nrow=n, ncol=(m - n)))
    inv <- rbind(inv, matrix(0, nrow=(m-n), ncol=m))
  }

  # compute the Moore-Penrose inverse
  inv <- s$v %*% inv %*% t(s$u)

  # set very small values to zero
  inv[abs(inv) < .Machine$double.eps] <- 0 # Use eps?


  return(structure(inv, m = m, rank = n, singular.values = s$d))
}
