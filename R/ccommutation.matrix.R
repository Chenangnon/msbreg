# commutation matrix generation
# All codes from package 'matrixcalc' version 1.0-6 (Novomestky, 2022)

matrixcalcK.matrix <- function (r, c = r) {
  if (missing(r))
    stop("argument r is missing")
  if (!is.numeric(r))
    stop("argument r is not numeric")
  if (r != trunc(r))
    stop("argument r is not an integer")
  if (r < 2)
    stop("argument r is less than 2")
  if (!is.numeric(c))
    stop("argument c is not numeric")
  if (c != trunc(c))
    stop("argument c is not an integer")
  if (c < 2)
    stop("argument c is less than 2")
  H <- matrixcalcH.matrices(r, c)
  p <- r * c
  K <- matrix(0, nrow = p, ncol = p)
  for (i in 1:r) {
    for (j in 1:c) {
      Hij <- H[[i]][[j]]
      K <- K + (Hij %x% t(Hij))
    }
  }
  return(K)
}

matrixcalcH.matrices <- function (r, c = r) {
  if (missing(r))
    stop("argument r is missing")
  if (!is.numeric(r))
    stop("argument r is not numeric")
  if (r != trunc(r))
    stop("argument r is not an integer")
  if (r < 2)
    stop("argument r is less than 2")
  if (!is.numeric(c))
    stop("argument c is not numeric")
  if (c != trunc(c))
    stop("argument c is not an integer")
  if (c < 2)
    stop("argument c is less than 2")
  Ir <- diag(rep(1, r))
  Ic <- diag(rep(1, c))
  H <- list()
  for (i in 1:r) {
    H[[i]] <- list()
    for (j in 1:c) {
      H[[i]][[j]] <- Ir[i, ] %o% Ic[j, ]
    }
  }
  return(H)
}
