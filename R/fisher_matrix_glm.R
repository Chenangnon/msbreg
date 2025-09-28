#' @rdname fisher.matrix
#' @exportS3Method msbreg::fisher.matrix
fisher.matrix.lm <- function (object,
                              sigma = NULL, ...) {

  # Recover the design matrix from the lm fit (assumes that the fit is not empty!)
  Qr <- object$qr %||% # Code from: stats:::qr.lm
    stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
  x <- qr.X (Qr, complete = FALSE)

  if (is.null(sigma)) {
    weights <- object$weights
    if (is.null(weights))
      weights <- 1
    rss <- sum(weights * object$residuals^2)
    sigma <- sqrt(rss / object$df.residual)
  }
  else {
    stopifnot(is.numeric(sigma))
    stopifnot(sigma > 0)
  }

  Ematrix <- t(x) %*% x
  Ematrix <- Ematrix / (sigma^2)

  return(structure(Ematrix,
                   sigma = sigma,
                   class = c("fisher.matrix",
                             class(Ematrix))))
}
### setMethod commented out to prevent the warning:
### Undocumented S4 methods:
###      generic 'fisher.matrix' and siglist 'lm'
#setMethod ("fisher.matrix",
#           signature = "lm",
#           definition = fisher.matrix.lm)

#' @rdname fisher.matrix
#' @exportS3Method msbreg::fisher.matrix
fisher.matrix.glm <- function (object,
                               sigma = NULL, ...) {

  # Recover the design matrix from the glm fit (assumes that the fit is not empty!)
  Qr <- object$qr %||% # Code from: stats:::qr.lm
    stop("lm object does not have a proper 'qr' component.\n Rank zero or should not have used lm(.., qr=FALSE).")
  x <- qr.X (Qr, complete = FALSE)

  fam <- object$family
  if (is.null(sigma)) {
    df.r <- object$df.residual
    sigma <- if (!is.null(fam$dispersion) && !is.na(fam$dispersion))
      fam$dispersion
    else if (fam$family %in% c("poisson", "binomial"))
      1
    else if (df.r > 0) {
      if (any(object$weights == 0))
        warning("observations with zero weight not used for calculating dispersion")
      sum((object$weights * object$residuals^2)[object$weights >
                                                  0])/df.r
    }
    else {
      NaN
    }
  }
  else {
    stopifnot(is.numeric(sigma))
    stopifnot(sigma > 0)

    if (fam$family %in% c("poisson", "binomial")) {
      warning(paste0("non-unit dispersion supplied for the ",
                     fam$family, " family"))
    }
  }

  Ematrix <- t(x) %*% x
  Ematrix <- Ematrix / (sigma^2)

  return(structure(Ematrix,
                   sigma = sigma,
                   class = c("fisher.matrix",
                             class(Ematrix))))
}
### setMethod commented out to prevent the warning:
### Undocumented S4 methods:
###      generic 'fisher.matrix' and siglist 'glm'
#setMethod ("fisher.matrix",
#           signature = "glm",
#           definition = fisher.matrix.glm)
