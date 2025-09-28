# Core function to test for identifiability
.identifiability.x.y <- function (x, y, weights = 1, offset = 0,
                                  method = "iterative_rectifier",
                                  ...) {

  # Remove NAs from x, y, weights, and offset (if any)
  eval(drop.NAs.from.x.y.weights.offset())

  # Check matrix rank
  QRx <- qr(x, tol = 1e-08)
  xrank <- QRx$rank # matrix.rank (x)
  singular <- xrank < NCOL(x)
  xrank10 <- c(matrix.rank (x[y == 1, , drop = FALSE]),
               matrix.rank (x[y == 0, , drop = FALSE]))

  # Remove any aliased column in 'x'
  coeffs <- qr.coef(QRx, y) # y + c(1:nobs)
  keep <- !is.na(coeffs)
  if (!all(keep)) {
    if (any(keep)) {
      x <- x[, keep, drop = FALSE]
    }
    else {
      stop("should not get here (at least one column should remain): something went wrong")
    }
  }
  aliased <- which(!keep)

  # Check separation for 'y' as a binary or non-negative count
  separable <- test.sep (y = y, x = x,
                         weights = weights,
                         offset = offset,
                         binary = TRUE,
                         method = method[1],
                         ...)

  # Main output
  identifiable <- (!singular) & (!separable)

  # Attributes
  separated <- attr(separable, "separated")
  if (length(separated) & length(dnaaction)) {
    nseparated <- logical(nobs)
    nseparated[rnames %in% dnaaction] <- NA
    nseparated[!(rnames %in% dnaaction)][separated] <- TRUE
    separated <- which(nseparated)
  }

  ## REPORT WHEN y = constant after dropping separated observations?

  summary <- c(singular = singular,
               rank.all       = xrank,
               rank.successes = xrank10[1],
               rank.failures  = xrank10[2],
               separable = as.vector(separable)[1],
               attr(separable, "convergence"))
  range.resid <- attr(separable, "resid")
  epsilon <- attr(separable, "epsilon")

  # Final result
  return(structure(identifiable,
                   summary = summary,
                   aliased = if (length(aliased)) list(aliased),
                   separated = if (length(separated)) list(separated),
                   resid = range.resid,
                   epsilon = epsilon,
                   class = 'id.analysis'))
}
