
# Get Linear Predictor
#
# Compute linear predictor from a vector of regression parameters and
# a matrix of model terms.
#
#

get.eta <- function (beta, design.matrix,
                     intercept = attr(design.matrix, 'intercept'),
                     offset = attr(design.matrix, 'offset')) {
  if (is.null(intercept))
    intercept <- 0
  stopifnot(length(beta) == NCOL(design.matrix) + intercept)

  if (intercept) {
    eta <- beta[1] + c(design.matrix %*% beta[-1])
  }
  else {
    eta <- c(design.matrix %*% beta)
  }

  if (is.null(offset))
    offset <- 0

  return(eta + offset)
}
