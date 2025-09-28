
#' Matrix of predictors
#'
#' Simulate a matrix of numerical predictors from uniform, normal or
#' shifted log-normal variates.
#'
#' @param n description
#'
#' @param xdistr character vector, indicates the distributions to sample from.
#'
#' @param xmin,xmax numeric vectors. These indicate the minimum/maximum of the
#' uniform distribution to sample for \code{xdistr = 'unif'}. For the other
#' two distributions, the samples \code{x} are re-scaled to
#' \code{xmin + (x + 2) * (xmax - xmin) / 10} for \code{xdistr = 'norm'} and
#' to \code{xmin + x * (xmax - xmin) / 10} for \code{xdistr = 'lnorm'}.
#'
#' @param cor.structure character, specifies a correlation structure for
#' the sampled random vector. Defaults to \code{cor.structure = NULL} for
#' independent components when the dimension of the sampled random vector is
#' two or more.
# Currently, the only available alternative is \code{'Gaussian'}:
# the centered version of the random vector is multiplied by the cholesky
# matrix \code{cor.par}.
#'
# @param cor.par description
#'
#' @export sim.predictors
#'
sim.predictors <- function (n = 100, xdistr = 'unif',
                            xmin = -2, xmax = 8,
                            cor.structure = NULL) {
  p <- max(sapply(list(xdistr, xmin, xmax), FUN = length))
  xdistr <- rep(xdistr, length.out = p)
  xmin <- rep(xmin, length.out = p)
  xmax <- rep(xmax, length.out = p)

  if (!is.null(cor.structure) & p > 1) {
    warning("ignoring 'cor.structure'")
  }

  if (p == 1) {
    x <- cbind(sim.onepredictor (n = n, xdistr = xdistr[1],
                                 xmin = xmin, xmax = xmax))
  }
  else {
    x <- sapply(1:p,
                FUN = function(j) {
                  c(sim.onepredictor (n = n, xdistr = xdistr[j],
                                    xmin = xmin[j], xmax = xmax[j]))
                })
  }

  colnames(x) <- paste0('x', 1:p)

  return(x)
}

sim.onepredictor <- function (n, xdistr, xmin = -2, xmax = 8) {
  switch(tolower(xdistr),
         unif = {
           x <- runif(n, min = xmin, max = xmax)
         },
         norm = {
           x <- rnorm(n, mean = 3, sd = 1.6)
           x <- xmin + (x + 2) * (xmax - xmin) / 10 # Rescale to the interval [xmax, xmin]
         },
         lnorm = {
           x <- rlnorm(n, meanlog = 1.1, sdlog = 0.5)
           x <- xmin + (x - 0) * (xmax - xmin) / 10 # Rescale to the interval [xmax, xmin]
         })

  return (x)
}
