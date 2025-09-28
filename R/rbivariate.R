#' Bivariate Random Variates
#'
#' Sample from common bivariate distributions with mean vector \code{mu},
#' standard deviation vector \code{sigma} and correlation \code{rho}.
#'
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#' to be the required number.
#'
#' @param mu  numeric vector of two means.
#' @param sigma numeric vector of two standard deviations.
#' @param rho  numeric scalar in \code{(-1, 1)}, correlation.
#' Note that for log-normal variates, the range of feasible
#' correlation may be limited by the mean values in \code{mu},
#' especially negative correlations are generally limited to a
#' threshold way above \code{-1}.
#'
#' @param family character, the family of bivariate distributions to sample
#' from.
#'
#' Currently, one of \code{'unif'} for the uniform (the default),
#' \code{'norm'} for the normal (Gaussian), and \code{'lnorm'}
#' for the log-normal family of distributions.
#'
#' @param a,b numeric vectors, minima (\code{a}) and maxima (\code{b}).
#' Only used when \code{family = 'unif'} and both \code{a} and \code{b}
#' are specified: indicate the ranges of the values of the two uniform
#' variables: in that case, \code{mu} and \code{sigma} are ignored.
#'
#' @param verbose logical, should correlation adjustment for
#' log-normal distribution be reported?
#'
#' @param ... further arguments used by additional families.
#' Currently none is used.
#'
#' @details
#' Bivariate normal variates \eqn{(X_1, X_2)} are simulated by
#' first sampling \eqn{X_1} from a normal distribution, and then
#' sampling \eqn{X_2} conditional on \eqn{X_1}. Specifically, \eqn{X_1}
#' is sampled as:
#'
#' \code{x1 <- rnorm (n = n, mean = mu[1], sd = sigma[1])}.
#'
#' Then, the conditional mean of \eqn{X_2} given \eqn{X_1 =}\code{ x1} is
#' \code{mu_X_2.X_1 = mu[2] + rho * sigma[2] * (x1 - mu[1])/sigma[1]}
#' and the conditional standard deviation is
#' \code{sigma_X_2.X_1 = sqrt((1 - rho^2)) * sigma[2])}. Thus,
#' \eqn{X_2} is sampled as:
#'
#' \code{x2 <- rnorm (n = n, mean = mu_X_2.X_1, sd = sigma_X_2.X_1}.
#'
#' Bivariate uniform is simulated via \code{Beta(a,1)} after
#' \insertCite{demirtas2014generating;textual}{msbreg}.
#'
#' Bivariate log normal is simulated after
#' \insertCite{mielke1977covariance;textual}{msbreg}.
#'
#' @export rbivariate
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' ## Bivariate Uniform Variates
#' # Generate a matrix of observations
#' BIVU <- rbivariate (100000,
#'                     mu = c(0.5, 2),
#'                     sigma = c(2, 2),
#'                     rho = -0.7,
#'                     family = "unif")
#' # Sample mean vector
#' colMeans (BIVU)
#'
#' # Sample covariance matrix
#' cov (BIVU)
#'
#' # Sample correlation
#' cor (BIVU)[1,2]
#'
#' # Scatter plot
#' plot (BIVU)
#'
#' ## Bivariate Normal Variates
#' # Generate a matrix of observations
#' BIVN <- rbivariate (100000,
#'                     mu = c(0.5, 2),
#'                     sigma = c(2, 2),
#'                     rho = -0.7,
#'                     family = "norm")
#' # Sample mean vector
#' colMeans (BIVN)
#'
#' # Sample covariance matrix
#' cov (BIVN)
#'
#' # Sample correlation
#' cor (BIVN)[1,2]
#'
#' # Scatter plot
#' plot (BIVN)
#'
#' ## Bivariate Log normal Variates
#' # Generate a matrix of observations
#' BIVL <- rbivariate (100000,
#'                     mu = c(0.5, 2),
#'                     sigma = c(2, 2),
#'                     rho = 0.7,
#'                     family = "lnorm")
#' # Sample mean vector
#' colMeans (BIVL)
#'
#' # Sample covariance matrix
#' cov (BIVL)
#'
#' # Sample correlation
#' cor (BIVL)[1,2]
#'
#' # Scatter plot
#' plot (BIVL)
#'
rbivariate <- function(n,               # Sample size
                       mu = c(0, 0),    # Expectations
                       sigma = c(1, 1), # Standard deviations (truth)
                       rho = 0,         # Correlation (truth)
                       family = c("unif", "norm", "lnorm"),
                       a, b, verbose = TRUE, # For family = 'unif'; only used if both are supplied: mu and sigma are then ignored
                       ...) {
  family <- match.arg(family)

  out <- switch(family,
                unif = {
                  if (missing(a) | missing(b))
                    rbivunif (n = n, mu = mu, sigma = sigma, rho = rho)
                  else
                    rbivunif (n = n, rho = rho, a = a, b = b)
                },
                norm = {
                  rbivnorm (n = n, mu = mu, sigma = sigma, rho = rho)
                },
                lnorm = {
                  rbivlnorm (n = n, mu = mu, sigma = sigma, rho = rho, verbose = verbose)
                })

  structure(out, call = match.call())

}


# Simulate n samples of a bivariate normal vector
rbivnorm <- function (n,               # Sample size
                      mu = c(0, 0),    # Expectations
                      sigma = c(1, 1), # Standard deviations (truth)
                      rho = 0) {       # Correlation (truth)
  # Adjust/check arguments
  if(length(n) > 1)
    n <- length(n)
  tol <- .Machine$double.eps^0.5
  if (n <= 0 | abs(n - round(n)) > tol) {
    stop("'n' must be a positive integer.")
  }
  n <- max(ceiling(n), 1)
  mu <- rep(mu, length.out = 2)
  stopifnot(all(sigma >= 0))
  sigma <- rep(sigma, length.out = 2)
  rho <- rho[1]
  stopifnot(abs(rho) <= 1)

  # Generate normal samples
  x1 <- rnorm (n = n, mean = mu[1], sd = sigma[1])
  x2 <- rnorm (n = n, mean = mu[2] + rho * sigma[2] * (x1 - mu[1])/sigma[1],
               sd = sqrt((1 - rho^2)) * sigma[2])

  return(cbind(x1, x2))
}

# Simulate n samples of a bivariate log-normal vector
rbivlnorm <- function (n,               # Sample size
                       mu = c(1, 1),    # Expectations
                       sigma = c(1, 1), # Standard deviations (truth)
                       rho = 0,         # Correlation (truth)
                       verbose = TRUE) {# Print a warning when rho < -1 for the equivalent normal variates
  # Adjust/check arguments
  if(length(n) > 1)
    n <- length(n)
  tol <- .Machine$double.eps^0.5
  if (n <= 0 | abs(n - round(n)) > tol) {
    stop("'n' must be a positive integer.")
  }
  n <- max(ceiling(n), 1)
  mu <- rep(mu, length.out = 2)
  sigma <- rep(sigma, length.out = 2)
  stopifnot(all(mu > 0))
  stopifnot(all(sigma > 0))
  rho <- rho[1]
  stopifnot(abs(rho) <= 1)

  # From log-normal to normal parameters (mean and variance)
  varlog <- log(1 + (sigma/mu)^2)
  sdlog <- sqrt(varlog)
  meanlog <- log(mu) - 0.5 * varlog

  # Minimum/Maximum correlation value
  rho_min <- (exp(- prod(sdlog)) - 1) * prod(mu/sigma)
  rho_max <- (exp(  prod(sdlog)) - 1) * prod(mu/sigma)

  # Correlation at normal scale
  if (rho <= rho_min) {
    if (rho < rho_min & verbose) {
      warning(paste0("the required correlation 'rho' is less than the lower bound (" ,
                     round(rho_min, 4),
                     ") given the first two moments: using the lower bound"))
    }

    # Here rholog = -1
    x <- rbivnorm (n = n, mu = meanlog, sigma = sdlog, rho = -1)
  }
  else if (rho >= rho_max) {
    if (rho > rho_max & verbose) {
      warning(paste0("the required correlation 'rho' is greater than the upper bound (" ,
                     round(rho_max, 4),
                     ") given the first two moments: usng the upper bound"))
    }

    # Here rholog = 1
    x <- rbivnorm (n = n, mu = meanlog, sigma = sdlog, rho = 1)
  }
  else {
    rholog <- log(1 + rho * prod(sigma/mu)) / prod(sdlog)

    # Generate normal samples
    x <- rbivnorm (n = n, mu = meanlog, sigma = sdlog, rho = rholog)
  }

  # Transform to log-normal scale
  return(exp(x))
}

# Simulate n samples of a bivariate log-normal vector
# Code 'genbivunif.a' from package 'BivUnifBin'
rbivunif <- function (n,               # Sample size
                      mu = c(0, 0),    # Expectations
                      sigma = c(1, 1), # Standard deviations (truth)
                      rho = 0,         # Correlation (truth)
                      a, b) { # Only used if both are supplied: mu and sigma are then ignored
  # Adjust/check arguments
  if(length(n) > 1)
    n <- length(n)
  tol <- .Machine$double.eps^0.5
  if (n <= 0 | abs(n - round(n)) > tol) {
    stop("'n' must be a positive integer.")
  }
  n <- max(ceiling(n), 1)
  mu <- rep(mu, length.out = 2)
  sigma <- rep(sigma, length.out = 2)
  stopifnot(all(sigma > 0))
  rho <- rho[1]
  stopifnot(abs(rho) <= 1)

  # Simulate unit uniforms
  x1 <- runif(n)
  if (rho == 0) {
    x2 <- runif(n)
  }
  else if (abs(rho) == 1) {
    x2 <- sign(rho) * x1
  }
  else {
    v2 <- runif(n)
    shape <- (-5/2) + (1/2) * sqrt((rho + 49)/(rho + 1))
    u <- rbeta(n = n, shape1 = shape, shape2 = 1)
    x2 <- ifelse(v2 < 0.5, abs(u - x1), 1 - abs(1 - u - x1))
  }

  # Shit and re-scale if required
  if (missing(a) | missing(b)) {
    a <- mu - sqrt(3) * sigma
    b <- mu + sqrt(3) * sigma
  }
  else {
    a <- rep(a, length.out = 2)
    b <- rep(b, length.out = 2)
    stopifnot(all(a <= b))
  }

  x1 <- a[1] + (b[1] - a[1]) * x1
  x2 <- a[2] + (b[2] - a[2]) * x2

  return(cbind(x1, x2))
}


# Uniform distribution with centered parametrization
rrunif <- function (n, mu, sigma, a, b) {

  # Standard uniform
  x <- runif(n = n[1], min = 0, max = 1)

  # Shit and re-scale if required
  if (missing(a) | missing(b)) {
    a <- mu[1] - sqrt(3) * sigma[1]
    b <- mu[1] + sqrt(3) * sigma[1]
  }
  else {
    a <- a[1]
    b <- b[1]
    stopifnot(all(a <= b))
  }
  x <- a + (b - a) * x

  return(x)

}
rcunif <- rrunif

# Log normal distribution with centered parametrization
rrlnorm <- function (n,         # Sample size
                     mu = 1,    # Expectations
                     sigma = 1) { # Standard deviations (truth)
  # Adjust/check arguments
  if(length(n) > 1)
    n <- length(n)
  tol <- .Machine$double.eps^0.5
  if (n <= 0 | abs(n - round(n)) > tol) {
    stop("'n' must be a positive integer.")
  }
  n <- max(ceiling(n), 1)
  mu <- mu[1]
  sigma <- sigma[1]
  stopifnot(mu > 0)
  stopifnot(sigma > 0)

  # From log-normal to normal parameters (mean and variance)
  varlog <- log(1 + (sigma/mu)^2)
  sdlog <- sqrt(varlog)
  meanlog <- log(mu) - 0.5 * varlog

  # Normal samples
  x <- rnorm (n = n, mean = meanlog, sd = sdlog)

  # Transform to log-normal scale
  return(exp(x))

}
rclnorm <- rrlnorm
