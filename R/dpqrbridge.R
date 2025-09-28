#'
#' The Logistic Bridge Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the logistic bridge distribution with mean equal to \code{location} and scale
#' parameter equal to \code{scale}.
#'
#' @param x,q vector of quantiles.
#'
#' @param p vector of probabilities.
#'
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#' to be the number required (\code{n = length(n)}).
#'
#' @param location vector of means.
#'
#' @param scale vector of precision parameters in the range \code{(0, 1]} (the
#' larger \code{scale} is, the lower the variance).
#'
#' @param log,log.p logical; if \code{TRUE}, densities and probabilities \code{p}
#' are given as \code{log(p)}.
#'
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#' \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
#'
#' @export dbridge
#' @export pbridge
#' @export qbridge
#' @export rbridge
#'
#' @aliases pbridge
#' @aliases qbridge
#' @aliases rbridge
#' @aliases bridge
#'
#' @usage
#' dbridge (x, location = 0, scale = 1/sqrt(2), log = FALSE)
#'
#' pbridge (q, location = 0, scale = 1/sqrt(2), lower.tail = TRUE, log.p = FALSE)
#'
#' qbridge (p, location = 0, scale = 1/sqrt(2), lower.tail = TRUE, log.p = FALSE)
#'
#' rbridge (n, location = 0, scale = 1/sqrt(2))
#'
#' @details
#' The logistic bridge distribution with \code{scale = }\eqn{\varphi} developed
#' by \insertCite{wang2003matching;textual}{msbreg} and extended to include a
#' location parameter \code{location = } \eqn{\mu} has density function
#'
#' \deqn{f(x) = \frac{1}{2\pi} \frac{\sin (\varphi\pi)}{\cosh (\varphi (x - \mu)) + \cos (\varphi\pi)},}
#'
#' distribution function
#'
#' \deqn{F(x) = 1 - \frac{1}{\pi\varphi} \left[\frac{\pi}{2} - \arctan{\left( \frac{\exp(\varphi (x - \mu)) + \cos (\varphi\pi)}{\sin (\varphi\pi)} \right) } \right],}
#'
#' and quantile function
#'
#' \deqn{F^{-1}(p) = \frac{1}{\varphi} \log \left[ \frac{\sin (\varphi\pi p)}{\sin (\varphi \pi (1 - p))} \right] + \mu.}
#'
#' It is a long-tailed distribution (as compared to the Gaussian) with mean
#' \eqn{\mu} and variance \eqn{\pi^2 \left(\varphi^{-2} - 1\right) / 3} \insertCite{wang2003matching}{msbreg}.
#'
#' When \code{location} or \code{scale} are not specified they assume the default
#' values of \eqn{0} and \eqn{1/\sqrt{2}}, respectively, which correspond to
#' zero mean and the same variance as the unit scale logistic distribution
#' (\eqn{\pi^2/ 3}).
#'
#' @return \code{dbridge} gives the density, \code{pbridge} gives the distribution
#' function, \code{qbridge} gives the quantile function, and \code{rbridge}
#' generates random variates.
#'
#' The length of the result is determined by \code{n} for \code{rbridge}, and
#' is the maximum of the lengths of the numerical arguments for the other
#' functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result.
#' Only the first elements of the logical arguments are used.
#'
#' Negative or null \code{scale} (\eqn{\varphi \leq 0}) is an error and returns
#' \code{NaN}. The same holds for \code{scale > 1}. For \code{scale = 1}, the
#' functions consider the limit as \eqn{\varphi} approaches \code{1}, i.e. a
#' point mass at \eqn{\mu} (\code{location}).
#'
#' @references
#' \insertAllCited{}
#'
#' @source
#' \code{[dpq]bridge} are calculated directly from the definitions.
#'
#' \code{rbridge} uses inversion (see \eqn{F^{-1}(p)} in the section \code{'Details'}).
#'
#' @examples
#' ## Draw the Bridge density with various scales
#' curve(dbridge(x, scale = 0.25), -10, 10, col = "black", n = 200,
#'       lty = 1, ylim = c(0, 0.4), xlab = "x value", ylab = "Density at x")
#'
#' curve(dbridge(x, scale = 0.5), -10, 10, col = "black", add = TRUE, n = 200,
#'       lty = 2)
#'
#' curve(dbridge(x, scale = 0.75), -10, 10, col = "black", add = TRUE, n = 200,
#'       lty = 3)
#'
#' legend(5, 0.4, legend = paste0("scale = ", c(0.25, 0.50, 0.75)),
#'       lty = 1:3)
#'
#' ## Draw the cumulative distribution functions of
#' # Gaussian, Logistic and Bridge distributions
#' curve(pnorm, -6, 6, col = "black", ylim = c(0, 1), n = 200,
#'      xlab = "x value", ylab = "P(X < x)")
#'
#' curve(plogis(x, scale = sqrt(3)/pi), -6, 6, col = "blue",
#'       add = TRUE, n = 200)
#'
#' curve(pbridge(x, scale = pi/sqrt(pi^2 + 3)), -6, 6, col = "red",
#'       add = TRUE, n = 200)
#'
#' legend(2, 0.8, legend = c("Gaussian", "Logistic", "Bridge"),
#'       col = c("black", "blue", "red"), lty = c(1,1,1))
#'
#'
#' ## Random variate generation
#' # Add density curve to histogram
#' hist(rbridge(500000, scale = 0.6), freq = FALSE,
#'      ylim = c(0, 0.25), xlim = c(-10, 10),
#'      breaks = 200, xlab = "x value", ylab = "Density at x")
#'
#' curve(dbridge(x, scale = 0.6), -10, 10, n = 200, add = TRUE, col = "blue")
#'
#' # Variance
#' var(rbridge(500000, 0, scale = 1/sqrt(2)))
#' pi^2/3 # True variance

# Density function
dbridge <- function(x, location = 0, scale = 1/sqrt(2), log = FALSE) {

  # Check numerical arguments
  stopifnot(is.numeric(x), is.numeric(location), is.numeric(scale))
  NonValid <- (scale <= 0) | (scale > 1)
  NonValid[is.na(NonValid)] <- FALSE
  scale[NonValid] <- NA

  # Use the max of the dimensions of x, location, scale (numerical arguments) for the result
  finalres <- pickArgMaxLength(x, location, scale)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(x),     length.out = n),
                      rep(c(location),  length.out = n),
                      rep(c(scale), length.out = n))
  NonValid <- is.na(rowSums(xlocationscale))
  finalres[NonValid] <- NA

  if (!all(NonValid)) {
    x <- (xlocationscale[!NonValid,1] - xlocationscale[!NonValid,2])
    scale <- xlocationscale[!NonValid,3]
    unitphi <- scale == 1

    if (all(unitphi)) {
      logpdf <- logb((x == 0) + 0)
    }
    else {
      logpdf <- unitphi
      if (any(unitphi))
        logpdf[unitphi] <- logb((x[unitphi] == 0) + 0)

      logpdf[!unitphi] <- - (logb(2) + 2*lgamma(.5)) + logb(sinpi(scale[!unitphi])) -
        logb(cosh(scale[!unitphi] * x[!unitphi]) + cospi(scale[!unitphi]))
    }

    finalres[!NonValid] <- if (log[1]) logpdf else exp(logpdf)
  }

  return(finalres)
}

# Distribution function
pbridge <- function(q, location = 0, scale = 1/sqrt(2),
                    lower.tail = TRUE, log.p = FALSE) {

  # Check numerical arguments
  stopifnot(is.numeric(q), is.numeric(location), is.numeric(scale))
  NonValid <- (scale <= 0) | (scale > 1)
  NonValid[is.na(NonValid)] <- FALSE
  scale[NonValid] <- NA

  # Use the max of the dimensions of q, location, scale (numerical arguments) for the result
  finalres <- pickArgMaxLength(q, location, scale)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(q),     length.out = n),
                          rep(c(location),  length.out = n),
                          rep(c(scale), length.out = n))
  NonValid <- is.na(rowSums(xlocationscale))
  finalres[NonValid] <- NA

  if (!all(NonValid)) {
    q <- xlocationscale[!NonValid,1] - xlocationscale[!NonValid,2]
    scale <- xlocationscale[!NonValid,3]
    unitphi <- scale == 1

    # cdf = P(X > x)
    if (all(unitphi)) {
      cdf <- (q < 0) + 0
    }
    else {
      cdf <- unitphi
      if (any(unitphi))
        cdf[unitphi] <- (q[unitphi] < 0) + 0

      cdf[!unitphi] <- ((pi/2) -
                             atan((exp(scale[!unitphi] * q[!unitphi]) +
                                     cospi(scale[!unitphi]))/
                                    sinpi(scale[!unitphi]))) /
        (scale[!unitphi] * pi)
    }

    if (lower.tail[1])
      cdf <- 1 - cdf

    finalres[!NonValid] <- if (log.p[1]) logb(cdf) else cdf
  }

  return(finalres)
}

# Quantile function
qbridge <- function(p, location = 0, scale = 1/sqrt(2),
                    lower.tail = TRUE, log.p = FALSE) {

  # Check numerical arguments
  stopifnot(is.numeric(p), is.numeric(location), is.numeric(scale))
  NonValid <- (scale <= 0) | (scale > 1)
  NonValid[is.na(NonValid)] <- FALSE
  scale[NonValid] <- NA

  # Use the max of the dimensions of p, location, scale (numerical arguments) for the result
  finalres <- pickArgMaxLength(p, location, scale)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(p),     length.out = n),
                          rep(c(location),  length.out = n),
                          rep(c(scale), length.out = n))
  zerop <- xlocationscale[,1] == 0
  onep <- xlocationscale[,1] == 1
  NonValid0 <- is.na(rowSums(xlocationscale))
  NonValid <- NonValid0 | (zerop | onep)
  finalres[NonValid] <- NA

  if (!all(NonValid)) {
    p <- xlocationscale[!NonValid,1]
    if (log.p[1])
      p <- exp(p)
    if (!lower.tail[1])
      p <- 1 - p

    location <- xlocationscale[!NonValid,2]
    scale <- xlocationscale[!NonValid,3]
    unitphi <- scale == 1

    if (all(unitphi)) {
      z <- 0
    }
    else {
      z <- unitphi
      if (any(unitphi))
        z[unitphi] <- numeric(sum(unitphi))

      z[!unitphi] <- logb(sinpi(scale[!unitphi] * p[!unitphi]) /
                              sinpi(scale[!unitphi] * (1 - p[!unitphi]))) / (scale[!unitphi])
    }

    finalres[!NonValid] <- z + location
  }

  if (any(zerop & !NonValid0))
    finalres[zerop & !NonValid0] <- -Inf

  if (any(onep & !NonValid0))
    finalres[onep & !NonValid0] <- Inf

  return(finalres)
}

# Random variate generator
rbridge <- function(n, location = 0, scale = 1/sqrt(2)) {
  # Check numerical arguments
  stopifnot(is.numeric(location), is.numeric(scale))
  NonValid <- (scale <= 0) | (scale > 1)
  NonValid[is.na(NonValid)] <- FALSE
  scale[NonValid] <- NA

  # Use length(n) for the result if length(n) > 1
  if (length(n) > 1)
    n <- length(n)

  # Re-cycle numerical vectors if required
  locationscale <- cbind(rep(c(location),  length.out = n),
                          rep(c(scale),     length.out = n))

  NonValid <- is.na(rowSums(locationscale))
  x <- numeric(n)
  x[NonValid] <- NA

  # Use inversion method
  if (any(!NonValid))
    x[!NonValid] <- qbridge(runif(sum(!NonValid)),
                            location = locationscale[!NonValid,1],
                            scale    = locationscale[!NonValid,2])

  return(x)
}

pickArgMaxLength <- function(...) {
  Args <- list(...)
  Args[[which.max(sapply(Args, FUN = length))[1]]]
}
