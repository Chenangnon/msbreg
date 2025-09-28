#'
#' The Symmetric Truncated Cauchy Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the symmetric truncated Cauchy distribution with expectation \code{location}
#' and scale parameter \code{scale}.
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
#' @param scale vector of positive dispersion parameters (the larger
#' \code{scale} is, the larger the variance of the distribution).
#'
#' Setting \code{scale = } \eqn{\sqrt{\frac{\arctan(\tau)}{\tau - \arctan(\tau)}}} with
#' \eqn{\tau} \code{= tau} gives a unit variance distribution.
#'
#' The special value \code{scale = NULL} is allowed. This corresponds to
#' \code{scale } = \eqn{\pi\sqrt{\frac{1}{3} \frac{\arctan(\tau)}{\tau - \arctan(\tau)}}}
#' which produces a distribution with the same variance as the unit scale
#' logistic distribution (\eqn{\pi^2/ 3}).
# The corresponding distribution is concentrated around
# the location, but has long and heavy tails.
#'
#' @param log,log.p logical; if \code{TRUE}, densities and probabilities \code{p}
#' are given as \code{log(p)}.
#'
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#' \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
#'
#' @param tau positive real, truncation parameter \eqn{\tau} of the distribution.
#' For \code{location = }\eqn{\mu} and \code{scale = }\eqn{\tau}, only real
#' values in the range \eqn{[-\sigma \tau + \mu, \sigma\tau + \mu]} have
#' positive probability densities.
#'
#' The default value is \code{stats::qcauchy(.9999)} which truncates the
#' standard Cauchy distribution to the central region with \eqn{99.99\%}
#' probability mass (the approximate truncation region for the standard
#' distribution is \eqn{[-3183, 3183]}).
#'
# is system dependent, but generally evaluates to \code{1433540284805664.75}.
# This is meant to approach the non-truncated distribution as close as possible,
# while having finite moments (unlike the Cauchy distribution itself).
# Note that \code{.Machine$double.eps} is the smallest positive floating-point
# number \code{a} such that \code{1 + a != 1}. The default truncation
# parameter thus allows to capture \code{1 - 2 * .Machine$double.eps} of
# the probability under the Cauchy distribution (that is *almost* 100%).
#'
#' @export dtcauchy
#' @export ptcauchy
#' @export qtcauchy
#' @export rtcauchy
#'
#' @aliases ptcauchy
#' @aliases qtcauchy
#' @aliases rtcauchy
#' @aliases tCauchy
#'
# @usage
# dtcauchy (x, location = 0, scale = 1, log = FALSE,
#           tau = stats::qcauchy(.9999))
#
# ptcauchy (q, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE,
#           tau = stats::qcauchy(.9999))
#
# qtcauchy (p, location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE,
#           tau = stats::qcauchy(.9999))
#
# rtcauchy (n, location = 0, scale = 1,
#           tau = stats::qcauchy(.9999))
#
#' @details
#' The symmetric truncated Cauchy distribution with location parameter
#' \code{location = } \eqn{\mu}, scale parameter \code{scale = }\eqn{\sigma}
#' and truncation parameter \code{tau = } \eqn{\tau > 0}
#' has density function for
#' \eqn{x \in \left[-\sigma\tau + \mu, \sigma\tau+\mu \right]}
#' \insertCite{johnson1995continuous}{msbreg}
# vila2023closed
#'
#' \deqn{f(x) = \frac{1}{\pi\sigma P} \left[1 + \frac{(x - \mu)^2}{\sigma^2} \right]^{-1} ,}
#'
#' where \eqn{P = \frac{2}{\pi} \arctan(\tau)} is the probability
#' in the region \eqn{\left[-\tau, \tau\right]} under the standard Cauchy
#' density, and \eqn{f(x) = 0} otherwise. The distribution has expectation
#' \eqn{\mu} and variance \eqn{\sigma^2\left(\frac{\tau}{\arctan(\tau)} - 1\right)}.
#'
#' When \code{location} or \code{scale} are not specified they assume the default
#' values of \eqn{0} and \eqn{1} respectively, which correspond to zero mean and
#' variance \eqn{\frac{\tau}{\arctan(\tau)} - 1} (see
#' \insertCite{walck1996hand;nobrackets}{msbreg}, page 30).
#'
#' The logistic distribution is the basis for the widely use \code{logit} link
#' function. If using the truncated Cauchy distribution as an alternative
#' heavy tail link function, it can be useful to set the variance of link
#' to the variance as the unit scale logistic distribution (\eqn{\pi^2/ 3})
#' so as to make coefficients estimates comparable to logistic regression
#' coefficients.
#' Setting \code{scale = NULL} produces a truncated Cauchy distribution with
#' variance \eqn{\pi^2/ 3}. Note that the exact distribution depends on the
#' truncation parameter \code{tau}. For large \code{tau}, the distribution has
#' long but light tail.
# An alternative route to obtain a variance of \eqn{\pi^2/ 3} while keeping
# \code{scale = 1} is to set \code{tau} to approximately
# \code{6.0339587602477751105}.
# As such, \code{dtcauchy(x, scale = NULL, tau = 6.0339587602477751105)}
# is equivalent to \code{dtcauchy(x, scale = 1, tau = 6.0339587602477751105)}.
# This corresponds to a distribution with shorter but heavier tails than
# the logistic.
# \code{dtcauchy(x, scale = 1, tau = 9.23542058228320073)} has 99.99%
# probability in the same central region as the standard logistic.
#'
#' @return \code{dtcauchy} gives the density, \code{ptcauchy} gives the
#' distribution function, \code{qtcauchy} gives the quantile function,
#' and \code{rtcauchy} generates random variates.
#'
#' The length of the result is determined by \code{n} for \code{rtcauchy}, and
#' is the maximum of the lengths of the numerical arguments for the other
#' functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of
#' the result. Only the first elements of the logical arguments are used.
#'
#' Negative \code{scale} (\eqn{\sigma < 0}) is an error and returns
#' \code{NaN}. For \code{scale = 0}, the functions consider the limit
#' as \eqn{\sigma} approaches \code{0}, i.e. a point mass at
#' \eqn{\mu} (\code{location}).
#'
#' @references
#' \insertAllCited{}
#'
#' @source
#' \code{[dpq]tcauchy} are calculated using their relation to the standard
#' Cauchy distribution (see \link[stats]{Cauchy}).
#'
#' \code{rtcauchy} uses inversion with \code{qtcauchy}.
#'
#' @examples
#' ## Draw the truncated Cauchy density (in -10, 10) with various scales
#' Maxz <- 10
#' ylim <- max(dtcauchy(0, scale = 1, tau = Maxz),
#'             dtcauchy(0, scale = NULL, tau = Maxz))
#'
#' curve(dtcauchy(x, scale = NULL, tau = Maxz), -Maxz, Maxz, col = "red",
#'       n = 200, ylim = c(0, ylim), xlab = "x value",
#'       ylab = "Density at x", lty = 4)
#'
#' curve(dtcauchy(x, scale = 1, tau = Maxz), -Maxz, Maxz, col = "black",
#'       add = TRUE, n = 200, lty = 1)
#'
#' curve(dtcauchy(x, scale = 2, tau = Maxz), -Maxz, Maxz, col = "blue",
#'       add = TRUE, n = 200, lty = 2)
#'
#' curve(dtcauchy(x, scale = 5, tau = Maxz), -Maxz, Maxz, col = "brown3",
#'       add = TRUE, n = 200, lty = 3)
#'
#' legend(5, 0.32, legend = paste0("scale = ", c('NULL', 1, 2, 5)),
#'        col = c("red", "black", "blue", "brown3"), lty = c(4, 1:3))
#'
#' ## Draw the cumulative distribution function of the truncated Cauchy
#' #  compared with Logistic, Logistic Bridge, and Gaussian distributions
#' scale1 <- sqrt(atan(Maxz) / (Maxz - atan(Maxz))) # scale to get unit variance
#'
#' curve(ptcauchy(x, scale = scale1), -Maxz, Maxz, col = "red", n = 200)
#'
#' curve(plogis(x, scale = sqrt(3)/pi), -Maxz, Maxz, col = "brown3",
#'       add = TRUE, n = 200)
#'
#' curve(pbridge(x, scale = pi/sqrt(pi^2 + 3)), -Maxz, Maxz, col = "blue",
#'       add = TRUE, n = 200)
#'
#' curve(pnorm, -Maxz, Maxz, col = "black", ylim = c(0, 1), add = TRUE, n = 200,
#'      xlab = "x value", ylab = "P(X < x)")
#'
#' legend(2, 0.8, legend = c("tCauchy", "Logistic", "Bridge", "Gaussian"),
#'       col = c("red", "brown3", "blue", "black"), lty = c(1,1,1))
#'
#' # The truncated Cauchy has heavier tails.
#'
#' ## Random variate generation
#' # Add density curve to histogram
#' rtcvalues <- rtcauchy(1e+6, scale = 1, tau = Maxz)
#' hist(rtcvalues, freq = FALSE,
#'      ylim = c(0, dtcauchy(0, scale = 1, tau = Maxz)),
#'      xlim = c(-Maxz, Maxz), breaks = 200,
#'      xlab = "x value", ylab = "Density at x")
#'
#' curve(dtcauchy(x, scale = 1, tau = Maxz), -Maxz, Maxz, n = 200,
#'       add = TRUE, col = "red")
#'
#' # Variance
#' var(rtcvalues)
#' Maxz / atan(Maxz) - 1 # True variance
#'
# # Density comparisons
# curve(dtcauchy(x, tau = 6.0339587602477751105), -7, 7, ylab = "Density at x",
#       ylim = c(0, 0.35), n = 200, col = 'red')
#
# curve(dtcauchy(x, scale = 1, tau = 9.23542058228320073), -7, 7,
#      add = TRUE, col = 'blue', n = 200)
#
# curve(dlogis(x), -7, 7, ylim = c(0, 0.35), col = 'black', n = 200,
#      add = TRUE, lwd = 2)

#'

# Density function
dtcauchy <- function(x, location = 0, scale = 1, log = FALSE,
                     tau = stats::qcauchy(.9999)) { # -stats::qcauchy(.Machine$double.eps)

  # range = [-tau, tau] for the standard Cauchy before truncation

  # Check numerical arguments
  stopifnot(is.numeric(x), is.numeric(location), is.numeric(tau))
  tau <- abs(tau)
  stopifnot(all(tau > 0))
  if (is.null(scale)) {
    scale <- pi * sqrt((1 / 3) * atan(tau) / (tau - atan(tau)))
  }
  else {
    stopifnot(is.numeric(scale))
    NonValid <- scale < 0
    NonValid[is.na(NonValid)] <- FALSE
    scale[NonValid] <- NA
  }

  # Use the max of the dimensions of x, location, scale (numerical arguments) for the result
  finalres <- pickArgMaxLength(x, location, scale, tau)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(x),     length.out = n),
                          rep(c(location),  length.out = n),
                          rep(c(scale), length.out = n),
                          rep(c(tau),  length.out = n))
  NonValid <- is.na(rowSums(xlocationscale))
  finalres[NonValid] <- NA

  if (!all(NonValid)) {
    x <- (xlocationscale[!NonValid,1] - xlocationscale[!NonValid,2])
    scale <- xlocationscale[!NonValid,3]
    zeroscale <- scale == 0
    x[!zeroscale] <- x[!zeroscale] / scale[!zeroscale]
    tau <- xlocationscale[!NonValid,4]

    if (all(zeroscale)) {
      logpdf <- logb((x == 0) + 0)
    }
    else {
      logpdf <- zeroscale
      if (any(zeroscale))
        logpdf[zeroscale] <- logb((x[zeroscale] == 0) + 0)

      outrange <- (x < - tau) | (x > tau)
      if (any(outrange))
        logpdf[outrange] <- -Inf

      xok <- !zeroscale & !outrange

      if (any(xok)) {
        logpdf[xok] <- stats::dcauchy(x = x[xok], location = 0, scale = 1,
                                      log = TRUE)  -
          logb(1 - 2 * stats::pcauchy (q = tau[xok], location = 0, scale = 1,
                                       lower.tail = FALSE))

        logpdf[xok] <- logpdf[xok] - logb(scale[xok])
      }
    }

    finalres[!NonValid] <- if (log[1]) logpdf else exp(logpdf)
  }

  return(finalres)
}

# Distribution function
#' @rdname dtcauchy
ptcauchy <- function(q, location = 0, scale = 1,
                    lower.tail = TRUE, log.p = FALSE,
                    tau = stats::qcauchy(.9999)) {

  # range = [-tau, tau] for the standard Cauchy before truncation

  # Check numerical arguments
  stopifnot(is.numeric(q), is.numeric(location))
  stopifnot(is.numeric(tau))
  tau <- abs(tau)
  stopifnot(all(tau > 0))
  if (is.null(scale)) {
    scale <- pi * sqrt((1 / 3) * atan(tau) / (tau - atan(tau)))
  }
  else {
    stopifnot(is.numeric(scale))
    NonValid <- scale < 0
    NonValid[is.na(NonValid)] <- FALSE
    scale[NonValid] <- NA
  }

  # Use the max of the dimensions of q, location, scale (numerical arguments) for the result
  finalres <- pickArgMaxLength(q, location, scale, tau)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(q),     length.out = n),
                          rep(c(location),  length.out = n),
                          rep(c(scale), length.out = n),
                          rep(c(tau),  length.out = n))
  NonValid <- is.na(rowSums(xlocationscale))
  finalres[NonValid] <- NA

  if (!all(NonValid)) {
    q <- xlocationscale[!NonValid,1] - xlocationscale[!NonValid,2]
    scale <- xlocationscale[!NonValid,3]
    zeroscale <- scale == 0
    q[!zeroscale] <- q[!zeroscale] / scale[!zeroscale]
    tau <- xlocationscale[!NonValid,4]

    # cdf = P(X > x)
    if (all(zeroscale)) {
      cdf <- (q < 0) + 0
    }
    else {
      cdf <- zeroscale
      if (any(zeroscale))
        cdf[zeroscale] <- (q[zeroscale] < 0) + 0

      leftrange <- q < - tau
      if (any(leftrange))
        cdf[leftrange] <- 0

      rightrange <- q > tau
      if (any(rightrange))
        cdf[rightrange] <- 1

      xok <- !zeroscale & !leftrange & !rightrange

      if (any(xok)) {
        cdf[xok] <- stats::pcauchy(q = q[xok], location = 0, scale = 1,
                                   lower.tail = TRUE, log = FALSE)

        plower <- stats::pcauchy (q = tau[xok], location = 0, scale = 1,
                                  lower.tail = FALSE, log = FALSE)

        cdf[xok] <- (cdf[xok] - plower) / (1 - 2 * plower)
      }
    }

    if (!lower.tail[1])
      cdf <- 1 - cdf

    finalres[!NonValid] <- if (log.p[1]) logb(cdf) else cdf
  }

  return(finalres)
}

# Quantile function q
#' @rdname dtcauchy
qtcauchy <- function(p, location = 0, scale = 1,
                    lower.tail = TRUE, log.p = FALSE,
                    tau = stats::qcauchy(.9999)) {

  # Check numerical arguments
  stopifnot(is.numeric(p), is.numeric(location), is.numeric(tau))
  tau <- abs(tau)
  stopifnot(all(tau > 0))
  if (is.null(scale)) {
    scale <- pi * sqrt((1 / 3) * atan(tau) / (tau - atan(tau)))
  }
  else {
    stopifnot(is.numeric(scale))
    NonValid <- scale < 0
    NonValid[is.na(NonValid)] <- FALSE
    scale[NonValid] <- NA
  }

  # Use the max of the dimensions of p, location, scale (numerical arguments) for the result
  finalres <- pickArgMaxLength(p, location, scale, tau)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(p),     length.out = n),
                          rep(c(location),  length.out = n),
                          rep(c(scale), length.out = n),
                          rep(c(tau),  length.out = n))
  zerop <- xlocationscale[,1] == 0
  onep <- xlocationscale[,1] == 1
  NonValid0 <- is.na(rowSums(xlocationscale))
  NonValid <- NonValid | (zerop | onep)
  finalres[NonValid] <- NA

  if (!all(NonValid)) {
    p <- xlocationscale[!NonValid,1]
    if (log.p[1])
      p <- exp(p)
    if (!lower.tail[1])
      p <- 1 - p

    location <- xlocationscale[!NonValid,2]
    scale <- xlocationscale[!NonValid,3]
    zeroscale <- scale == 0
    tau <- xlocationscale[!NonValid,4]

    if (all(zeroscale)) {
      z <- 0
    }
    else {
      z <- zeroscale
      if (any(zeroscale))
        z[zeroscale] <- numeric(sum(zeroscale))

      plower <- finalres[!NonValid] + NA
      plower[!zeroscale] <- stats::pcauchy (q = tau[!zeroscale], location = 0,
                                            scale = 1, lower.tail = FALSE,
                                            log = FALSE)

      xok <- !zeroscale
      if (any(xok)) {
        cdf <- (1 - 2 * plower[xok]) * p[xok] + plower[xok]
        z[xok] <- stats::qcauchy(p = cdf, location = 0, scale = 1,
                                 lower.tail = TRUE, log = FALSE)

        z[xok] <- z[xok] * scale[xok]
      }
    }

    finalres[!NonValid] <- z + location
  }

  if (any(zerop & !NonValid0))
    finalres[zerop & !NonValid0] <- -tau

  if (any(onep & !NonValid0))
    finalres[onep & !NonValid0] <- tau

  return(finalres)
}

# Random variate generator

#' @rdname dtcauchy
rtcauchy <- function(n, location = 0, scale = 1,
                     tau = stats::qcauchy(.9999)) {
  # Check numerical arguments
  stopifnot(is.numeric(location), is.numeric(tau))
  tau <- abs(tau)
  stopifnot(all(tau > 0))
  if (is.null(scale)) {
    scale <- pi * sqrt((1 / 3) * atan(tau) / (tau - atan(tau)))
  }
  else {
    stopifnot(is.numeric(scale))
    NonValid <- scale < 0
    NonValid[is.na(NonValid)] <- FALSE
    scale[NonValid] <- NA
  }

  # Use length(n) for the result if length(n) > 1
  if (length(n) > 1)
    n <- length(n)

  # Re-cycle numerical vectors if required
  locationscale <- cbind(rep(c(location),  length.out = n),
                         rep(c(scale),     length.out = n),
                         rep(c(tau),      length.out = n))

  NonValid <- is.na(rowSums(locationscale))
  x <- numeric(n)
  x[NonValid] <- NA

  # Use inversion method
  if (any(!NonValid))
    x[!NonValid] <- qtcauchy(runif(sum(!NonValid)),
                             location = locationscale[!NonValid,1],
                             scale    = locationscale[!NonValid,2],
                             tau     = locationscale[!NonValid,3])

  return(x)
}
