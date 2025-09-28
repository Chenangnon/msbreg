# TO BE COMPLETED with density, quantile and random deviate generator
#
#' The normal location mixture of logistic distributions
#'
#' The normal location mixture of logistic distribution is the mixed logistic
#' distribution obtained from continuously mixing logit probabilities using
#' a normal mixing distribution.
#'
#' @param x description
#'
#' @param location description
#'
#' @param scale description
#'
#' @param mix.sd description
#'
#' @param lower.tail description
#'
#' @param log.p description
#'
#' @param N description
#'
#' @details
#' When \code{mix.sd = 0}, \code{pmixlogis} is the common logistic cumulative
#' distribution function (cdf) \link{plogis}. For \code{mix.sd > 0}, the function
#' uses one of two methods: either an adaptive numerical integration with
#' \link{integrate}, or a truncated series approximation to the cdf. Numerical
#' integration is used when \code{N = NULL} or \code{N = 0} (the default).
#' When a numerical value is supplied for \code{N}, this value is used as the
#' truncation point of the series representation to the cdf.
#'
#' @export pmixlogis
#'
pmixlogis <- function (x, location = 0, scale = 1, mix.sd = 0,
                    lower.tail = TRUE, log.p = FALSE, N = NULL) {
  # Pre-processing
  x0 <- x
  x <- c(x)
  x <- cbind(c(x), c(location), c(scale), c(mix.sd))
  location <- x[,2]
  scale <- x[,3]
  mix.sd <- x[,4]
  x <- x[,1]

  # Initialized a vector to store the results
  p <- numeric(length(x))

  # Special case of a mass point mixing distribution
  zeromsigma <- mix.sd == 0
  if (any(zeromsigma)) {
    p[zeromsigma] <- plogis (x[zeromsigma], location = location[zeromsigma],
                             scale = scale[zeromsigma], lower.tail = lower.tail,
                             log.p = log.p)
    if (all(zeromsigma)) {
      x0[] <- p
      return(x0)
    }
  }

  # Scale x
  x <- (x[!zeromsigma] - location[!zeromsigma]) / scale[!zeromsigma]
  mix.sd <- mix.sd[!zeromsigma]

  # Use numerical integration if required
  if (is.null(N)) {
    p[!zeromsigma] <- apply(cbind(x, mix.sd), MARGIN = 1, FUN = pmlogit)
  }
  else {
    if (N == 0) {
      p[!zeromsigma] <- apply(cbind(x, mix.sd), MARGIN = 1, FUN = pmlogit)
    }
    else {
      # Set the truncation born for the series representation
      N <- max(30, N)

      # Get the cdf values for each individual point in x
      p[!zeromsigma] <- apply(cbind(x, mix.sd), MARGIN = 1, FUN = pmlogit.series, N = N)
    }
  }

  # post processing
  if (!lower.tail[1])
    p[!zeromsigma] <- 1 - p[!zeromsigma]
  if (log.p)
    p[!zeromsigma] <- log(p[!zeromsigma])

  x0[] <- p
  return(x0)
}

# Internal routine for 'pmixlogis'
# 'z' is a vector of two element, the first is a real, the second is assumed to be a non negative real
pmlogit.series <- function(z, N = 100) {
  eta <- z[1]

  # If eta = 0, pmlogit = 0.5 as for the logistic distribution
  if (eta == 0) {
    return(.5)
  }

  # If a mass point mixing distribution
  z <- z[2]
  if (z <= 0) {
    return(plogis (eta, location = 0, scale = 1,
                   lower.tail = TRUE, log.p = FALSE))
  }

  # Series only used for eta != 0
  N <- ceiling(abs(N))
  N <- max(N, 30)
  k <- 0:N
  p <- sum(
    choose(-1, k) * (exp((k+1) * eta + .5 * ((k+1)*z)^2 + pnorm(-eta/z - (k+1)*z, mean = 0, sd = 1,
                                                                lower.tail = TRUE, log.p = TRUE)) +
                       exp(-k * eta + .5 * (k*z)^2 + pnorm(eta/z - k*z, mean = 0, sd = 1,
                                                           lower.tail = TRUE, log.p = TRUE)))
  )

  return(p)
}

pmlogit <- function(z, abs.tol = 1e-15, N = 10000) {
  eta <- z[1]

  # If eta = 0, pmlogit = 0.5 as for the logistic distribution
  if (eta == 0) {
    return(.5)
  }

  # If a mass point mixing distribution
  z <- z[2]
  if (z <= 0) {
    return(plogis (eta, location = 0, scale = 1,
                   lower.tail = TRUE, log.p = FALSE))
  }

  # If eta == 0, return 0.5

  # If eta > 0, use -eta, and return 1 - p

  p <- catch.conditions({
    integrate (f = function(x) {
      plogis (x, location = 0, scale = 1,
              lower.tail = TRUE, log.p = FALSE) *
        dnorm(x, mean = eta, sd = z)
    },
    lower = -Inf, upper = Inf, abs.tol = abs.tol)$value
  })$value

  if (any(class(p) %in% c("simpleError", "error", "condition"))) {
    p <- catch.conditions({
      integrate (f = function(u) {
        plogis (sd * u, location = 0, scale = 1,
                lower.tail = TRUE, log.p = FALSE) *
          dnorm(u, mean = eta/sd)
      },
      lower = -Inf, upper = Inf, abs.tol = abs.tol)$value
    })$value
  }

  if (any(class(p) %in% c("simpleError", "error", "condition"))) {
    if (identical(p$message, " oups ")) {
      if (abs.tol < 1e-8) { # (.Machine$double.eps^0.25)^2
        p <- catch.conditions({
          integrate (f = function(x) {
            plogis (x, location = 0, scale = 1,
                    lower.tail = TRUE, log.p = FALSE) *
              dnorm(x, mean = eta, sd = z)
          },
          lower = -Inf, upper = Inf, abs.tol = 1e-8)$value
        })$value
        if (any(class(p) %in% c("simpleError", "error", "condition"))) {
          p <- catch.conditions({
            integrate (f = function(x) {
              plogis (x, location = 0, scale = 1,
                      lower.tail = TRUE, log.p = FALSE) *
                dnorm(x, mean = eta, sd = z)
            },
            lower = -Inf, upper = Inf, abs.tol = .Machine$double.eps^0.25)$value
          })$value
        }
        if (any(class(p) %in% c("simpleError", "error", "condition"))) {
          p <- pmlogit.series (c(eta, z), N = N)
        }
      }
    }
    else {
      p <- pmlogit.series (c(eta, z), N = N)
    }
  }

  return(p)
}

