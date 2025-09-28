#'
#' The Asymmetric Normal Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the asymmetric normal distribution with location equal to \code{location},
#' scale parameters equal to \code{sigma} (left) and \code{tau} (right).
#'
#' @param x,q vector of quantiles.
#'
#' @param p vector of probabilities.
#'
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#' to be the number required (\code{n = length(n)}).
#'
#' @param location vector of location parameters.
#'
#' @param scale1,scale2 vectors of positive dispersion parameters (the larger
#' each scale parameter is, the larger the variance).
#'
#' @param log,log.p logical; if \code{TRUE}, densities and probabilities \code{p}
#' are given as \code{log(p)}.
#'
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#' \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
#'
#' @export danorm
#' @export panorm
#' @export qanorm
#' @export ranorm
#'
#' @aliases panorm
#' @aliases qanorm
#' @aliases ranorm
#' @aliases anorm
#'
#' @details
#' The asymmetric normal distribution \eqn{AN(m, \sigma, \tau)} was
#' proposed by \insertCite{sicard2013loinormale;textual}{msbreg}. For
#' \eqn{m} \code{ = location}, \eqn{\sigma} \code{= scale1} and
#' \eqn{\tau} \code{= scale2}, it is the normal distribution with standard
#' deviation \eqn{\sigma} for quantiles below the point \eqn{m}, and standard
#' deviation \eqn{\tau} for quantiles above \eqn{m}. For \eqn{\tau = \sigma},
#' the \eqn{AN(m, \sigma, \tau)} distribution is thus reduced to the Gaussian
#' distribution with mean \eqn{m} and standard deviation \eqn{\sigma} which is
#' symmetric. For \eqn{\sigma \neq \tau}, \eqn{AN(m, \sigma, \tau)} is a skewed
#' distribution with negative skewness (longer left tail) when
#' \eqn{\sigma > \tau} and positive skewness (longer right tail) where
#' \eqn{\sigma < \tau}. The distribution is also known as
#' the epsilon–skew–normal distribution
#' \insertCite{mudholkar2000epsilon}{msbreg}.
#'
#' The asymmetric normal distribution has density function
#'  \insertCite{sicard2013loinormale}{msbreg}
#'
#' \deqn{f(x) = \frac{\sqrt{2}}{\sqrt{\pi} (\sigma + \tau)} \exp \left( -\frac{(x - m)^2}{2 \sigma_x^2} \right),}
#'
#' where \eqn{\sigma_x = \sigma} if \eqn{x \leq m} and \eqn{\sigma_x = \tau}
#' if \eqn{x > m}. For \eqn{X \sim AN(m, \sigma, \tau)}, the probability of
#' the event \eqn{X \leq m} is
#' \eqn{\alpha = P\left[X \leq m \right] = \frac{\sigma}{\sigma + \tau}}.
#' Setting \eqn{\omega = \sigma + \tau} so that \eqn{\sigma = \alpha \omega},
#' \eqn{\tau = (1 - \alpha) \omega}, we can rewrite the density function as
#'
#' \deqn{f(x) = \left\lbrace \begin{matrix} \quad \quad 2 \alpha \phi \left(x|m, \alpha \omega\right) & \mathrm{ if } \quad x \leq m \\  \\ 2 (1 - \alpha)  \phi \left(x|m, (1 - \alpha) \omega\right)  & \mathrm{ if } \quad x > m \end{matrix} \right.}
#'
#' where \eqn{\phi(\cdot|m,s)} denotes the density of the normal distribution
#' with mean \eqn{m} and standard deviation \eqn{s}.
#' The cumulative distribution function (cdf) \eqn{F(x)=P\left[X \leq x \right]}
#' is given by
#'
#' \deqn{F(x) = \left\lbrace \begin{matrix} \quad \quad 2 \alpha \Phi \left( \frac{x - m}{\sigma} \right) & \mathrm{ if } \quad x \leq m \\  \\ \alpha + (1 - \alpha) \left[ 2 \Phi \left( \frac{x - m}{\tau} \right) - 1 \right]  & \mathrm{ if } \quad x > m \end{matrix} \right.}
#'
#' where \eqn{\Phi} denotes the cdf of the standard normal distribution.
#' The quantile function of \eqn{X} is for \eqn{p \in (0, 1)}:
#'
#' \deqn{F^{-1}(p) = m + \sigma_x \Phi^{-1}(w)}
#'
#' where \eqn{\left\lbrace \begin{matrix} \sigma_x = \sigma & \mathrm{ and } \quad w = \frac{p}{2 \alpha} & \mathrm{ if } \quad p \leq \alpha, \\  \\ \sigma_x = \tau & \mathrm{ and } \quad w = \frac{p - \alpha}{2(1 - \alpha)} + \frac{1}{2}  & \mathrm{ if } \quad p > \alpha. \end{matrix} \right.}
#'
#' The \eqn{AN(m, \sigma, \tau)} distribution has expectation
#' \eqn{m+\sqrt{\frac{2}{\pi}} \left(\tau-\sigma\right)} and variance
#' \eqn{\sigma\tau + (1 - \frac{2}{\pi}) \left(\tau-\sigma\right)^2}.
#' More properties of a variable \eqn{X \sim AN(m, \sigma, \tau)} can
#' be derived from the following stochastic representation:
#'
#' \deqn{X \overset{d}{=} -U |X_{\sigma}| + (1 - U) |X_{\tau}| + m}
#'
#' where \eqn{\overset{d}{=}} means "is equal in distribution to",
#' \eqn{U} is a Bernoulli variable with success probability
#' \eqn{\alpha = \frac{\sigma}{\sigma + \tau}}, \eqn{|X_{a}|} is, for
#' \eqn{a \in \left\lbrace \sigma, \tau \right\rbrace}, the absolute value
#' of a Gaussian variable with mean zero and standard deviation \eqn{a},
#' and \eqn{U}, \eqn{X_{\sigma}} and \eqn{X_{\tau}} are mutually
#' independent \insertCite{sicard2013loinormale}{msbreg}.
#'
#' When \code{location} or \code{scale1} or \code{scale2} are not specified
#' they respectively assume the default values of \eqn{0}, \eqn{1} and \eqn{1},
#' which correspond to the standard normal distribution.
#'
#' @return \code{danorm} gives the density, \code{panorm} gives the distribution
#' function, \code{qanorm} gives the quantile function, and \code{ranorm}
#' generates random variates.
#'
#' The length of the result is determined by \code{n} for \code{ranorm}, and
#' is the maximum of the lengths of the numerical arguments for the other
#' functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result.
#' Only the first elements of the logical arguments are used.
#'
#' A negative scale (\eqn{\sigma + \tau \leq 0}) is an error and returns
#' \code{NaN}. For a null scale, the functions consider the limit as
#' \eqn{\sigma + \tau} approaches \code{0},
#' i.e. a point mass at \eqn{m} (\code{location}).
#'
#' @references
#' \insertAllCited{}
#'
#' @source
#' \code{[dpq]anorm} are calculated directly from the definitions using the
#' normal distribution.
#'
#' \code{ranorm} uses the stochastic representation described in the section
#' **Details**. Note that if \eqn{\tau = \sigma}, the function \link[stats]{rnorm}
#' should be used. For \eqn{\tau \neq \sigma}, only two of the three variables
#' in the stochastic representation are needed: \eqn{U} is first generated
#' using \link[stats]{rbinom}, if \eqn{U = 0}, then \eqn{X_{\tau}} is generated
#' using \link[stats]{rnorm} (while \eqn{X_{\sigma}} is taken to be zero),
#' otherwise (when \eqn{U = 1}), \eqn{X_{\sigma}} is generated using
#' \link[stats]{rnorm} (while \eqn{X_{\tau}} is taken to be zero).
#'
#' @seealso See \link[msbreg]{snorm} for a mean-parametrized version
#' of the asymmetric normal distribution.
#'
#' @examples
#' ## Draw the asymmetric normal density with various scales
#' curve(danorm(x, scale1 = 1, scale2 = 1), -8, 8, col = "black", n = 200,
#'       lty = 1, ylim = c(0, 0.4), xlab = "x value", ylab = "Density at x", lwd = 2)
#'
#' curve(danorm(x, scale1 = 2.5, scale2 = 0.5), -8, 8, col = "blue", add = TRUE,
#'       n = 200, lty = 2, lwd = 2)
#'
#' curve(danorm(x, scale1 = 1.5, scale2 = 0.5), -8, 8, col = "red", add = TRUE,
#'       n = 200, lty = 3, lwd = 1)
#'
#' curve(danorm(x, scale1 = 0.5, scale2 = 1.5), -8, 8, col = "blue", add = TRUE,
#'       n = 200, lty = 4, lwd = 1)
#'
#' curve(danorm(x, scale1 = 0.5, scale2 = 2.5), -8, 8, col = "red", add = TRUE,
#'       n = 200, lty = 5, lwd = 2)
#'
#' legend(2, 0.4, legend = paste0("s1 = ", c(1.0, 2.5, 1.5, 0.5, 0.5),
#'                                  ", s2 = ", c(1.0, 0.5, 0.5, 1.5, 2.5)),
#'       lty = 1:5, col = c("black", "blue", "red", "blue", "red"),
#'       lwd = c(2, 2, 1, 1, 2))
#'
#' ## Draw the cumulative distribution functions of
#' # Normal, and ANormal distributions with the same variance
#' SD <- sqrt(0.5 * 1.5 + (1-2/pi))
#' curve(pnorm(x, sd = SD), -5, 5, col = "black",
#'       ylim = c(0, 1), n = 200, xlab = "x value", ylab = "P(X < x)", lwd = 1)
#'
#' curve(panorm(x, scale1 = 0.5, scale2 = 1.5), -5, 5, col = "blue",
#'       add = TRUE, n = 200, lwd = 2, lty = 2)
#'
#' curve(panorm(x, scale1 = 1.5, scale2 = 0.5), -5, 5, col = "red",
#'       add = TRUE, n = 200, lwd = 2, lty = 3)
#'
#' legend(-4, 1, legend = c("Normal(0, SD)", "ANormal(0, 0.5, 1.5)",
#'                          "ANormal(0, 1.5, 0.5)"),
#'       col = c("black", "blue", "red"), lty = 1:3, lwd = c(1,2,2))
#'
#'
#' ## Random variate generation
#' # Compare density curve to histogram
#' rvalues <- ranorm(500000, scale1 = 0.5, scale2 = 1.5)
#' hist(rvalues, freq = FALSE,
#'      ylim = c(0, 0.4), xlim = c(-2, 6), main = "",
#'      breaks = 200, xlab = "x value", ylab = "Density at x")
#'
#' curve(danorm(x, scale1 = 0.5, scale2 = 1.5), -2, 6, n = 200, add = TRUE, col = "blue")
#'
#' # Variance
#' var(rvalues)
#' SD^2 # True variance

# Density function
danorm <- function(x, location = 0, scale1 = 1, scale2 = 1, log = FALSE) {

  # Check numerical arguments
  stopifnot(is.numeric(x), is.numeric(location), is.numeric(scale1), is.numeric(scale2))
  NonValid <- scale1 < 0 | scale2 < 0
  NonValid[is.na(NonValid)] <- FALSE
  scale1[NonValid] <- scale2[NonValid] <- NA

  # Use the max of the dimensions of x, location, scale1 and scale2 (numerical arguments) for the result
  finalres <- pickArgMaxLength(x, location, scale1, scale2)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(x),     length.out = n),
                          rep(c(location),  length.out = n),
                          rep(c(scale1), length.out = n),
                          rep(c(scale2), length.out = n))
  NonValid <- is.na(rowSums(xlocationscale))
  finalres[NonValid] <- NA

  if (!all(NonValid)) {
    x <- (xlocationscale[!NonValid,1] - xlocationscale[!NonValid,2])
    scale1 <- xlocationscale[!NonValid,3]
    scale2 <- xlocationscale[!NonValid,4]
    scale <- scale1 + scale2
    zeroscale <- scale == 0

    if (all(zeroscale)) {
      logpdf <- logb((x == 0) + 0)
    }
    else {
      logpdf <- zeroscale
      if (any(zeroscale))
        logpdf[zeroscale] <- logb((x[zeroscale] == 0) + 0)

      scale_x <- scale1[!zeroscale]
      rightside <- x[!zeroscale] > 0
      scale_x[rightside] <- scale2[!zeroscale][rightside]

      logpdf[!zeroscale] <- stats::dnorm(x=x[!zeroscale], mean = 0,
                                         sd = scale_x, log = TRUE) +
        logb(2) - logb(scale[!zeroscale]) + logb(scale_x)
    }

    finalres[!NonValid] <- if (log[1]) logpdf else exp(logpdf)
  }

  return(finalres)
}

# Distribution function
#' @rdname danorm
panorm <- function(q, location = 0, scale1 = 1, scale2 = 1,
                   lower.tail = TRUE, log.p = FALSE) {

  # Check numerical arguments
  stopifnot(is.numeric(q), is.numeric(location), is.numeric(scale1), is.numeric(scale2))
  NonValid <- scale1 < 0 | scale2 < 0
  NonValid[is.na(NonValid)] <- FALSE
  scale1[NonValid] <- scale2[NonValid] <- NA

  # Use the max of the dimensions of q, location, scale1 and scale2 (numerical arguments) for the result
  finalres <- pickArgMaxLength(q, location, scale1, scale2)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(q),     length.out = n),
                          rep(c(location),  length.out = n),
                          rep(c(scale1), length.out = n),
                          rep(c(scale2), length.out = n))
  NonValid <- is.na(rowSums(xlocationscale))
  finalres[NonValid] <- NA

  if (!all(NonValid)) {
    q <- xlocationscale[!NonValid,1] - xlocationscale[!NonValid,2]
    scale1 <- xlocationscale[!NonValid,3]
    scale2 <- xlocationscale[!NonValid,4]
    scale <- scale1 + scale2
    zeroscale <- scale == 0

    # cdf = P(X > x)
    if (all(zeroscale)) {
      cdf <- (q < 0) + 0
    }
    else {
      cdf <- zeroscale
      if (any(zeroscale))
        cdf[zeroscale] <- (q[zeroscale] >= 0) + 0

      alpha <- scale1[!zeroscale] / scale[!zeroscale]

      rightside <- q[!zeroscale] > 0
      scale_x <- scale1[!zeroscale]
      scale_x[rightside] <- scale2[!zeroscale][rightside]

      cdf[!zeroscale][!rightside] <- stats::pnorm(q[!zeroscale][!rightside]/scale_x[!rightside],
                                                  mean = 0, sd = 1, lower.tail = TRUE, log = FALSE) *
        alpha[!rightside] * 2

      cdf[!zeroscale][rightside] <- (stats::pnorm(q[!zeroscale][rightside]/scale_x[rightside],
                                                  mean = 0, sd = 1, lower.tail = TRUE, log = FALSE) - 0.5) *
        (1 - alpha[rightside]) * 2 + alpha[rightside]
    }

    if (!lower.tail[1])
      cdf <- 1 - cdf

    finalres[!NonValid] <- if (log.p[1]) logb(cdf) else cdf
  }

  return(finalres)
}

# Quantile function
#' @rdname danorm
qanorm <- function(p, location = 0, scale1 = 1, scale2 = 1,
                   lower.tail = TRUE, log.p = FALSE) {

  # Check numerical arguments
  stopifnot(is.numeric(p), is.numeric(location), is.numeric(scale1), is.numeric(scale2))
  NonValid <- scale1 < 0 | scale2 < 0
  NonValid[is.na(NonValid)] <- FALSE
  scale1[NonValid] <- scale2[NonValid] <- NA

  # Use the max of the dimensions of p, location, scale1 and scale2 (numerical arguments) for the result
  finalres <- pickArgMaxLength(p, location, scale1, scale2)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(p),     length.out = n),
                          rep(c(location),  length.out = n),
                          rep(c(scale1), length.out = n),
                          rep(c(scale2), length.out = n))
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
    scale1 <- xlocationscale[!NonValid,3]
    scale2 <- xlocationscale[!NonValid,4]
    scale <- scale1 + scale2
    zeroscale <- scale == 0

    if (all(zeroscale)) {
      z <- 0
    }
    else {
      z <- zeroscale
      if (any(zeroscale))
        z[zeroscale] <- numeric(sum(zeroscale))

      alpha <- scale1[!zeroscale] / scale[!zeroscale]
      rightside <- p[!zeroscale] > alpha

      scale_x <- scale1[!zeroscale]
      scale_x[rightside] <- scale2[!zeroscale][rightside]

      z[!zeroscale][!rightside] <- stats::qnorm(p[!zeroscale][!rightside] / (alpha[!rightside] * 2), mean = 0,
                                                sd = 1, lower.tail = TRUE, log = FALSE)


      z[!zeroscale][rightside] <- stats::qnorm((p[!zeroscale][rightside] - alpha[rightside]) /
                                                 ((1 - alpha[rightside]) * 2) + 0.5,
                                               mean = 0, sd = 1, lower.tail = TRUE, log = FALSE)

      z[!zeroscale] <- z[!zeroscale] * scale_x
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
#' @rdname danorm
ranorm <- function(n, location = 0, scale1 = 1, scale2 = 1) {
  # Check numerical arguments
  stopifnot(is.numeric(location), is.numeric(scale1), is.numeric(scale2))
  NonValid <- scale1 < 0 | scale2 < 0
  NonValid[is.na(NonValid)] <- FALSE
  scale1[NonValid] <- scale2[NonValid] <- NA

  # Use length(n) for the result if length(n) > 1
  if (length(n) > 1)
    n <- length(n)

  # Re-cycle numerical vectors if required
  locationscale <- cbind(rep(c(location),  length.out = n),
                          rep(c(scale1), length.out = n),
                          rep(c(scale2), length.out = n))

  NonValid <- is.na(rowSums(locationscale))
  x <- numeric(n)
  x[NonValid] <- NA

  # Use the stochastic representation with absolute normals
  if (any(!NonValid)) {
    location <- locationscale[!NonValid,1]
    scale1 <- locationscale[!NonValid,2]
    scale2 <- locationscale[!NonValid,3]
    scale <- scale1 + scale2
    zeroscale <- scale == 0

    if (any(zeroscale)) {
      x[!NonValid][zeroscale] <- location[zeroscale]
    }

    if (!all(zeroscale)) {
      alpha <- scale1[!zeroscale] / scale[!zeroscale]
      leftside <- rbinom(sum(!zeroscale), size = 1, prob = alpha) == 1

      scale_x <- scale1[!zeroscale]
      scale_x[!leftside] <- scale2[!zeroscale][!leftside]

      if (any(leftside))
        x[!NonValid][leftside] <- - abs(stats::rnorm(sum(leftside), mean = 0, sd = scale_x[leftside]))

      if (any(!leftside)) {
        x[!NonValid][!leftside] <- abs(stats::rnorm(sum(!leftside), mean = 0, sd = scale_x[!leftside]))
      }

      x[!NonValid][!zeroscale] <- x[!NonValid][!zeroscale] + location[!zeroscale]
    }
  }

  return(x)
}
