#'
#' The Sicard Skew Normal Distribution
#'
#' Density, distribution function, quantile function and random generation for
#' the Sicard skew normal distribution with mean equal to \code{mean},
#' standard deviation \code{sd} and skewness parameter \code{lambda}.
#'
#' @param x,q vector of quantiles.
#'
#' @param p vector of probabilities.
#'
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#' to be the number required (\code{n = length(n)}).
#'
#' @param mean vector of means (real).
#'
#' @param sd vectors of standard deviations (positive).
#'
#' @param lambda vector of real skewness parameter (real).
#'
#' @param log,log.p logical; if \code{TRUE}, densities and probabilities \code{p}
#' are given as \code{log(p)}.
#'
#' @param lower.tail logical; if \code{TRUE} (default), probabilities are
#' \eqn{P[X \leq x]} otherwise, \eqn{P[X > x]}.
#'
#' @export dsnorm
#' @export psnorm
#' @export qsnorm
#' @export rsnorm
#'
#' @aliases psnorm
#' @aliases qsnorm
#' @aliases rsnorm
#' @aliases snorm
#'
#' @details
#' The asymmetric normal distribution \eqn{AN(m, \sigma, \tau)} was proposed
#' by \insertCite{sicard2013loinormale;textual}{msbreg}. The distribution is
#' also known as the epsilon–skew–normal distribution
#' \insertCite{mudholkar2000epsilon}{msbreg}. The skew normal
#' distribution \eqn{SN(\mu, \sigma, \lambda)} presented here is a
#' mean-parametrization of the asymmetric normal distribution.
#' For \eqn{x \in \mathbb{R}}, a \eqn{SN(\mu, \sigma, \lambda)} variable
#' \eqn{X} has density function
#'
#' \deqn{f(x) = \left\lbrace \begin{matrix} \quad \quad 2 \alpha \phi \left(x|m, \alpha \omega\right) & \mathrm{ if } \quad x \leq m \\  \\ 2 (1 - \alpha)  \phi \left(x|m, (1 - \alpha) \omega\right)  & \mathrm{ if } \quad x > m \end{matrix} \right.}
#'
#' where \eqn{\alpha = 1/\left( 1 + e^{2\lambda} \right)} is the
#' probability mass under the left side of the mode \eqn{m} of the density
#' function, \eqn{m = \mu - \sqrt{\frac{2}{\pi}} (1 - 2 \alpha) \omega},
#' \eqn{\omega = \sigma /\sqrt{\alpha (1 - \alpha) + (1 - \frac{2}{\pi}) (1 - 2 \alpha)^2}},
#' \eqn{\phi(\cdot|m,s)} denotes the density of the normal distribution
#' with mean \eqn{m} and standard deviation \eqn{s}.
#' The cumulative distribution function (cdf) \eqn{F(x)=P\left[X \leq x \right]}
#' is given by
#'
#' \deqn{F(x) = \left\lbrace \begin{matrix} \quad \quad 2 \alpha \Phi \left( \frac{x - m}{\alpha \omega} \right) & \mathrm{ if } \quad x \leq m \\  \\ \alpha + (1 - \alpha) \left[ 2 \Phi \left( \frac{x - m}{(1 - \alpha) \omega} \right) - 1 \right]  & \mathrm{ if } \quad x > m \end{matrix} \right.}
#'
#' where \eqn{\Phi} denotes the cdf of the standard normal distribution.
#' The quantile function of \eqn{X} is for \eqn{p \in (0, 1)}:
#'
#' \deqn{F^{-1}(p) = m + \sigma_x \Phi^{-1}(v)}
#'
#' where \eqn{\left\lbrace \begin{matrix} \sigma_x = \alpha \omega & \mathrm{ and } \quad v = \frac{p}{2 \alpha} & \mathrm{ if } \quad p \leq \alpha, \\  \\ \sigma_x = (1 - \alpha) \omega & \mathrm{ and } \quad v = \frac{p - \alpha}{2(1 - \alpha)} + \frac{1}{2}  & \mathrm{ if } \quad p > \alpha. \end{matrix} \right.}
#'
#' The \eqn{SN(\mu, \sigma, \lambda)} distribution has expectation
#' \eqn{\mu} and variance \eqn{\sigma^2}.
#' More properties of a variable \eqn{X \sim SN(\mu, \sigma, \lambda)} can
#' be derived from the following stochastic representation:
#'
#' \deqn{X \overset{d}{=} -U |X_{\alpha}| + (1 - U) |X_{1-\alpha}| + m}
#'
#' where \eqn{\overset{d}{=}} means "is equal in distribution to",
#' \eqn{U} is a Bernoulli variable with success probability
#' \eqn{\alpha}, \eqn{|X_{\alpha}|} is, for
#' \eqn{\alpha \in (0, 1)}, the absolute value
#' of a Gaussian variable with mean zero and standard deviation \eqn{\alpha \omega},
#' and \eqn{U}, \eqn{X_{\alpha}} and \eqn{X_{1-\alpha}} are mutually
#' independent \insertCite{sicard2013loinormale}{msbreg}. For instance, it
#' follows that the \eqn{k}th moment of the distribution is given when
#' \eqn{m = 0} by
#'
#' \deqn{E \left[ X^k \right] = - \alpha E \left[ |X_{\alpha}|^k \right] + (1 - \alpha) E \left[ |X_{1-\alpha}|^k \right]}
#'
#' where \eqn{E \left[ |X_{\alpha}|^k \right]} is the \eqn{k}th absolute moment
#' of a centered normal variable (this can be easily obtained from
#' \link[msbreg]{absmnorm}). Moments for general \eqn{m} values can
#' be obtained using a binomial expansion.
#'
#' When \code{mean} or \code{sd} or \code{lambda} are not specified
#' they respectively assume the default values of \eqn{0}, \eqn{1} and \eqn{0},
#' which correspond to the standard normal distribution.
#'
#' @return \code{dsnorm} gives the density, \code{psnorm} gives the distribution
#' function, \code{qsnorm} gives the quantile function, and \code{rsnorm}
#' generates random variates.
#'
#' The length of the result is determined by \code{n} for \code{rsnorm}, and
#' is the maximum of the lengths of the numerical arguments for the other
#' functions.
#'
#' The numerical arguments other than \code{n} are recycled to the length of the result.
#' Only the first elements of the logical arguments are used.
#'
#' A negative standard deviation (\eqn{\sigma \leq 0}) is an error and returns
#' \code{NaN}. For a null standard deviation, the functions consider the limit as
#' \eqn{\sigma} approaches \code{0},
#' i.e. a point mass at \eqn{\mu}.
#'
#' @references
#' \insertAllCited{}
#'
#' @source
#' \code{[dpq]snorm} are calculated directly from the definitions using the
#' normal distribution.
#'
#' \code{rsnorm} uses the stochastic representation described in the section
#' **Details**. Note that if \eqn{\lambda = 0}, the function \link[stats]{rnorm}
#' should be used. For \eqn{\lambda \neq 0}, only two of the three variables
#' in the stochastic representation are needed: \eqn{U} is first generated
#' using \link[stats]{rbinom}, if \eqn{U = 0}, then \eqn{X_{1-\alpha}} is generated
#' using \link[stats]{rnorm} (while \eqn{X_{\alpha}} is taken to be zero),
#' otherwise (when \eqn{U = 1}), \eqn{X_{\alpha}} is generated using
#' \link[stats]{rnorm} (while \eqn{X_{1-\alpha}} is taken to be zero).
#'
#' @seealso See \link[msbreg]{anorm} for the original parametrization
#' of the asymmetric normal distribution.
#'
#' @examples
#' ## Draw the asymmetric normal density with various scales
#' curve(dsnorm(x, sd = 1, lambda = 0), -5, 5, col = "black", n = 200,
#'       lty = 1, ylim = c(0, 0.44), xlab = "x value", ylab = "Density at x", lwd = 2)
#'
#' curve(dsnorm(x, sd = 1, lambda = -0.5), -5, 5,, col = "blue", add = TRUE,
#'       n = 200, lty = 2, lwd = 2)
#'
#' curve(dsnorm(x, sd = 1, lambda = -1), -5, 5,, col = "red", add = TRUE,
#'       n = 200, lty = 3, lwd = 1)
#'
#' curve(dsnorm(x, sd = 1, lambda = 0.5), -5, 5, col = "blue", add = TRUE,
#'       n = 200, lty = 4, lwd = 1)
#'
#' curve(dsnorm(x, sd = 1, lambda = 1), -5, 5,, col = "red", add = TRUE,
#'       n = 200, lty = 5, lwd = 2)
#'
#' legend(-4.8, 0.4, legend = paste0("lambda = ", c(0, -0.5, -1, 0.5, 1)),
#'       lty = 1:5, col = c("black", "blue", "red", "blue", "red"),
#'       lwd = c(2, 2, 1, 1, 2))
#'
#' ## Draw the cumulative distribution functions of
#' # Normal, and ANormal distributions with the same variance
#' curve(pnorm(x, sd = 1), -4, 4, col = "black",
#'       ylim = c(0, 1), n = 200, xlab = "x value", ylab = "P(X < x)", lwd = 1)
#'
#' curve(psnorm(x, sd = 1, lambda = 2), -4, 4, col = "blue",
#'       add = TRUE, n = 200, lwd = 2, lty = 2)
#'
#' curve(psnorm(x, sd = 1.5, lambda = -2), -4, 4, col = "red",
#'       add = TRUE, n = 200, lwd = 2, lty = 3)
#'
#' legend(-3.9, 1, legend = c("Normal(0, 1)", "ANormal(0, 1, 2)",
#'                          "ANormal(0, 1, -2)"),
#'       col = c("black", "blue", "red"), lty = 1:3, lwd = c(1,2,2))
#'
#'
#' ## Random variate generation
#' # Compare density curve to histogram
#' rvalues <- rsnorm(500000, sd = 1, lambda = 1)
#' hist(rvalues, freq = FALSE,
#'      ylim = c(0, 0.5), xlim = c(-2, 5), main = "",
#'      breaks = 200, xlab = "x value", ylab = "Density at x")
#'
#' curve(dsnorm(x, sd = 1, lambda = 1), -2, 5, n = 200,
#'       add = TRUE, col = "red", lwd = 2)
#'
#' # Variance
#' var(rvalues)
#' 1 # True variance

# Density function
dsnorm <- function(x, mean = 0, sd = 1, lambda = 0, log = FALSE) {

  # Check numerical arguments
  stopifnot(is.numeric(x), is.numeric(mean), is.numeric(sd), is.numeric(lambda))
  NonValid <- sd < 0
  NonValid[is.na(NonValid)] <- FALSE
  sd[NonValid] <- NA

  # Use the max of the dimensions of x, mean, sd (numerical arguments) for the result
  finalres <- pickArgMaxLength(x, mean, sd, lambda)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(x),     length.out = n),
                          rep(c(mean),  length.out = n),
                          rep(c(sd), length.out = n),
                          rep(c(lambda), length.out = n))
  NonValid <- is.na(rowSums(xlocationscale))
  finalres[NonValid] <- NA

  if (!all(NonValid)) {
    x <- (xlocationscale[!NonValid,1] - xlocationscale[!NonValid,2])
    sd <- xlocationscale[!NonValid,3]
    lambda <- xlocationscale[!NonValid,4]
    zeroscale <- sd == 0

    if (all(zeroscale)) {
      logpdf <- logb((x == 0) + 0)
    }
    else {
      logpdf <- zeroscale
      if (any(zeroscale))
        logpdf[zeroscale] <- logb((x[zeroscale] == 0) + 0)

      if (all(lambda == 0)) {
        logpdf[!zeroscale] <- stats::dnorm (x = x[!zeroscale], mean = 0,
                                            sd = sd[!zeroscale], log = TRUE)
      }
      else {
        # Shape alpha
        #alpha <- (tanh(-lambda[!zeroscale]) + 1)/2 # Would be faster? No idea
        alpha <- 1 / (1 + exp(2*lambda[!zeroscale]))

        # Scale omega
        omega <- sd[!zeroscale] / sqrt(alpha * (1 - alpha) + (1 - 2/pi) * (1 - 2 * alpha)^2)

        # Location m: subtract m from x: c <- sqrt(2/pi); m <- - c * omega * (1 - 2 * alpha)
        x[!zeroscale] <- x[!zeroscale] + sqrt(2/pi) * omega * (1 - 2 * alpha)

        # Standard deviation weight
        rightside <- x[!zeroscale] > 0
        weight <- alpha
        if (any(rightside))
          weight[rightside] <- (1 - alpha[rightside])

        logpdf[!zeroscale] <- stats::dnorm(x = x[!zeroscale], mean = 0,
                                           sd = weight * omega, log = TRUE) +
          logb(2) + logb(weight)
      }
    }

    finalres[!NonValid] <- if (log[1]) logpdf else exp(logpdf)
  }

  return(finalres)
}

# Distribution function
#' @rdname dsnorm
psnorm <- function(q, mean = 0, sd = 1, lambda = 0,
                    lower.tail = TRUE, log.p = FALSE) {

  # Check numerical arguments
  stopifnot(is.numeric(q), is.numeric(mean), is.numeric(sd), is.numeric(lambda))
  NonValid <- sd < 0
  NonValid[is.na(NonValid)] <- FALSE
  sd[NonValid] <- NA

  # Use the max of the dimensions of q, mean, sd (numerical arguments) for the result
  finalres <- pickArgMaxLength(q, mean, sd, lambda)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(q),     length.out = n),
                          rep(c(mean),  length.out = n),
                          rep(c(sd), length.out = n),
                          rep(c(lambda), length.out = n))
  NonValid <- is.na(rowSums(xlocationscale))
  finalres[NonValid] <- NA

  if (!all(NonValid)) {
    q <- xlocationscale[!NonValid,1] - xlocationscale[!NonValid,2]
    sd <- xlocationscale[!NonValid,3]
    lambda <- xlocationscale[!NonValid,4]
    zeroscale <- sd == 0

    # cdf = P(X > x)
    if (all(zeroscale)) {
      cdf <- (q < 0) + 0
    }
    else {
      cdf <- zeroscale

      if (any(zeroscale))
        cdf[zeroscale] <- logb((q[zeroscale] == 0) + 0)

      if (all(lambda[!zeroscale] == 0)) {
        cdf[!zeroscale] <- stats::pnorm (q = q[!zeroscale], mean = 0,
                                         sd = sd[!zeroscale], log = FALSE)
      }
      else {
        # Shape alpha
        #alpha <- (tanh(-lambda[!zeroscale]) + 1)/2 # Would be faster? No idea
        alpha <- 1 / (1 + exp(2*lambda[!zeroscale]))

        # Scale omega
        omega <- sd[!zeroscale] / sqrt(alpha * (1 - alpha) + (1 - 2/pi) * (1 - 2 * alpha)^2)

        # Location m: subtract m from x: c <- sqrt(2/pi); m <- - c * omega * (1 - 2 * alpha)
        q[!zeroscale] <- q[!zeroscale] + sqrt(2/pi) * omega * (1 - 2 * alpha)

        # Standard deviation weight
        rightside <- q[!zeroscale] > 0
        if (any(!rightside)) {
          cdf[!zeroscale][!rightside] <- stats::pnorm(q = q[!zeroscale][!rightside] / (alpha[!rightside] * omega[!rightside]),
                                                      mean = 0, sd = 1, log = FALSE) * alpha[!rightside] * 2

        }

        if (any(rightside)) {
          cdf[!zeroscale][rightside] <- (stats::pnorm(q = q[!zeroscale][rightside] / ((1 - alpha[rightside]) * omega[rightside]),
                                                      mean = 0, sd = 1, log = FALSE) - 0.5) * (1 - alpha[rightside]) * 2 +
            alpha[rightside]
        }
      }
    }

    if (!lower.tail[1])
      cdf <- 1 - cdf

    finalres[!NonValid] <- if (log.p[1]) logb(cdf) else cdf
  }

  return(finalres)
}

# Quantile function
#' @rdname dsnorm
qsnorm <- function(p, mean = 0, sd = 1, lambda = 0,
                    lower.tail = TRUE, log.p = FALSE) {

  # Check numerical arguments
  stopifnot(is.numeric(p), is.numeric(mean), is.numeric(sd), is.numeric(lambda))
  NonValid <- sd < 0
  NonValid[is.na(NonValid)] <- FALSE
  sd[NonValid] <- NA

  # Use the max of the dimensions of p, mean, sd (numerical arguments) for the result
  finalres <- pickArgMaxLength(p, mean, sd, lambda)

  # Re-cycle numerical vectors if required
  n <- length(finalres)
  xlocationscale <- cbind(rep(c(p),     length.out = n),
                          rep(c(mean),  length.out = n),
                          rep(c(sd), length.out = n),
                          rep(c(lambda), length.out = n))
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

    sd <- xlocationscale[!NonValid,3]
    lambda <- xlocationscale[!NonValid,4]
    zeroscale <- sd == 0

    if (all(zeroscale)) {
      z <- 0
    }
    else {
      z <- zeroscale
      if (any(zeroscale))
        z[zeroscale] <- numeric(sum(zeroscale))

      if (all(lambda[!zeroscale] == 0)) {
        z[!zeroscale] <- stats::qnorm (p = p[!zeroscale], mean = 0,
                                         sd = sd[!zeroscale], log = FALSE)
      }
      else {
        # Shape alpha
        alpha <- 1 / (1 + exp(2*lambda[!zeroscale]))

        # Scale omega
        omega <- sd[!zeroscale] / sqrt(alpha * (1 - alpha) + (1 - 2/pi) * (1 - 2 * alpha)^2)

        # Quantile
        rightside <- p[!zeroscale] > alpha
        if (any(!rightside)) {
          z[!zeroscale][!rightside] <- stats::qnorm (p = 0.5 * p[!zeroscale][!rightside] / alpha[!rightside],
                                                     mean = 0, sd = 1, log = FALSE)

          z[!zeroscale][!rightside] <- z[!zeroscale][!rightside] * alpha[!rightside] * omega[!rightside]
        }

        if (any(rightside)) {
          z[!zeroscale][rightside] <- stats::qnorm (p = 0.5 * ((p[!zeroscale][rightside] - alpha[rightside]) / (1 - alpha[rightside]) + 1),
                                                     mean = 0, sd = 1, log = FALSE)

          z[!zeroscale][rightside] <- z[!zeroscale][rightside] * (1 - alpha[rightside]) * omega[rightside]
        }

        # Location m: subtract m from z: c <- sqrt(2/pi); m <- - c * omega * (1 - 2 * alpha)
        z[!zeroscale] <- z[!zeroscale] - sqrt(2/pi) * omega * (1 - 2 * alpha)
      }
    }

    finalres[!NonValid] <- z + xlocationscale[!NonValid,2]
  }

  if (any(zerop & !NonValid0))
    finalres[zerop & !NonValid0] <- -Inf

  if (any(onep & !NonValid0))
    finalres[onep & !NonValid0] <- Inf

  return(finalres)
}

# Random variate generator
#' @rdname dsnorm
rsnorm <- function(n, mean = 0, sd = 1, lambda = 0) {
  # Check numerical arguments
  stopifnot(is.numeric(mean), is.numeric(sd), is.numeric(lambda))
  NonValid <- sd < 0
  NonValid[is.na(NonValid)] <- FALSE
  sd[NonValid] <- NA

  # Use length(n) for the result if length(n) > 1
  if (length(n) > 1)
    n <- length(n)

  # Re-cycle numerical vectors if required
  locationscale <- cbind(rep(c(mean),  length.out = n),
                         rep(c(sd), length.out = n),
                         rep(c(lambda), length.out = n))
  NonValid <- is.na(rowSums(locationscale))
  x <- numeric(n)
  x[NonValid] <- NA

  # Use the stochastic representation with absolute normals
  if (any(!NonValid)) {

    sd <- locationscale[!NonValid,2]
    lambda <- locationscale[!NonValid,3]
    zeroscale <- sd == 0
    xz <- numeric(sum(!NonValid))

    if (!all(zeroscale)) {
      if (all(lambda[!zeroscale] == 0)) {
        xz[!zeroscale] <- stats::rnorm (n = sum(!zeroscale), mean = 0,
                                        sd = sd[!zeroscale])
      }
      else {
        # Shape alpha
        alpha <- 1 / (1 + exp(2*lambda[!zeroscale]))

        # Scale omega
        omega <- sd[!zeroscale] / sqrt(alpha * (1 - alpha) + (1 - 2/pi) * (1 - 2 * alpha)^2)

        # Bernoulli variable
        leftside <- stats::rbinom(n = sum(!zeroscale), size = 1, prob = alpha) == 1

        if (any(leftside)) {
          xz[!zeroscale][leftside] <- - abs(stats::rnorm (n = sum(leftside), mean = 0,
                                                          sd = omega[leftside] * alpha[leftside]))
        }

        if (any(!leftside)) {
          xz[!zeroscale][!leftside] <- abs(stats::rnorm (n = sum(!leftside), mean = 0,
                                                         sd = omega[!leftside] * (1 - alpha[!leftside])))
        }

        xz[!zeroscale] <- xz[!zeroscale] - sqrt(2/pi) * omega * (1 - 2 * alpha)
      }
    }

    x[!NonValid] <- xz + locationscale[!NonValid,1]
  }

  return(x)
}
