#' Coefficient of Skewness
#'
#' Estimation of the coefficient of skewness.
#'
#' @param x a numeric \code{vector}, \code{matrix} or \code{data.frame}.
#'
#' @param na.rm logical, should missing values be removed?
#'
#' @param type character, type of skewness to compute. One of \code{'G1'}
#' (the adjusted Fisher-Pearson standardized moment coefficient which uses both
#' unbiased variance and third central moment estimates),
#' \code{'b1'} (which uses unbiased variance estimate, but biased third central
#' moment estimate), \code{'g1'} (the traditional Fisher-Pearson
#' coefficient of skewness, which uses both biased variance and third central
#' moment estimates), and \code{'all'} (all the three previous types).
#'
#' @details
#' The skewness is a measure of the asymmetry of a distribution. For a random
#' variable \eqn{X}, the skewness \eqn{\gamma_1} (the third standardized moment)
#' is defined as \insertCite{pearson1895contributions,doane2011measuring}{msbreg}
#'
#' \deqn{\gamma_1 = \frac{E\left[ (X - \mu)^3 \right]}{E\left[ (X - \mu)^2 \right]^{3/2}}}
#'
#' where \eqn{\mu = E\left[ X \right]} and all expectations are assumed finite.
#' The population skewness \eqn{\gamma_1} is zero for a *symmetric* distribution,
#' negative for a *left-skewed* distribution, and positive for a *right-skewed*
#' distribution.
#'
#' The three provided sample coefficients of skewness are the common
#' estimators, but are generally biased estimates of \eqn{\gamma_1}, unless
#' \eqn{X} follows a normal distribution.
#'
#' The function \code{skewness} is a simple alias for \code{skew}.
#'
#' @return a vector or matrix of sample coefficients of skewness.
#'
#' If argument \code{x} is a \code{vector}, the output is a vector of three skewness
#' coefficient estimates if \code{type = 'all'}, otherwise, the output is
#' a scalar, the specified \code{type} estimate of the skewness coefficient.
#'
#' If argument \code{x} is a \code{matrix} or a \code{data.frame}, the output
#' is a matrix where each column is the result from skewness applied to each
#' column of \code{x} separately.
#'
#' In all cases, the returned object has attributes \code{'se'} (estimated
#' asymptotic standard error(s) of the skewness estimates, vector or matrix of
#' the same size as the output) and \code{'n'} (vector of the size(s) of
#' samples used for each column of \code{x}).
#'
#' @references
#' \insertAllCited{}
#'
#' @export skewness
#' @export skew
#' @aliases skewness
#'
#' @examples
#' ## Gaussian samples
#' rGaussian <- rnorm(500000, sd = 1)
#'
#' skew (rGaussian, type = 'G1')
#'
#' ## Asymmetric Gaussian samples
#' rAGaussian <- ranorm(500000, scale1 = 0.5, scale2 = 1.5)
#'
#' skew (rAGaussian, type = 'all')
#'
#' skew (cbind(rGaussian, rAGaussian), type = 'all')
#'
skew <- function(x, na.rm = FALSE, type = c('G1', 'b1', 'g1', 'all')) {

  type <- match.arg(type)

  if (is.matrix(x) | is.data.frame(x)) {
    skews <- apply(x, MARGIN = 2, FUN = function(x) {
      ret <- skewness (x, na.rm = na.rm, type = type)
      c(ret, attr(ret, 'se'), attr(ret, 'n'))
    })

    p <- NCOL(x)

    if (p > 1) {
      if (identical(type, 'all')) {
        out <- skews[1:3,]

        attr(out, 'se') <- skews[4:6,]

        attr(out, 'n') <- skews[7,]
      }
      else {
        out <- skews[1,]

        attr(out, 'se') <- skews[2,]

        attr(out, 'n') <- skews[3,]
      }
    }
    else {
      if (identical(type, 'all')) {
        out <- skews[1:3]

        attr(out, 'se') <- skews[4:6]

        attr(out, 'n') <- skews[7]
      }
      else {
        out <- skews[1]

        attr(out, 'se') <- skews[2]

        attr(out, 'n') <- skews[3]
      }
    }

    return(out)
  }

  stopifnot(is.numeric(x))

  m1 <- mean(x, na.rm = FALSE)
  m3 <- x - m1
  m2 <- mean(m3^2, na.rm = FALSE)
  m3 <- mean(m3^3, na.rm = FALSE)

  g1 <- as.numeric(m3 / m2^1.5)

  if (na.rm)
    n <- length(na.omit(x))
  else
    n <- length(x)

  if (n > 2)
    se <- sqrt(6 * (n - 2) / ((n + 1) * (n + 3)))
  else
    se <- NA

  switch(type,
         g1 = {
           out <- structure(g1, se = se, n = n)
         },
         G1 = {
           #G1 <- g1 * sqrt(n * (n - 1)) / (n - 2)
           kg <- sqrt(n * (n - 1)) / (n - 2)
           out <- structure(g1 * kg, se = se, n = n)

           if (n > 2)
             se <- se * kg
         },
         b1 = {
           #b1 <- g1 * ((1 - 1 / n))^1.5
           kb <- ((1 - 1 / n))^1.5
           out <- structure(g1 * kb, se = se, n = n)

           if (n > 2)
             se <- se * kb
         },
         all = {
           kg <- sqrt(n * (n - 1)) / (n - 2)
           kb <- ((1 - 1 / n))^1.5

           out <- structure(c(G1 = g1 * kg,
                              b1 = g1 * kb,
                              g1 = g1),
                            se = se * c(G1 = kg, b1 = kb, g1 = 1),
                            n = n)
         })

  return(out)
}

#' @rdname skew
skewness <- skew
