#' Absolute Normal Moments
#'
#' Calculate integer order moments of the absolute value of a centered normal
#' distribution.
#'
#' @param n vector of integers; order(s) of the moment(s).
#' @param sd vector of standard deviations.
#'
#' @details
#' The \eqn{n}th absolute moment of a centered normal variable \eqn{Z} with
#' standard deviation \eqn{\sigma} is given by
#' \insertCite{walck1996hand}{msbreg}.
#'
#'  \deqn{E\left[|Z|^{n}\right] = \left\lbrace \begin{matrix} (n - 1)!! \sigma^n & \mathrm{ for } \ n \ even, \\  \\ \sqrt{\frac{2}{\pi}} 2^{k} k! \sigma^{n} & \mathrm{ for \ odd} \ n = 2k + 1. \end{matrix} \right.}
#'
#' For a normal random variable \eqn{X} with mean \eqn{\mu}
#' and standard deviation \eqn{\sigma}, the general moments are given
#' in terms of the moments of the standard normal variable
#' \eqn{Z = \frac{X - \mu}{\sigma}} as
#'
#' \deqn{E\left[X^{n}\right] = \sum_{k = 0}^n {k \choose n} \mu^n \sigma^{n - k} E\left[Z^{n - k}\right].}
#'
#' The standard normal moments \eqn{E\left[Z^{n - k}\right]} are zero for odd \eqn{n}
#' and can be computed using the function \code{absmnorm} for even \eqn{n}.
#'
#' @return a numeric vector of the absolute moments.
#'
#' @export absmnorm
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' # First six order absolute moments of the standard normal
#' absmnorm (1:6)
#'
absmnorm <- function(n, sd = 1) {
  n <- cbind(n, sd)
  sd <- n[,2]
  n <- n[,1]

  out <- numeric(length(n))
  validn <- n >= 0
  out[!validn] <- NaN
  out[n == 0] <- 1
  if (all(!validn)) {
    return(out)
  }


  even <- n[validn]/2
  k <- floor(even)
  even <- k == even
  if (any(even)) {
    out[validn][even] <- DoubleFactorial(n[validn][even] - 1, log = FALSE)
  }

  if (any(!even)) {
    out[validn][!even] <- exp(lgamma(k[!even] + 1) + (k[!even] + .5) * logb(2) - lgamma(.5))
  }

  out[validn] <- out[validn] * (sd[validn] ^ n[validn])

  return(out)
}
