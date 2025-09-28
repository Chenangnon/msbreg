#' Double Factorial
#'
#' Calculate the (logarithm of the) double factorial of a number.
#'
#' @param n Integers vector.
#'
#' @param log logical, should the logarithm of the double factorial be returned?
#'
#' @details
#' The double factorial of an integer \eqn{n} is defined as
#' \insertCite{gould2012double}{msbreg}
#'
#' \deqn{n!! = n \cdot (n - 2) \cdot (n - 4) \cdot \cdots \cdot 3 \cdot 1} if \eqn{n} is odd, and
#'
#' \deqn{n!! = n \cdot (n - 2) \cdot (n - 4) \cdot \cdots \cdot 4 \cdot 2}
#'
#' if \eqn{n} is even. The double factorial satisfies
#'
#' \deqn{n!! = n \cdot (n - 2)!!}
#'
#' with the convention \eqn{(-1)!! = 0!! = 1}.
#'
#' @return A numeric vector of the (logarithm of the) double factorials of
#' \code{n} values.
#'
#' @source \link[msbreg]{DoubleFactorials} and \link[msbreg]{LogDoubleFactorials}
#' are used for \code{n} values below \code{300} and \code{50000} respectively.
#' The direct computation uses formulae from
#' \insertCite{davis1972gamma;textual}{msbreg}.
#'
#' @references
#' \insertAllCited{}
#'
#' @export DoubleFactorial
#'
DoubleFactorial <- function (n, log = FALSE) {
  n[n < 2] <- 1
  nok <- n <= 300

  if (all(nok)) {
    if (log)
      return(LogDoubleFactorials[n])
    else
      return(DoubleFactorials[n])
  }

  if (!log) {
    warning("n!! is too large to represent for n > 300. Set 'log = TRUE' to obtain the logarithm.")

    out <- rep(Inf, length(n))
    out[nok] <- DoubleFactorials[n[nok]]

    return(out)
  }

  if (all(n <= 50000)) {
    return(LogDoubleFactorials[n])
  }

  even <- n/2
  k <- ceiling(even)
  even <- k == even

  out <- numeric(length(n))

  if (any(even)) {
    out[even] <- lgamma(k[even] + 1) + k[even] * logb(2)
  }

  if (any(!even)) {
    out[!even] <- lgamma(k[!even] + .5) + k[!even] * logb(2) - lgamma(.5)
  }

  return(out)
}


#' Internal Objects
#'
#' Numeric vectors for internal use.
#'
#' @details
#' \code{DoubleFactorials} is a vector of pre-calculated double
#' factorials up to 300!!.
#'
#' \code{LogDoubleFactorials} has the
#' logarithms of double factorial values up to 50 000!!.
#'
#' @docType data
#'
#' @aliases LogDoubleFactorials
#'
#' @usage
#' data(DoubleFactorials)
#'
#' data(LogDoubleFactorials)
#'
#' @format An object of class \code{"numeric"}. It is of length 300
#' for \code{DoubleFactorials} and of length 50 000 for
#' \code{LogDoubleFactorials}.
#'
#' @keywords datasets
#'
#' @source \insertCite{Smith2019TreeTools;textual}{msbreg}
#'
#' @seealso \link[msbreg]{DoubleFactorial}.
#'
#' @references
#' \insertAllCited{}
#'
"DoubleFactorials"

#save(DoubleFactorials, file = "data/DoubleFactorials.rda", version = 2,
#     compress = 'xz', compression_level = 9)
#
#save(LogDoubleFactorials, file = "data/LogDoubleFactorials.rda", version = 2,
#     compress = 'xz', compression_level = 9)

# usethis::use_data(DoubleFactorials,
#                  internal = FALSE,
#                  overwrite = TRUE,
#                  compress = "bzip2",
#                  version = 2,
#                  ascii = FALSE)
