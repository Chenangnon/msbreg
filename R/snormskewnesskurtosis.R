
# ' Internal Objects
# '
# ' Numeric matrix for internal use.
# '
# ' @details
# ' The object \code{skt} is a matrix (of class \code{'data.frame'}) of
# ' pre-calculated skewness and excess of kurtosis of the standard Sicard skew
# ' normal distribution (see \link[msbreg]{snorm}).
# '
# ' \itemize{
# ' \item \code{skt$lambda}: values of the skewness parameter \eqn{\lambda} of the Sicard
# ' skew normal distribution with mean 0 and standard
# ' deviation 1. The \eqn{\lambda} values range from \code{-10} to 0, with a
# ' step of \code{0.00001}.
# '
# ' \item \code{skt$skewness}: values of the skewness \eqn{\gamma_1}, i.e. the third
# ' standardized moment of the standard Sicard skew normal distribution
# ' with skewness parameter \eqn{\lambda}: \eqn{\gamma_1 = E \left[ X^3 \right]}
# ' where \eqn{X \sim SN(0,1,\lambda)}.
# '
# ' \item \code{skt$kurtosis}: values of the excess of kurtosis \eqn{\gamma_2}, i.e. the
# ' excess of fourth standardized moment of the standard Sicard skew normal
# ' distribution with skewness parameter \eqn{\lambda} over \code{3} (which is
# ' the fourth standardized moment of the standard normal distribution):
# ' \eqn{\gamma_2 = E \left[ X^4 \right] - 3} where \eqn{X \sim SN(0,1,\lambda)}.
# ' }
# '
# ' This gives the maps from the parameter \eqn{\lambda} to the skewness
# ' \eqn{\gamma_1(\lambda)} and to the excess of kurtosis
# ' \eqn{\gamma_2(\lambda)} values.
# ' Note that for a positive \eqn{\lambda} value, the skewness can be
# ' obtained through \eqn{\gamma_1(\lambda) = - \gamma_1(-\lambda)}, and
# ' the excess of kurtosis through \eqn{\gamma_2(\lambda) = \gamma_2(-\lambda)}.
# '
# ' @docType data
# '
# ' @format An object of class \code{"data.frame"}, with \code{1000001} rows
# ' and \code{3} columns \code{'lambda'}, \code{'skewness'} and \code{'kurtosis'}.
# ' See section **Details**.
# '
# ' @keywords datasets
# '
# ' @source Calculated using \link{integrate} on \link[msbreg]{dsnorm}.
# '
# '
# 'skt'

##lvalues <- seq(-10, 0, 0.00001)
##skvals <- sapply(lvalues,
##                 FUN = function(l) {
##                   integrate(f = function(x) dsnorm(x, lambda = l)*x^3,
##                             lower = -Inf, upper = Inf)$value
##                 })
##
##ktvalues <- sapply(lvalues,
##                   FUN = function(l) {
##                     integrate(f = function(x) dsnorm(x, lambda = l)*x^4,
##                               lower = -Inf, upper = Inf)$value
##                   }) - 3
##
##skewkurt <- as.data.frame(cbind(lambda = lvalues, skewness = skvals, kurtosis = ktvalues))
##
##save(skewkurt, file = "data/skt.rda", version = 2,
##     compress = 'xz', compression_level = 9)
