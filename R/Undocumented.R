#' Undocumented objects
#'
#' Page for undocumented objects in \code{msbreg} (included to avoid warning during R CMD check).
#'
#' @details
#' The current undocumented objects are (the list changes as we add/remove
#' features to/from the package):
#' \describe{
#' \item{\code{skt}}{ numeric matrix for internal use.
#' The object \code{skt} is a matrix (of class \code{'data.frame'}) of
#' pre-calculated skewness and excess of kurtosis of the standard Sicard skew
#' normal distribution (see \link[msbreg]{snorm}).}
#' }
#'
#' @export undocumented
#' @aliases skt
undocumented <- function() {
  c('skt')
}
