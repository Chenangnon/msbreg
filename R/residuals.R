#' Extract Model Residuals
#'
#' Methods for extracting model residuals from objects of class \code{'msbm'}.
#'
#' @aliases residuals.msbm
#'
#' @param object an object of class \code{'msbm'}, typically
#' the result of a call to \link[msbreg]{msbreg}.
#'
#' @param type the type of residuals to be returned, the
#' alternative choices being \code{"deviance"} (default),
#' \code{"pearson"}, and \code{"response"}.
#' Can be abbreviated.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' This is a \link[stats]{residuals} method for \code{'msbm'} class objects,
#' mimicking the corresponding methods for "\link[stats]{glm}" class
#' objects (see \link[stats]{residuals.glm}).
#' The motivation and definitions for the three types of residuals
#' (and more) can be found in \insertCite{davison1991residuals;textual}{msbreg}.
#'
#' @return Residuals extracted from the fit \code{object}.
#'
#' @seealso
#' See \link[msbreg]{summary.msbm} for summarizing
#' a multistage binomial model fit.
#'
#' @exportS3Method stats::residuals
#'
#' @references
#' \insertAllCited{}
#'
residuals.msbm <- function (object,
                            type = c("deviance",
                                     "pearson",
                                     "response"),
                            ...) {
  type <- match.arg(type)
  y.resid <- object$y.resid
  mu <- object$fitted.values
  wts <- object$weights

  res <- switch(type,
                deviance = if (object$df.residual > 0) object$deviance.resid * sign(y.resid)
                else rep.int(0, length(mu)),
                pearson = y.resid * sqrt(wts)/sqrt(mu * (1 - mu)),
                response = y.resid)

  if (!is.null(object$na.action))
    res <- naresid(object$na.action, res)

  res
}
#setMethod ("residuals",
#           signature = "msbm",
#           definition = residuals.msbm)
