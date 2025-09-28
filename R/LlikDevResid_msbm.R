#' Assessing Multistage Binomial Model Fits
#'
#' Methods for assessing class \code{msbm} objects.
#'
#' @aliases logLik.msbm deviance.msbm
#' @aliases coef.msbm link.msbm
#' @aliases rsquared.msbm
#
#' @param object an object of class \code{msbm}, typically
#' the result of a call to \link[msbreg]{msbreg}.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @param method character indicating the desired pseudo R squared statistic(s).
#' Any of the following alternatives:
#'
#' \itemize{
#' \item \code{'KL'} for Kullback–Leibler divergence based
#' \eqn{R^2} (default), equivalent to the deviance reduction ratio achieved
#' by the predictors in a fitted model as compared to a model with no predictor;
#'
#' \item \code{'Nagelkerke'} for Nagelkerke's \eqn{R^2};
#'
#' \item \code{'COR'} for the squared linear/monotonic correlation
#' coefficient (computed by \link[stats]{cor}) between the observed
#' responses and the corresponding fitted values.
#' }
#'
#' @param cor.method character indicating which correlation coefficient
#' to consider (only used when \code{method = 'COR'}), this is passed
#' to \link[stats]{cor} as the \code{method} argument
#' (\code{method = cor.method}).
#'
#' @param adjust.size logical, should the pseudo R-squared statistic be adjusted
#' for model complexity (number of free parameters)? If \code{TRUE}, the returned
#' pseudo R-squared statistic is \eqn{1 - \frac{n-1}{n-p} (1 - R^2)} where
#' \eqn{R^2} is the unadjusted statistic, \eqn{n} is the sample size (number of
#' observations) and \eqn{p} is the number of free model parameters.
#'
#' @details
#' The functions \code{coef}, \code{deviance} and \code{logLik} mimick the
#' corresponding methods for "\link[stats]{glm}" class objects.
#' The \code{link} method mimicks the \link[stats]{family} function.
#' The function \code{rsquared} is a method for \link[msbreg]{rsquared}.
#'
#' @seealso
#' See \link[msbreg]{summary.msbm} for summarizing
#' a multistage binomial model fit.
#'
#' @exportS3Method stats::logLik
logLik.msbm <- function(object, ...) {
  out <- attr(object$fisher, "logLik") # %||% NA
  if (is.null(out)) {
    out <- object$dims$npars - object$aic/2
  }
  attr(out, "nobs") <- sum(!is.na(object$y.resid))
  attr(out, "df") <- object$dims$npars
  class(out) <- "logLik"

  return(out)
}
#setMethod ("logLik",
#           signature = "msbm",
#           definition = logLik.msbm)

#' @exportS3Method stats::logLik
logLik.glm2msbm <- function(object, ...) {
  logLik.msbm (object)
}
#setMethod ("logLik",
#           signature = "glm2msbm",
#           definition = logLik.glm2msbm)

#' @rdname logLik.msbm
#' @exportS3Method stats::deviance
deviance.msbm <- function (object, ...) {
  object$deviance
}
#setMethod ("deviance",
#           signature = "msbm",
#           definition = deviance.msbm)

#' @rdname logLik.msbm
#' @exportS3Method stats::coef
coef.msbm <- function (object, ...) {
  object$coef
}
#setMethod ("coef",
#           signature = "msbm",
#           definition = coef.msbm)

#' @rdname logLik.msbm
#' @exportS3Method msbreg::link
link.msbm <- function (object, ...)
  object$link

### setMethod commented out to prevent the warning:
### Undocumented S4 methods:
###      generic 'link' and siglist 'msbm'
#setMethod ("link",
#           signature = "msbm",
#           definition = link.msbm)

#' @rdname logLik.msbm
#' @exportS3Method msbreg::rsquared
rsquared.msbm <- function(object, ...,
                          method = 'KL',
                          cor.method = "pearson",
                          adjust.size = FALSE) {
  # Save the call (for naming purposes)
  Call <- match.call()

  out <- rsquared_default_core(object = object,
                               ...,
                               method = method,
                               cor.method = cor.method,
                               adjust.size = adjust.size)

  rownames(out) <- as.character(Call[-1L])[1:NROW(out)]

  return(out)
}
