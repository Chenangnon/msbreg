#' Empirical Estimating Functions and Bread Matrix
#'
#' Functions to build matrix components of sandwich covariance
#' estimators: \code{estfun} extracts empirical estimating functions
#' (EEFs) and \code{bread} extracts an estimator of the scaled expected
#' information matrix from fitted Multistage Binomial (MSB)
#' model objects.
#'
#' @param x an object of class \code{"msbm"}, a fitted MSB model.
#'
#' @param coefficients a numeric vector of all MSB model parameters.
#' When supplied, \code{coefficients} must have finite values,
#' and the same length as \code{object$coefficients}.
#'
#' @param fisher.matrix logical, should the Fisher (expected) information
#' matrix be used to obtain the bread matrix? Defaults to \code{TRUE}.
#' The alternative is the observed information matrix.
#'
#' @param adjust logical, should the EEFs (for \code{estfun}) or the
#' bread matrix (for \code{bread}) be adjusted for any penalty
#' added to the log-likelihood function during model fitting?
#' Defaults to \code{TRUE}.
#'
#' @param ... further arguments passed to or from other methods.
#' None is currently used.
#'
#' @details
#' These are methods for the package \code{sandwich}'s generics
#' \link[sandwich]{estfun} and \link[sandwich]{bread}, to handle
#' objects of class \code{"msbm"}. The methods allow to use the
#' \link[sandwich]{bread}, \link[sandwich]{meat}, and
#' \link[sandwich]{sandwich} functions and thus construct any
#' modified sandwich covariance estimator for objects of class
#' \code{"msbm"} (see for instance
#' \insertCite{carroll1998sandwich,skene2010analysis;textual}{msbreg}).
#'
#' The argument \code{adjust} is by default set to \code{adjust = TRUE}
#' to ensure that the exact estimating functions are returned, in line
#' with the default \link[sandwich]{estfun} and \link[sandwich]{bread}
#' methods. However, \code{adjust = FALSE} is required to obtain the
#' estimating functions used by default by \link[msbreg]{varcov}.
#' In general, if likelihood inference is the focus, \code{adjust = FALSE}
#' should be used since a penalty in the likelihood surface generally
#' leads to an estimate with a reduced curvature as compared to the
#' maximum likelihood point (if any it exists). The resulting lost of
#' efficiency can be masked when using \code{adjust = TRUE}.
#'
#' @return For \code{estfun}, the empirical estimating functions, an
#' \eqn{n \times p} matrix where \eqn{n} is the number of observations
#' and \eqn{p} is the length of \code{coefficients} (the columns are named
#' after the \code{$coefficients} component of \code{x}). These are the
#' scores of the individual contributions to the (penalized) log-likelihood
#' (derivatives with respect to model parameters).
#'
#' For \code{bread}, a \eqn{p \times p} matrix where \eqn{p} is the length
#' of \code{coefficients}: an estimate of the expectation of the the negative
#' hessian of the (penalized) log-likelihood function (Jacobian of the
#' estimating functions with respect to model parameters).
#' If \code{fisher.matrix = FALSE}, the result is the observed hessian while
#' for \code{fisher.matrix = TRUE}, the expected hessian is returned.
#' The rows and columns are named after the \code{$coefficients} component of
#' \code{x}.
#'
#' @aliases bread.meat
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \link[msbreg]{varcov} for variance-covariance matrix estimation.
#'
#' @examples
#' data("test1data")
#'
#' MSBfit = msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                  weights = Total,
#'                  data = test1data)
#'
#' round(vcov (MSBfit, sandwich = TRUE, adjust = TRUE), 5)
#'
# do.call(require, list("sandwich", character.only=TRUE))
#' \dontrun{
#' require (sandwich)
#' round(sandwich (MSBfit), 5)
#' }
#'
#' @exportS3Method sandwich::bread
# fisher.matrix = TRUE means we use the strict definition of 'bread': expected information
# adjust = TRUE means that we account for the Jeffrey penalty
bread.msbm <- function (x, coefficients = x$coefficients,
                        fisher.matrix = TRUE, adjust = TRUE, ...) {

  if(missing(coefficients))
    covmat <- msbm_varcov (x,
                           fisher.matrix = fisher.matrix,
                           adjust = adjust, sigma = 1, ...)
  else
    covmat <- msbm_varcov (x, coefficients = coefficients,
                           fisher.matrix = fisher.matrix,
                           adjust = adjust, sigma = 1, ...)

  covmat <- covmat * x$dims$nobs

  dimnames(covmat) <- list(names(x$coefficients),
                           names(x$coefficients))

  return(covmat)
}

#' @rdname bread.msbm
#' @exportS3Method sandwich::estfun
estfun.msbm <- function (x, coefficients = x$coefficients,
                         adjust = TRUE, ...) {

  if (!missing(coefficients)) {
    okay <- all(x$coefficients == coefficients)
  }
  else {
    okay <- TRUE
  }

  adjust <- adjust & !identical(x$control$criterion, "ML")
  if (okay & !adjust) {
    scores <- attr(x$fisher, "score.i")
  }
  else {
    if (!adjust) {
      scores <- mutate.params (coefficients,
                               frame = model.frame.msbm (x),
                               link = link(x),
                               score = TRUE,
                               information = FALSE,
                               observed = FALSE)

      scores <- scores$score.i
    }
    else {
      frame <- x$frame %||% model.frame(x)
      dev_correction <- function (theta) {
        mutatedtheta <- mutate.params (theta, frame = frame,
                                       link = x$link,
                                       information = TRUE,
                                       observed = FALSE)

        # Penalty is the sum for all observations
        neg_devvalue <- 0.5 * determinant(mutatedtheta$EFinfo, logarithm = TRUE)$modulus[1]

        return(neg_devvalue)
      }

      scores <- numDerivjacobian (func = dev_correction, x = coefficients)
      scores <- c(scores) / x$dims$nobs

      mutation <- mutate.params (coefficients,
                                 frame = model.frame.msbm (x),
                                 link = link(x),
                                 score = TRUE,
                                 information = FALSE,
                                 observed = FALSE)

      scores <- t(t(mutation$score.i) + scores)
    }
  }

  colnames(scores) <- names(x$coefficients)

  return(scores)
}

# Code adapted from 'sandwich::meat'
meat.msbm <- function (x, coefficients = x$coefficients,
                       adjust = TRUE, ...) {
  if (!is.null(x$na.action))
    class(x$na.action) <- "omit"

  if (missing(coefficients)) {
    rval <- crossprod(estfun.msbm(x, adjust = adjust, ...)) / x$dims$nobs
  }
  else {
    rval <- crossprod(estfun.msbm(x, adjust = adjust,
                                  coefficients = coefficients, ...)) / x$dims$nobs
  }

  if (adjust)
    rval <- (x$dims$nobs/(x$dims$nobs - x$dims$npars)) * rval

  dimnames(rval) <- list(names(x$coefficients),
                         names(x$coefficients))

  rval
}
