#' Variance-Covariance Matrix for \code{'msbm'} Objects
#'
#' Estimate the variance-covariance matrix of a multistage
#' binomial (MSB) model fit.
#'
#' @param object an object of class \code{'msbm'} (a fitted MSB model as
#' returned by \link[msbreg]{msbreg}) or class \code{'summary.msbm'} (as
#' returned by \link[msbreg]{summary.msbm}).
#'
#' @param coefficients a numeric vector of all MSB model parameters.
#' When supplied, \code{coefficients} must have finite value,
#' and the same length as \code{object$coefficients}.
#'
#' @param fisher.matrix logical, should the Fisher (expected) information
#' matrix be used to obtain the covariance matrix? Defaults to \code{TRUE}.
#' The alternative is the observed information matrix.
#'
#' @param sigma a numeric scalar (positive real), optional
#' dispersion parameter, defaults to \code{sigma = 1}.
#'
#' @param sandwich logical, should the sandwich estimator be used?
#' Defaults to \code{FALSE}.
#'
#' @param adjust logical, should the covariance matrix be adjusted for
#' any penalty added to the log-likelihood function during model fitting?
#' Only relevant when \code{fisher.matrix = FALSE} (Fisher information matrix
#' is independent of any penalty not involving the response).
#' Defaults to \code{FALSE}.
#'
#' When \code{sandwich = TRUE}, \code{adjust = TRUE} also means a finite sample
#' adjustment made by multiplying the result by \eqn{n/(n−k)} where \eqn{n} is
#' the number of observations and \eqn{k} the number of estimated parameters.
#'
#' @param complete logical indicating if the full variance-covariance
#' matrix should be returned even when some coefficients are aliased (undefined)
#' and \link[stats]{coef}\code{(.)} contains \code{NA}s correspondingly.
#' Using \code{complete = TRUE} ensures that \code{vcov(.)} is compatible
#' with \code{coef(.)}. Setting \code{complete = FALSE} can be useful
#' to drop rows and columns corresponding to aliased coefficients.
#'
#' @param ... further arguments passed to or from other methods.
#' None is currently used.
#'
#' @details
#' These functions mimic and extend the \link[stats]{vcov} method for
#' "\link[stats]{glm}" class objects to handle \code{'msbm'} objects.
#' The \code{vcov} method returns a variance-covariance matrix
#' estimate, by default, the inverse of the information matrix.
#'
#' For an object of class \code{'msbm'}, the function allows to substitute
#' the vector of estimated parameters by a user-supplied vector
#' (\code{coefficients}). It also allows the use of an observed instead
#' of expected information, and adjusts, if desired, for any  penalty
#' added to the log-likelihood function during model fitting.
#'
#' When \code{sandwich = TRUE}, the sandwich covariance matrix is
#' constructed by multiplying **bread** and **meat** matrices
#' \insertCite{zeileis2006object}{msbreg}.
#' The *bread* matrix is the inverse information matrix times the
#' number of observations, and the *meat* matrix is the average cross
#' product of the empirical estimating functions (matrix of scores, that
#' is, first order  derivatives of the log-likelihood contributions
#' with respect to model parameters, evaluated at \code{coefficients}).
#' The implementation follows \insertCite{zeileis2020various;textual}{msbreg}.
#' See \link[msbreg]{bread.meat} for extracting these elementary matrices.
#'
#' @return A matrix of the estimated variances-covariances for the parameter
#' in \code{coefficients}.
#' This has row and column names corresponding to the parameter names
#' given by the \link[stats]{coef} method.
#'
#' The returned matrix has the attribute \code{alpha.lambda.deficit}, a
#' logical indicating whether the information matrix is rank deficient
#' and the rank deficit is due to the \eqn{\alpha} and/or \eqn{\lambda}
#' component(s) of the MSB model.
#' If \code{alpha.lambda.deficit = TRUE}, the rows and columns corresponding
#' to the \eqn{\alpha} and/or \eqn{\lambda} component(s) of the model fit
#' in the returned matrix are all set to zero.
#'
#' @references
#' \insertAllCited{}
#'
#' @aliases varcov
#' @importFrom Rdpack reprompt
#'
#' @seealso See \link[msbreg]{bread.meat} for \code{"msbm"} class methods for
#' computing empirical estimating functions and their covariance matrix.
#'
#' @examples
#' data("test1data")
#'
#' MSBfit = msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                  weights = Total,
#'                  data = test1data)
#' summary(MSBfit)
#'
#' sqrt (diag (vcov (MSBfit)))
#'
#' round(vcov (MSBfit), 5)
#'
#' round(vcov (MSBfit, sandwich = TRUE), 5)
#'
#' @exportS3Method stats::vcov
vcov.msbm <- function (object, coefficients = object$coefficients,
                       fisher.matrix = TRUE, adjust = FALSE,
                       sigma = NULL, sandwich = FALSE,
                       complete = TRUE, ...) {

  covmat <- msbm_varcov (object, coefficients = coefficients,
                         fisher.matrix = fisher.matrix, adjust = adjust,
                         sigma = sigma, ...)

  alpha.lambda.deficit <- attr(covmat, "alpha.lambda.deficit") %||% FALSE

  if (!any(is.na(covmat)) & sandwich) {
    vmeat <- meat.msbm (object, coefficients = coefficients, adjust = adjust, ...)
    covmat <- covmat %*% vmeat %*% covmat
    covmat <- covmat * object$dims$nobs
  }
  else if (any(is.na(covmat))) {
    warning(paste0(covmat[1], ": returning NAs"))
    covmat <- matrix(NA, nrow = object$dims$npars,
                     ncol = object$dims$npars)
  }

  dimnames(covmat) <- list(names(object$coefficients),
                           names(object$coefficients))

  if (!complete & any((aliased <- is.na(coefficients)))) {
    covmat <- covmat[!aliased, !aliased]
  }

  return(structure(covmat,
                   alpha.lambda.deficit = alpha.lambda.deficit,
                   class = c('vcov', class(covmat))))
}

# Basic 'vcov'
msbm_varcov <- function (object, coefficients = object$coefficients,
                        fisher.matrix = TRUE, adjust = FALSE,
                        sigma = NULL, try_empirical = FALSE, ...) {

  # Default sigma if required
  eval(toget.msbm.sigma())

  # If the Fisher matrix is already there
  okay <- !is.null(object$fisher)
  if (!missing(coefficients)) {
    okay <- okay & all(object$coefficients == coefficients)
  }

  adjust <- adjust & !identical(object$control$criterion, "ML")
  if (okay & !adjust) {
    covmat <- if (fisher.matrix)
      object$fisher
    else
      attr(object$fisher, "observed")
    covmat <- attr(attr(covmat, "qr"), "solve")

    if (any(is.na(covmat)) &
        object$rank == object$dims$npars - object$dims$d - object$dims$r) {
      rindex <- 1:object$rank
      covmatr <- if (fisher.matrix)
        object$fisher
      else
        attr(object$fisher, "observed")
      covmatr <- safe_solve(covmatr[rindex, rindex])

      if (!any(class(covmatr) %in% c("simpleError", "error",
                                     "condition", "try-error"))) {
        covmat <- matrix(0, nrow = object$dims$npars,
                         ncol = object$dims$npars)
        covmat[rindex, rindex] <- covmatr
        attr(covmat, "alpha.lambda.deficit") <- TRUE
      }
    }

    if (any(is.na(covmat))) {
      if (try_empirical) {
        vmeat <- meat.msbm(object, adjust = adjust)
        if (is.matrix(vmeat)) {
          covmat <- safe_solve(vmeat * object$dims$nobs)

          if (any(class(covmat) %in% c("simpleError", "error",
                                       "condition", "try-error"))) {
            covmat <- NA
          }
        }
      }

      if (any(is.na(covmat))) {
        covmat <- matrix(NA, nrow = object$dims$npars,
                         ncol = object$dims$npars)
      }
    }
  }
  else {
    if (!adjust) {
      mutation <- mutate.params (coefficients,
                                 frame = model.frame.msbm (object),
                                 link = link(object),
                                 information = TRUE,
                                 observed = !fisher.matrix)

      covmat <- catch.conditions({
        if (fisher.matrix) solve(mutation$EFinfo)
        else solve(mutation$OFinfo)
      })$value

      if (any(class(covmat) %in% c("simpleError", "error",
                                   "condition", "try-error")) &
          object$rank == object$dims$npars - object$dims$d - object$dims$r) {
        rindex <- 1:object$rank
        covmatr <- if (fisher.matrix)
          mutation$EFinfo
        else
          mutation$OFinfo
        covmatr <- safe_solve(covmatr[rindex, rindex])

        if (!any(class(covmatr) %in% c("simpleError", "error",
                                       "condition", "try-error"))) {
          covmat <- matrix(0, nrow = object$dims$npars,
                           ncol = object$dims$npars)
          covmat[rindex, rindex] <- covmatr
          attr(covmat, "alpha.lambda.deficit") <- TRUE
        }
      }

      if (any(class(covmat) %in% c("simpleError", "error",
                                   "condition", "try-error"))) {

        warning(paste0(covmat[1], ": returning NAs"))
        covmat <- matrix(NA, nrow = object$dims$npars,
                         ncol = object$dims$npars)
      }
    }
    else if (okay & !fisher.matrix) {
      covmat <- attr(object$criterion, "hessian")

      covmat <- catch.conditions({
        solve(covmat)
      })$value

      if (any(class(covmat) %in% c("simpleError", "error",
                                   "condition", "try-error")) &
          object$rank == object$dims$npars - object$dims$d - object$dims$r) {
        rindex <- 1:object$rank
        covmatr <- attr(object$criterion, "hessian")
        covmatr <- safe_solve(covmatr[rindex, rindex])

        if (!any(class(covmatr) %in% c("simpleError", "error",
                                       "condition", "try-error"))) {
          covmat <- matrix(0, nrow = object$dims$npars,
                           ncol = object$dims$npars)
          covmat[rindex, rindex] <- covmatr
          attr(covmat, "alpha.lambda.deficit") <- TRUE
        }
      }
    }
    else {
      if (fisher.matrix) {
        infomat <- mutation$EFinfo
      }
      else {
        frame <- object$frame %||% model.frame(object)
        dev_correction <- function (theta) {
          mutatedtheta <- mutate.params (theta, frame = frame,
                                         link = object$link,
                                         information = TRUE,
                                         observed = FALSE)

          neg_devvalue <- 0.5 * determinant(mutatedtheta$EFinfo, logarithm = TRUE)$modulus[1]

          return(-neg_devvalue)
        }

      # Penalty related correction to the information matrix (always observed)
      infomat <- numDerivhessian (func = dev_correction, x = coefficients)

      mutation <- mutate.params (coefficients,
                                 frame = model.frame.msbm (object),
                                 link = link(object),
                                 information = TRUE,
                                 observed = !fisher.matrix)

      # Using observed penalty information, even if fisher.matrix = TRUE
      infomat <- infomat + mutation$OFinfo
      }

      covmat <- catch.conditions({
        solve(infomat)
      })$value

      if (any(class(covmat) %in% c("simpleError", "error",
                                   "condition", "try-error")) &
          object$rank == object$dims$npars - object$dims$d - object$dims$r) {
        rindex <- 1:object$rank
        covmatr <- safe_solve(infomat[rindex, rindex])

        if (!any(class(covmatr) %in% c("simpleError", "error",
                                       "condition", "try-error"))) {
          covmat <- matrix(0, nrow = object$dims$npars,
                           ncol = object$dims$npars)
          covmat[rindex, rindex] <- covmatr
          attr(covmat, "alpha.lambda.deficit") <- TRUE
        }
      }
    }
  }

  if (any(class(covmat) %in% c("simpleError", "error",
                                      "condition", "try-error"))) {

    warning(paste0(covmat[1], ": returning NAs"))
    covmat <- matrix(NA, nrow = object$dims$npars,
                     ncol = object$dims$npars)
  }
  else if (object$dispersion != 1) {
    covmat <- covmat / object$dispersion^2
  }

  covmat <- (sigma^2) * covmat

  dimnames(covmat) <- list(names(object$coefficients),
                           names(object$coefficients))

  attr(covmat, "alpha.lambda.deficit") <- attr(covmat, "alpha.lambda.deficit") %||% FALSE

  return(covmat)
}

#'
#' @rdname vcov.msbm
#' @exportS3Method stats::vcov
vcov.summary.msbm <- function (object, complete = TRUE, ...)
  .vcov.aliased(object$aliased, object$cov.scaled, complete = complete)

#' @exportS3Method base::print
print.vcov <- function(x, n = Inf,
                              digits = getOption("digits") - 3, ...) {
  out <- as.data.frame(as.matrix(x))
  if (is.finite(n)) {
    n <- min(n, NCOL(out))
    out <- out[1:n, 1:n]
  }

  print(out, digits = digits, ...)

  if (!is.null(attr(x, "sigma"))) {
    cat(" *Dispersion (sigma): ",
        round(attr(x, "sigma"), digits, ...),
        "\n")
  }
}
