# For bold text: \strong{.}, alternative to **.**
#' Summarizing a Multistage Binomial Model Fits
#'
#' Compute summary statistics for a multistage binomial (MSB) model fit.
#' This is a method for class \code{"msbm"}.
#'
#' @param object an object of class "msbm",
#' usually, the result of a call to \link[msbreg]{msbreg}.
#
#' @param fisher.matrix logical, should the Fisher (expected) information
#' matrix be used to obtain standard errors? Defaults to \code{TRUE}.
#' The alternative is the observed information matrix.
#'
#' @param adjust logical, should the covariance matrix be adjusted for
#' any penalty added to the log-likelihood function during model fitting?
#' Defaults to \code{FALSE}.
#'
#' @param sandwich logical, should the sandwich variance-covariance
#' estimator be used? Defaults to \code{FALSE}.
#'
#' @param se.values logical, should estimated standard errors of
#' the \code{mu}, \code{alpha}, \code{lambda}, and
#' \code{alpha.lambda} components (of the input \code{object}) be
#' calculated and included in the output? Defaults to \code{FALSE}.
#'
#' @param IC logical, should the fit with intercept-correction be considered.
#'
#' @param correlation,symbolic.cor Same as for \link[stats]{glm}, i.e.
#' \itemize{
#' \item \code{correlation}: a logical, if \code{TRUE}, the correlation matrix
#' of the estimated parameters is returned and printed;
#'
#' \item \code{symbolic.cor}: a logical, if \code{TRUE}, print the correlations
#' in a symbolic form (see \link[stats]{symnum}) rather than as numbers.
#' }
#'
# @param me logical, if \code{TRUE}, the measurement error component
# of the model is summarized and printed.
#'
#' @param rsquared.method character indicating which pseudo \eqn{R^2}
#' statistics to include in the model summary.
#' Passed as \code{method} to \link[msbreg]{rsquared}.
#'
#' @param cor.method character indicating which correlation coefficient
#' to consider (only used when \code{rsquared.method} includes \code{"COR"}).
#' Passed as \code{cor.method} to \link[msbreg]{rsquared}.
#'
#' @param frame logical, if \code{TRUE}, the model frame (component
#' \code{$frame} of the input) is returned/printed.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' This function mimics the similar \code{summary}
#' method for "\link[stats]{glm}" class objects.
#'
#' @return
#' The function returns an object of class \code{"summary.msbm"},
#' a list with components:
#'
# \describe{
#' \item{\code{call}, \code{link}, \code{deviance}, \code{aic},
#' \code{bic}, \code{hqc}, \code{df.residual}, \code{null.deviance},
#' \code{df.null}, \code{iter}, \code{converged}, \code{boundary},
#' \code{na.action}, \code{deviance.resid}, \code{y.resid},
#' \code{dispersion}, \code{dims}, \code{me}, \code{me.offset},
#' \code{me.sd} and \code{frame}}{ the corresponding components from
#' \code{object};}
#'
#' \item{\code{coefficients}}{ the matrix of coefficients,
#' standard errors, \code{z}-values and \code{p}-values
#' (aliased coefficients are omitted in the returned output,
#' but printed);}
#'
#' \item{\code{aliased}}{ named logical vector showing if the original coefficients are aliased;}
#'
#' \item{\code{df}}{ a 3-vector of the rank of the model and the number of
#' residual degrees of freedom, plus number of coefficients (including aliased ones);}
#'
#' \item{\code{cov.scaled}}{ the (scaled) estimate of the covariance matrix of
#' the estimated coefficients;}
#'
#' \item{\code{correlation}}{ if requested (\code{correlation = TRUE}),
#' the estimate of the matrix of correlations between the estimated
#' coefficients;}
#'
#' \item{\code{symbolic.cor}}{ if \code{correlation = TRUE},
#' the value of the argument \code{symbolic.cor};}
#'
#' \item{\code{mu}}{ numeric, the \code{$fitted.values} component of
#' \code{object}, that is, the fitted individual success probabilities
#' \eqn{\mu} of the binary responses;}
#'
#' \item{\code{alpha}, \code{lambda}}{ numeric, the \eqn{\alpha} and
#' \eqn{\lambda} components (fitted or known values/offsets) of the model:
#' the vector of individual \eqn{\alpha} or \eqn{\lambda} values is
#' reduced to a scalar if constant across all individuals;}
#'
#' \item{\code{se.alpha}, \code{se.lambda}}{ (vector of) estimated standard
#' errors of \eqn{\alpha} and \eqn{\lambda} value(s);}
#'
#' \item{\code{alpha.lambda}}{ numeric, product of \eqn{\alpha} and
#' \eqn{\lambda} components of the model: this gives the fitted
#' theoretical minimum success probability of the binary response
#' in the model;}
#'
#' \item{\code{se.alpha.lambda}}{ (vector of) estimated standard
#' errors of \code{alpha.lambda} value(s);}
#'
#' \item{\code{r.squared}}{ a numeric vector of the same length as
#' \code{rsquared.method} with pseudo-\eqn{R^2} values returned by
#' \link[msbreg]{rsquared};}
#'
#' \item{\code{criterion}}{ the \code{$control$criterion} component of
#' \code{object};}
#'
#' \item{\code{alpha.lambda.deficit}}{ logical, does the fit rank suggest that
#' the \eqn{\alpha} and/or \eqn{\lambda} components are over-parametrized?
#' Specifically, for a rank deficient information matrix, if the rank
#' deficit corresponds to the number of parameters in the \eqn{\alpha} and
#' \eqn{\lambda} components, these parameters are ignored and the rank of
#' the information for the remaining parameters checked: if the reduced
#' information is of full rank, then \code{alpha.lambda.deficit = TRUE}
#' and the variance components for \eqn{\alpha} and/or \eqn{\lambda} are
#' set to zero. Otherwise, \code{alpha.lambda.deficit = FALSE}.}
# }
#'
#' @seealso See \link[msbreg]{print.msbm} for printing a quick summary.
#' Also see \link[msbreg]{msbreg} for model fitting.
#'
#' @aliases summary.msbm
# @exportS3Method base::print
# @exportS3Method print msbm
# @export print.msbm
#'
#' @examples
#' data("test1data")
#'
#' MSBfit = msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                  alpha.formula = ~ 0,
#'                  lambda.formula = ~ 0,
#'                  weights = Total,
#'                  data = test1data)
#' MSBfit # Print coefficient estimates (descriptive)
#'
#' summary (MSBfit) # Print a more detailed summary of the fit (inference)
#'
#' @exportS3Method base::summary
summary.msbm <- function (object,
                          fisher.matrix = TRUE,
                          adjust = FALSE,
                          sandwich = FALSE,
                          se.values = FALSE,
                          IC = NULL,
                          correlation = FALSE,
                          symbolic.cor = FALSE,
                          rsquared.method = 'KL',
                          cor.method = "pearson",
                          frame = TRUE,
                          ...) {
  if (!is.null(IC)) {
    IC <- as.logical(IC)
    if (!IC & object$IC) {
      object <- update(object, action = "cancel.IC")
    }
    else if (IC & !object$IC) {
      object <- update(object, action = "apply.IC")
    }
  }

  if (is.null(object$frame)) {
    object$frame <-  model.frame.msbm(object)
  }

  coefs <- object$coefficients
  p <- object$dims$p
  d <- object$dims$d
  r <- object$dims$r
  npars <- p + d + r
  aliased <- is.na(coefs)
  mrank <- object$rank
  covmat <- vcov.msbm (object, fisher.matrix = fisher.matrix,
                       adjust = adjust, sigma = 1,
                       sandwich = sandwich)
  alpha.lambda.deficit <- attr(covmat, "alpha.lambda.deficit") %||% FALSE

  dn <- c("Estimate", "Std. Error")
  if (any(is.na(covmat)) |
      any(class(covmat) %in% c("simpleError", "error",
                               "condition", "try-error"))) {
    # What if (mrank < NCOL(object$fisher)) ?

    # Elaborate on that here

    warning("algorithm did not converge: no covariance matrix estimate available")
    covmat <- NA # matrix(NA, ncol = mrank, nrow = mrank)
    coef.table <- cbind(coefs, NaN, NaN, NaN)
    dimnames(coef.table) <- list(names(coefs),
                                 c(dn, "z value", "Pr(>|z|)"))

    mu.var <- NA

    if (d) {
      alpha.var <- NA
    }

    if (r) {
      lambda.var <- NA
    }

    if (d | r) {
      alpha.lambda.var <- NA
    }
  }
  else {

    dimnames(covmat) <- list(names(coefs), names(coefs))
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    zvalue <- coefs/s.err
    pvalue <- 2 * pnorm(-abs(zvalue))
    coef.table <- cbind(coefs, s.err, zvalue, pvalue)
    dimnames(coef.table) <- list(names(coefs),
                                 c(dn, "z value", "Pr(>|z|)"))

    if (se.values) {
      # Standard errors for mu
      if (p > 0) {
        dmu <- attr(object$fisher, "mu.coefficients")
        mu.var <- rowSums ((dmu %*% covmat) * dmu)
      }

      # Standard errors for alpha (from gamma parameters, if any) by the Delta method
      if (d) {
        alpha.index <- (p + 1):(p + d)
        if ((d - attr(object$frame$alpha.dictionary, "intercept")) == 0 &
            !attr(object$frame$alpha.dictionary, "offset")) {
          c.alpha <- attr(object$fisher, "c.alpha")[1]
          dalpha <- attr(object$fisher, "mu.coefficients")[1, alpha.index]
          dalpha <- dalpha / c.alpha
          alpha.var <- (dalpha^2) * covmat[alpha.index, alpha.index]
        }
        else {
          c.alpha <- attr(object$fisher, "c.alpha")
          dalpha <- attr(object$fisher, "mu.coefficients")[, alpha.index, drop = FALSE]
          dalpha <- dalpha / c.alpha
          alpha.var <- rowSums ((dalpha %*%
                                   covmat[alpha.index, alpha.index]) *
                                  dalpha)
        }
      }

      # Standard errors for lambda (from delta parameter, if any)
      if (r) {
        lambda.index <- (p + d + 1):npars

        if ((r - attr(object$frame$lambda.dictionary, "intercept")) == 0 &
            !attr(object$frame$lambda.dictionary, "offset")) {
          c.lambda <- attr(object$fisher, "c.lambda")[1]
          dlambda <- attr(object$fisher, "mu.coefficients")[1, lambda.index]
          dlambda <- dlambda / c.lambda
          lambda.var <- (dlambda^2) * covmat[lambda.index, lambda.index]
        }
        else {
          c.lambda <- attr(object$fisher, "c.lambda")
          dlambda <- attr(object$fisher, "mu.coefficients")[, lambda.index, drop = FALSE]
          dlambda <- dlambda / c.lambda
          lambda.var <- rowSums ((dlambda %*%
                                    covmat[lambda.index, lambda.index]) *
                                   dlambda)
        }
      }

      # Standard errors for alpha * lambda ()
      if (d | r) {
        if (d * r) {
          dalpha.lambda <- cbind(dalpha * object$lambda.values, dlambda * object$alpha.values)
          alpha.lambda.var <- rowSums ((dalpha.lambda %*%
                                          covmat[c(alpha.index, lambda.index),
                                                 c(alpha.index, lambda.index)]) *
                                         dalpha.lambda)
        }
        else if (d) {
          alpha.lambda.var <- alpha.var * (object$lambda.values^2)
        }
        else if (r) {
          alpha.lambda.var <- lambda.var * (object$alpha.values^2)
        }
      }
    }
  }

  keep <- match(c("call", "link", "deviance", "aic", "bic", "hqc",
                  "df.residual", "null.deviance", "df.null", "iter",
                  "converged", "boundary", "na.action",
                  "deviance.resid", "y.resid", "dispersion", "dims",
                  "me", "me.offset", "me.sd",
                  if (frame) "frame"),
                names(object), 0L)

  if (length(rsquared.method)) {
    r.squared <- as.numeric(rsquared(object, method = rsquared.method,
                                     cor.method = cor.method,
                                     adjust.size = FALSE)[,1:length(rsquared.method)])
    names(r.squared) <- rsquared.method
  }
  else {
    r.squared <- as.numeric(rsquared(object, method = "KL", adjust.size = FALSE)[1,])
    names(r.squared) <- "KL"
  }

  ans <- c(object[keep],
           list(coefficients = coef.table,
                aliased = aliased,
                df = c(rank = mrank,               # Fisher matrix rank
                       df.r = object$df.residual,  # Residual df
                       df.f = length(aliased)),    # Fit df
                cov.scaled = covmat,
                mu = object$fitted.values, se.mu = if (se.values) sqrt(mu.var),
                alpha.lambda.deficit = alpha.lambda.deficit,
                alpha = object$alpha.values, se.alpha = if (d & se.values) sqrt(alpha.var),
                lambda = object$lambda.values, se.lambda = if (r & se.values) sqrt(lambda.var),
                alpha.lambda = object$alpha.values * object$lambda.values,
                se.alpha.lambda = if ((d | r) & se.values) sqrt(alpha.lambda.var),
                r.squared = r.squared,
                criterion = object$control$criterion))

  # Add squared-correlation as a primary R^2 statistic?
  # Can use Pearson, Spearman, Kendal correlation,
  # or the Mutual Information Coefficient:
  # R(x,y) = [1 - exp(-2*MI(x,y))]^(0.5)
  # See also FNN::mutinfo(, k = 2) for Kraskov’s nearest neighbor estimator

  if (correlation && mrank > 0) {
    dd <- sqrt(diag(covmat))
    ans$correlation <- covmat/outer(dd, dd)
    ans$symbolic.cor <- symbolic.cor
  }

  class(ans) <- "summary.msbm"
  return(ans)
}

#setMethod ("summary",
#           signature = "msbm",
#           definition = summary.msbm)
