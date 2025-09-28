# DEBUG NAs handling
#
#' Detect Separation
#'
#' Test for data separation in a generalized linear regression model.
# Poisson or binomial
#'
#' @param object an object of a class with a \code{test.separation} method,
#' or for the default metjod, a numeric vector representing the dependent
#' variable in a Generalized Linear Model (GLM) framework.
#'
#' @param data,weights,subset,na.action,offset,contrasts same as for
#' \link[stats]{glm}.
#'
#' Currently, the only available \code{method} does not handle missing values.
#' Accordingly, the argument \code{na.action} should always be left as
#' \code{na.action = na.omit}.
#'
#' Argument \code{weights} is currently used only when \code{binary = TRUE}.
#' Otherwise, it is ignored by the currently available method.
#' Argument \code{offset} is currently not used at all.
#'
#' @param epsilon small positive real in the open \eqn{(0, 1)}, tolerance for
#' the convergence of residuals to zero. Defaults to \eqn{10^{-4}}.
#'
#' @param method character indicating the separation detection algorithm to
#' be used. Currently, only \code{"iterative_rectifier"} (also \code{'ir'})
#' is available.
#'
#' @param maxit integer, maximum number of iterations allowed.
#'
#' @param x,y numeric matrix/vector, i.e. \code{x} is a design
#' matrix of dimension \code{n * p}; and \code{y} is a vector of observations
#' of length \code{n}.
#'
#' @param binary logical, should the response in \code{formula} or \code{y}
#' be considered as a binomial response? Defaults to \code{TRUE}.
#'
#' @param sanitize logical, should columns corresponding to aliased
#' coefficients and rows corresponding to missing observations (with
#' \link{NA}s) be removed from the design matrix? Defaults to \code{TRUE}.
#'
#' @param ... further arguments passed to or from other methods.
#' Currently not used.
#'
#' @usage
#'
#' test.separation (object, ...)
#'
#' ## S3 method for class 'formula'
#' test.separation (object, binary = TRUE, data = environment(formula(object)),
#'                  weights, subset, na.action = stats::na.omit, offset = NULL,
#'                  contrasts = NULL, sanitize = TRUE, epsilon = 1e-04,
#'                  method = "iterative_rectifier", maxit = 500, ...)
#'
#' ## S3 method for class 'glm'
#' test.separation (object, data = object$data, subset = NULL,
#'                  na.action = stats::na.omit, sanitize = TRUE,
#'                  epsilon = 1e-04, method = "iterative_rectifier",
#'                  maxit = 500, ...)
#'
#' ## Routine for prepared response vector 'y' and design matrix 'x'
#' test.sep (y, x, weights = 1, offset = 0, binary = TRUE, epsilon = 1e-08,
#'           method = "iterative_rectifier", maxit = 500, ...)
#'
#' @details
#' \code{test.separation} is a generic function with a default method that
#' handles a 'numeric' vector \code{object} supplied with a numeric matrix \code{x}.
#'
#' Complete (or quasi) separation occurs in generalized linear models when one or
#' more of the independent variables can perfectly predict all (or some) outcomes
#' \insertCite{albert1984existence}{msbreg}. There a few routines for separation
#' detection in \code{R} (see e.g. \link[detectseparation]{detect_separation}),
#' generally based on linear programming algorithms
#' \insertCite{konis2007linear,albert1984existence}{msbreg}.
#' So far, the only \code{method} implemented here uses the \code{Iterative Rectifier}
#' algorithm of \insertCite{zylkin2019verifying;textual}{msbreg} (p. 39).
#' This is an iterative least squares algorithm that is very fast even in very
#' high dimensional datasets. As the \link[detectseparation]{detect_separation}
#' function, the \code{test.separation} routine is a *pre-fit* method,
#' that is, it does not (need to) estimate the model to detect separation
#' (presence of infinite maximum likelihood estimates) and identify the separated
#' observations.
#'
#' When \code{binary = TRUE} (the default), the binomial model is assumed (no
#' checking is however performed). In this case, the probit-Poisson transformation is
#' applied to the data before the algorithm is applied.
#'
#' When \code{sanitize = TRUE}, sanitation consists of:
#' \itemize{
#' \item removal of \link{NA}s by applying \link[stats]{na.omit} to the data
#' in the response vector (\code{y}) and the design matrix (\code{x}),
#' after their extraction from the \code{object} argument;
#' \item removal of *aliased* columns by using \link[base]{qr}
#' decomposition to find aliased coefficients (see \link[base]{qr.coef}).
#' As such, separation analysis is conditioned (as it should be) on
#' *non-singularity* when \code{sanitize = TRUE}.
#' }
#'
#' The general user-level routine is \code{test.separation}.
#'
#' For the advanced user with prepared response vector \code{y} and design
#' matrix \code{x}, the \code{test.sep} version may be preferred.
#' The latter is the workhorse of \code{test.separation} for the methods
#' for classes \code{formula} and \code{'glm'}, may
#' be substantially faster in large datasets where repeated calls are needed.
#'
#' @return a logical scalar, indicating separation.
#'
#' The returned value has the following attributes:
#'
#' \item{\code{convergence}}{ a named vector of information from the used
#' algorithm:}
#' \itemize{
#' \item \emph{converged}, a binary scalar, did the algorithm converge?
#' \item \emph{niter}, the number of iterations used by the employed algorithm;
#' }
#'
#' \item{\code{separated}}{ \code{NULL}, or an integer vector indicating
#' which observations (if any) are separated by the response;}
#' \item{\code{aliased}}{ \code{NULL}, or integer vector indicating columns of
#' the design matrix \code{x} with aliased coefficients (these columns are
#' removed before any separation related calculation, leading
#' to a non-singular reduced design matrix, used for the subsequent
#' separation analysis);}
#' \item{\code{resid}}{ a named vector of additional information on the
#' achieved convergence threshold, i.e. the absolute values of iterated
#' weighted least-square residuals:}
#'
#' \itemize{
#' \item \emph{min}, the observed minimum of the absolute residuals;
#' \item \emph{max}, the observed maximum of the absolute residuals
#' (convergence is determined by \code{max < epsilon});
#' }
#'
#' \item{\code{epsilon}}{ the used tolerance argument.}
#'
#' Note that the attribute \code{separated} will be non-null only when
#' the output is \code{TRUE}, i.e., when there is separation (otherwise,
#' no observation is separated, and \code{NULL} is returned).
#' Likewise, the attribute \code{aliased} will be non-null only when
#' there are aliased columns (linear dependency between columns of the
#' design matrix).
#'
#' For the algorithm \code{iterative_rectifier},
#' when the output has attribute \code{convergence[1] = 1}:
#'
#' \itemize{
#' \item a \code{TRUE} output indicates that the response separates the
#' linear predictor for some non-zero vector of regression coefficients;
#' \item conversely, a \code{FALSE} output indicates non-separation, i.e.
#' the maximum likelihood estimates of all regression coefficients are finite.
#' }
#'
#' If the output has attribute \code{convergence[1] = 0}, then the
#' attributes \code{resid} (\emph{min} and \emph{max}) indicate to
#' what extent the returned result is far from convergence.
#'
#' @export test.separation
#' @export test.sep
#' @aliases test.sep
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso See \link[msbreg]{identifiability} to test for the practical
#' identifiability of model parameters from observed data.
#'
#' @examples
#' ### Basic functioning
#' # Generate a dataset (binary response
#' # unrelated to the predictors)
#' set.seed(10)
#' df <- data.frame (y = rbinom(15, 1, 0.5),
#'                   # First continuous predictor
#'                   x1 = rnorm(15),
#'                   # Factor with three levels
#'                   f1 = rep(c('a', 'b', 'c'),
#'                            length.out = 15))
#'
#' # Test for separation
#' test.separation (y ~ x1 + f1, data = df)
#'
#
# ### Example from package 'detectseparation'
# data (endometrial)
#
# GLMfit <- glm(HG ~ NV + PI + EH, family = binomial(), data = endometrial)
# summary(GLMfit)
#
# # The huge estimate and standard error for the predictor 'NV' makes
# # separation a very likely issue. Check that:
# test.separation (GLMfit, epsilon = 1e-3)
#
# # The above output indicate that the response separates at least one of the
# # predictors. Indeed, when HG = 0, then all NV = 0, but for HG = 1,
# # NV is zero or one.
#
# The separated observations are: in rows 3 4, 13 32 57 and 67
# # of the dataset ('endometrial').
#
test.separation <- function(object, ...) {
  UseMethod("test.separation")
}

# @rdname test.separation
test.separation.default <- function(object, ...) {

  if (inherits(object, 'formula')) {
    out <- test.separation.formula (object, ...)
    attr(out, 'call') <- match.call()
    return(out)
  }

  if (inherits(object, 'glm')) {
    out <- test.separation.glm (object, ...)
    attr(out, 'call') <- match.call()
    return(out)
  }

  out <- test.sep (y = object, ...)
  attr(out, 'call') <- match.call()
  return(out)
}

setGeneric(name = "test.separation",
           def = test.separation.default)

# @rdname test.separation
#' @exportS3Method msbreg::test.separation
test.separation.formula <- function (object, binary = TRUE, data = environment(formula),
                                     weights, subset, na.action = stats::na.omit, offset = NULL,
                                     contrasts = NULL, sanitize = TRUE, epsilon = 1e-04,
                                     method = "iterative_rectifier", maxit = 500, ...) {
  stopifnot(epsilon > 0, epsilon < 1)
  stopifnot(length(binary) > 0, is.logical(as.logical(binary)))
  mcall <- match.call()

  # Buld a model frame from the formula
  formula <- formula (object)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  if (missing(na.action))
    mf$na.action <- quote(stats::na.omit)
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())

  # Get the output structure
  eval(sep.model.frame())

  attr(out, 'call') <- mcall

  return(out)
}

# @rdname test.separation
#' @exportS3Method msbreg::test.separation
test.separation.glm <- function (object, data = object$data, subset = NULL,
                                 na.action = stats::na.omit,
                                 sanitize = TRUE, epsilon = 1e-04,
                                 method = "iterative_rectifier", maxit = 500, ...) {

  mcall <- match.call()

  binary <- identical(object$family$family, 'binomial')
  contrasts <- NULL

  # Extract model frame from object
  mf <- model.frame(object, data = data, subset = subset,
                    na.action = na.action)

  # Get the output structure
  eval(sep.model.frame())

  attr(out, 'call') <- mcall

  return(out)
}

sep.model.frame <- function() {
  expression({
    if (any(is.na(mf))) {
      warning("NA detected in the model matrix: the current implementation of 'test.separation' cannot handle missing data")
      stop("use 'na.action = stats::na.omit' to discard any NA value")
    }

    # Extract response, design matrix ...
    Y <- stats::model.response(mf, "any")

    if (is.null(Y)) {
      stop('no response found in model frame')
    }

    if (length(dim(Y)) == 1L) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm))
        names(Y) <- nm
    }
    nobs <- NROW(Y)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights)) {
      if (!is.numeric(weights))
        stop("'weights' must be a numeric vector")
      if (any(weights < 0))
        stop("negative weights not allowed")
      if (length(weights) != nobs)
        stop(gettextf("number of weights is %d should equal %d (number of observations)",
                      length(weights), nobs), domain = NA)
    }
    else {
      weights <- 1
    }

    if (NCOL(Y) == 1) {
      if (is.factor(Y))
        Y <- Y != levels(Y)[1L]
      Y[weights == 0] <- 0
      if (any(na.omit(Y) < 0 | na.omit(Y) > 1))
        stop("y values must be 0 <= y <= 1")
      my <- weights * Y
      if (any(abs(na.omit(my - round(my))) > 0.001))
        warning(gettextf("non-integer #successes in a %s glm!",
                         "binomial"), domain = NA)
      Y <- Y * weights
      ntrials <- weights
    }
    else if (NCOL(Y) == 2) {
      if (any(abs(Y - round(Y)) > 0.001))
        warning(gettextf("non-integer counts in a %s glm!",
                         "binomial"), domain = NA)
      n <- (Y1 <- Y[, 1L]) + Y[, 2L]
      Y <- Y1
      if (any(n0 <- n == 0))
        Y[n0] <- 0
      Y <- Y * weights
      ntrials <- weights * n
    }
    else stop(gettextf("for the '%s' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is #successes and col 2 is #failures",
                       "binomial"), domain = NA)
    if (length(ntrials) > 1)
      names(ntrials) <- names(Y)

    X <- stats::model.matrix(mf,
                             data = data,
                             contrasts.arg = contrasts)
    if (NROW(X) != nobs) {
      if (!is.null(rownames(mf))) {
        keep <- rownames(X) %in% rownames(mf)
        X <- X[keep, , drop = FALSE]
        rownames(X) <- rownames(mf)[keep]
      }
      else {
        dnaaction <- attr(mf, "na.action")
        if (length(dnaaction)) {
          notkeep <- (1:NROW(X)) %in% dnaaction
          X <- X[!notkeep, , drop = FALSE]
        }
      }
    }

    offset <- as.vector(stats::model.offset(mf))
    if (!is.null(offset)) {
      if (length(offset) != nobs)
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(offset), nobs), domain = NA)
    }
    else {
      offset <- 0
    }

    #### Sanitation if required
    # Remove any NA
    n_y <- length(Y)
    if (sanitize) {
      rnames <- names(Y)
      if (is.null(rnames))
        rnames <- names(Y) <- 1:length(Y)
      yx <- na.omit(cbind(Y, ntrials, offset, X))
      dnaaction <- attr(yx, 'na.action')
      if (length(dnaaction)) {
        Y <- yx[,1]
        ntrials <- yx[, 2]
        offset <- yx[, 3]
        X <- yx[, -(1:3), drop = FALSE]
      }
      rm(yx)
    }
    else {
      dnaaction <- integer(0L)
    }

    # Remove aliased columns
    if (sanitize) {
      # Remove any aliased column in 'X'
      QRx <- base::qr(X, tol = min(epsilon, 1e-07))
      coeffs <- qr.coef(QRx, Y)
      keep <- !is.na(coeffs)
      if (!all(keep)) {
        if (any(keep)) {
          X <- X[, keep, drop = FALSE]
        }
        else {
          stop("should not get here (at least one column should remain): something went wrong")
        }
      }
      aliased <- which(!keep)
    }
    else {
      aliased <- integer(0L)
    }

    if (binary) {
      stopifnot(all(Y <= ntrials))
    }

    # Call the workhorse function
    out <- test.sep (y = Y, x = X,
                     weights = ntrials,
                     offset = offset,
                     epsilon = epsilon,
                     binary = binary,
                     method = method[1],
                     maxit = maxit[1])

    if (length(dnaaction)) {
      da.nm <- names(attr(out, "separated"))
      separated <- sapply(attr(out, "separated"),
                          FUN = function (j) {
                            j + sum(dnaaction < j)
                          })
      names(separated) <- da.nm
      attr(out, "separated") <- separated
    }

    if (length(aliased))
      attr(out, "aliased") <- aliased
  })
}
