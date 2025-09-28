#' Calculate Pseudo \eqn{R^2}
#'
#' Calculate generalized \eqn{R^2} statistics to assess the goodness of fit
#' (explanative power) of generalized (non)linear models.
#' EXTEND TO COVER ALL FIT RESULTS WITH an 'nobs' component, a 'rank' component
#' (or 'df.null' and 'df.residual' components), and a 'deviance' and a
#' 'deviance.null' components; a 'y' and a 'fitted' component for
#' correlation based measures.
#'
#' @aliases rsquared
#' @aliases rsquared.default
#' @aliases NagelkerkeR2
#'
#' @param object an object of (or inheriting from) a class with a
#' \code{rsquared} method. The default method handles objects of
#' class \code{"glm"}, or a list of such objects.
#'
#' @param ... optional further arguments to be passed to or from other methods.
#' Generally, more fitted model objects of a class with a \code{rsquared}
#' method.
#'
# For a \code{'glm'} object, the default \code{rsquared} method accepts
# any of the following arguments.
#'
#' @param method character indicating the desired pseudo R squared statistic(s).
#' Any of the following alternatives:
#'
#' \itemize{
#' \item \code{'KL'} for Kullback–Leibler divergence based
#' \eqn{R^2} (default), equivalent to the deviance reduction ratio achieved by the predictors
#' in a fitted model as compared to a model with no predictor;
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
#' @export rsquared
#' @export rsquared.default
#' @import Rdpack
#'
#' @details
#' Arguments after \code{...} must be matched exactly (by their names).
#'
#' The function \code{rsquared} computes generalized R-squared statistics
#' such as the deviance reduction ratio or the Nagelkerke's \eqn{R^2}.
#' This is a generic function for which methods can be added to deal
#' with specific class objects.
#'
#' The default method for \code{rsquared} handles objects of classes
#' \code{'glm'}. It computes by default (\code{method = 'KL'})
#' the deviance based pseudo-\eqn{R^2} of \insertCite{cameron1996r;textual}{msbreg}.
#' The alternative \code{method = 'Nagelkerke'} computes
#' \insertCite{nagelkerke1991note;textual}{msbreg}'s pseudo R-squared.
#'
#' The function \code{NagelkerkeR2} is an alias to compute Nagelkerke's
#' R-squared. It extends \link[fmsb]{NagelkerkeR2} of library \code{fmsb} to
#' handle objects of any class with an \code{rsquared} method (e.g.
#' \code{'msbm'}) in addition to \code{'glm'} class objects.
#'
#' @note
#' The original function \link[fmsb]{NagelkerkeR2} was written by Minato Nakazawa.
#'
#' @return A \code{'data.frame'} with rows corresponding to the supplied
#' \code{objects} (including \code{...} if any), and a column of \code{R2} for
#' each provided method, and two additional columns including \code{N} (number
#' of observations), \code{df} (number of free model parameters).
#' As such, \code{NagelkerkeR2} always returns a three-column
#' \code{data.frame}.
#'
#' @seealso \link[stats]{AIC}, \link[stats]{BIC} and \link[msbreg]{HQC} for
#' using information criteria to assess model fits.
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' ## Logistic regression fit
#' GLMres <- glm (case ~ spontaneous,
#'                data = infert,
#'                family = binomial())
#'
#' summary(GLMres)
#'
#' # Kullback–Leibler divergence based R^2 (the default)
#' rsquared (GLMres)
#'
#' ## MRB model fit with a maximum probability limit
#' require(msbreg)
#' MSBres <- msbreg (case ~ spontaneous,
#'                   lambda.formula = ~ 1,
#'                   data = infert,
#'                   start = c(as.numeric(GLMres$coefficients), 4),
#'                   criterion = 'ML')
#'
#' summary(MSBres)
#'
#' # Nagelkerke's R^2 for the two fits
#' NagelkerkeR2 (GLMres, MSBres)
#'
#' # Nagelkerke's R^2 and Kullback–Leibler divergence based R^2 for the two fits
#' rsquared (GLMres, MSBres, method = c('Nagelkerke', 'KL'))
#'
rsquared <- function(object, ...) {
  UseMethod("rsquared")
}

rsquared.default <- function(object, ...) {
  # Save the call (for naming purposes)
  Call <- match.call()
  out <- rsquared_default_core(object, ...)

  rownames(out) <- as.character(Call[-1L])[1:NROW(out)]

  return(out)
}

setGeneric(name = 'rsquared',
           def = rsquared.default)

rsquared_default_core <-
  function (object, ...,
            method = 'KL',
            cor.method = "pearson",
            adjust.size = FALSE) {
    # Save the call (for naming purposes)
    Call = match.call()

    # If many methods, return a named list with one call to 'rsquared' per supplied method
    method <- unique(method)
    if (length(method) > 1) {
      Lresults <- lapply(method,
                         FUN = function(method_i) {
                           out <- rsquared.default (object, ..., method = method_i, adjust.size = adjust.size)
                           return(out)
                         })
      out <- Lresults[[1]]
      out <- cbind(as.numeric(out[,1]),
                   do.call("cbind",
                           lapply(Lresults[-1], FUN = function (outi) {
                             as.numeric(outi[,1])
                           })),
                   as.numeric(out[,2]),
                   as.numeric(out[,3]))

      method[method %in% c('n', 'N')] <- 'Nagelkerke'
      colnames(out) <- c(toupper(method), 'N', 'df')
      rownames(out) <- as.character(Call[-1L])[1:NROW(out)]

      return(out)
    }

    # Select the function for a specified method
    get.R2VAL <- switch(toupper(method),
                        COR = function(object, n, rank, adjust.size) {
                          get.COR2VAL (object = object, n = n, rank = rank,
                                       adjust.size = adjust.size,
                                       method = cor.method)
                        },
                        KL = get.KLR2VAL,
                        NAGELKERKE = get.NR2VAL,
                        N = get.NR2VAL)

    # List of classes with a method
    availclasses <- c("glm", "msbm", "glm2msbm")

    # Argument adjust.size: coerce to be logical
    adjust.size <- as.logical(adjust.size[1])

    # List of additional object, if any
    objects <- list(...)

    if (length(objects) == 0) { # If no additional object provided through  ...

      if (any(class(object) %in% availclasses)) { # If object is of a class with an available method

        if (inherits(object, "glm")) {
          n <- nrow(object$model)
          rank <- object$df.null - object$df.residual + 1
        }
        else {
          n <- object$nobs
          rank <- object$rank
        }

        R2 <- get.R2VAL (object = object, n = n, rank = rank, adjust.size = adjust.size)

        RVAL <- data.frame(rbind(c(R2 = R2, N = n, df = rank)))
        colnames(RVAL) <- c(method, 'N', 'df')
        rownames(RVAL) <- make.unique(as.character(Call[-1L])[1])

        return(RVAL)
      }
      else {
        if (is.list(object)) {
          nm <- names(object)

          clobject <- sapply(object,
                             FUN = function(x) any(class(x) %in% availclasses))
          if (all(!clobject)) {
            cat("\n classes with a default 'rsquared' method are: \n")
            cat(availclasses)
            cat(" \n")
            stop("object(s) of a class not handled")
          }
          res <- matrix(NA, nrow = length(clobject), ncol = 3)

          RVAL <- sapply(object[clobject],
                         FUN = function(objectj) {
                           if (inherits(objectj, "glm")) {
                             nj <- length(objectj$residuals)
                             rankj <- objectj$df.null - objectj$df.residual + 1
                           }
                           else {
                             nj <- objectj$nobs
                             rankj <- objectj$rank
                           }

                           R2j <- get.R2VAL (object = objectj, n = nj, rank = rankj, adjust.size = adjust.size)

                           return(c(R2 = R2j, N = nj, df = rankj))
                         })
          res[clobject,] <- if (sum(clobject) > 0) t(RVAL) else RVAL
          RVAL <- data.frame(R2 = res[,1], N = res[,2], df = res[,3])
          colnames(RVAL) <- c(method, 'N', 'df')

          if (!is.null(nm))
            row.names(RVAL) <- nm

          return(RVAL)
        }

        stop(paste0(" object of a class not handled by 'rsquared.default'. Available classes are: ",
                    paste(availclasses, collapse = ', '), '.'))
      }
    }
    else {
      objects <- list(object, ...)
      out <- rsquared.default(objects, method = method, adjust.size = adjust.size)
      rownames(out) <- make.unique(as.character(Call[-1L])[1:NROW(out)])
      return(out)
    }
  }

get.KLR2VAL <- function (object, n, rank, adjust.size) {

  R2 <- 1 - (object$deviance/object$null.deviance) * (
    ifelse (adjust.size, (n - 1) / (n - rank), 1) )

  return(R2)
}

#' @importFrom stats cor
#' @importFrom stats cov.wt
#'
get.COR2VAL <- function (object, n, rank, adjust.size,
                         method = c("pearson", "kendall", "spearman")) {
  y <- object$y %||% model.frame(object)$y

  if (is.null(y))
    stop("no 'y' component in 'object'")

  weights <- object$weights
  if (is.null(weights))
    weights <- 1

  yhat <- fitted(object) * weights # useless to multiply by 'object$weights'
  sweights <- object$sample.weights

  if (length(sweights) > 1) {
    sweights <- n * sweights / sum(sweights, na.rm = TRUE)
  }
  else {
    sweights <- 1
  }

  method <- match.arg(method)

  goPearson <- expression({
    if (all(sweights == 1) & all(weights == 1)) {
      # From definition on https://en.wikipedia.org/wiki/Point-biserial_correlation_coefficient
      y1s <- y == y[1]
      Delta.yhat <- mean(yhat[y1s], na.rm = TRUE) -
        mean(yhat[!y1s], na.rm = TRUE)
      nnm <- sum(!is.na(yhat))
      var.yhat <- var(yhat, na.rm = TRUE) * (nnm - 1) / nnm
      pzeroes <- mean(y1s, na.rm = TRUE)
      R2 <- (Delta.yhat^2) * pzeroes * (1 - pzeroes) / var.yhat
    }
    else if (any(sweights != 1)) {
      R2 <- stats::cov.wt (cbind(yhat, y), w = sweights,
                           method = "unbiased", cor = TRUE,
                           na.rm = TRUE)$cor[1,2]^2
    }
    else {
      R2 <- stats::cor (x=yhat, y = y, method = method,
                        use = "complete.obs")^2
    }
  })
  switch(method,
         pearson = {
           eval(goPearson)
         },
         spearman = {
           Rank <- function(u) {
             if (length(u) == 0L)
               u
             else if (is.matrix(u)) {
               if (nrow(u) > 1L)
                 apply(u, 2L, rank, na.last = "keep")
               else row(u)
             }
             else rank(u, na.last = "keep")
           }
           y <- Rank(y)
           yhat <- Rank(yhat)
           eval(goPearson)
         },
         kendall = {
           if (all(sweights == 1)) {
             R2 <- stats::cor (x=yhat, y = y, method = method,
                               use = "complete.obs")^2
           }
           else {
             warning("Kendall tau not implemented for sample-weighted binary response")
             R2 <- NA
           }
         })

  if (adjust.size) {
    R2 <- 1 - R2 * (n - 1) / (n - rank)
  }

  return(R2)
}

#' @rdname rsquared
#' @exportS3Method msbreg::rsquared
rsquared.glm <- function(object, ...,
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
