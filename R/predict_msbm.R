#' Prediction Method for \code{msbm} Objects
#'
#' Extract various types of predictions from Multistage Binomial (MSB)
#' regression models, on the scale of either the response variable,
#' the success probability (in \code{[0, 1]}), or the linear predictors
#' in the MSB model.
#'
#' @param object fitted model object of class \code{"msbm"}, as returned by
#' \link[msbreg]{msbreg}.
#'
#' @param newdata optional, data frame in which to look for variables
#' with which to predict. If omitted, the original observations are used.
#'
#' @param input.me,stage.me optional named list of one sided
#' \link[stats]{formula}(s) (\code{input.me}) specifying measurement
#' errors in some numerical predictors, and a one sided
#' \link[stats]{formula} (\code{stage.me}) specifying measurement
#' errors in some stages in the formula of the fitted model
#' \code{object}. See \link[msbreg]{msbm.frame} for details.
#'
#' The defaults depend on the presence of argument \code{newdata}.
#'
#' \itemize{
#' \item \code{newdata} is missing or \code{NULL}: the default
#' \code{input.me} and \code{stage.me} are the corresponding components of
#' \code{object$call} (since the same data used for model fitting is
#' reconsidered, the default assumes that there are measurement errors,
#' if any).
#'
#' \item \code{newdata} is supplied and not \code{NULL}: the defaults
#' are \code{input.me = list()} and \code{stage.me = NULL} which correspond
#' to no measurement error in any predictor or stage.
#' }
#'
#' @param type character indicating the desired type of prediction. One of:
#'
#' \itemize{
#' \item \code{"response"}: vector of fitted means of binomial responses;
#'
#' \item \code{"prob"}: vector of fitted success probabilities;
#'
#' \item \code{"link"}: matrix of linear predictors for model components
#' (each model stage, and the \eqn{\alpha} and \eqn{\lambda} components);
#'
#' \item \code{"quantile"}: vector of fitted quantile(s) of the binomial
#' response distribution;
#' }
#'
#' The default is \code{type = "response"}.
#' The option \code{"prob"} is simply the fitted response divided by the
#' number of trials (weights). Thus, \code{type = "response"} and
#' \code{type = "prob"} are equivalent if the response is binary
#' (all weights are one).
#'
#' @param at numeric vector indicating the level(s) at which quantiles
#' should be predicted (only if \code{type = "quantile"}),
#' defaulting to the median (\code{at = 0.5)}.
#'
#' @param se.fit logical switch indicating if standard errors are desired
#' when  \code{type = "response"},  \code{type = "prob"},
#' or \code{type = "link"}.
#' Note that when \code{type = quantile}, no standard error is available.
#'
#' @param fisher.matrix logical, should the Fisher (expected) information
#' matrix be used to obtain standard errors? Defaults to \code{TRUE}.
#' Only used when \code{se.fit = TRUE}.
#'
#' @param na.action function determining what should be done with missing
#' values in \code{newdata}. The default is to predict \code{NA}.
#'
#' @param ... further arguments passed to or from other methods.
#' None is currently used.
#'
#' @details
#' Variables are first looked for in \code{newdata} and then searched for
#' in the usual way (which will include the environment of the formula
#' used in the fit). A warning will be given if the variables found are
#' not of the same length as those in \code{newdata} if it was supplied.
#'
#' If the argument \code{newdata} is omitted, predictions are based
#' on the data used for the fit. In that case, any measurement error
#' in the data used for the fit is reconsidered for the predictions
#' (note however that have no impact on link scale predictions, given
#' that measurement error distributions are symmetric and have zero mean).
#' To reuse the data used for model fitting but ignore measurement
#' error in (some of) the predictors, either set \code{input.me = list()}
#' and \code{stage.me = NULL}, or leave \code{input.me} and \code{stage.me}
#' unspecified while passing the dataset used for model fitting as
#' \code{newdata}.
#'
#' When \code{se.fit = TRUE}, standard error computation is based on
#' the delta-method, starting from standard error for model estimated
#' parameters. Since quantile extraction is not a continuous operation
#' under the binomial distribution, no standard error is available if
#' \code{type = "quantile"}, and \code{NA}s are returned as standard
#' errors in that case.
#'
#' @return If \code{se.fit = FALSE}, a vector or matrix of predictions.
#' For \code{type = 'link'}, this is a matrix with a column per model
#' component (each model stage, and the \eqn{\alpha} and \eqn{\lambda}
#' components).
#'
#' If \code{se.fit = TRUE}, a list with two components:
#'
#' \itemize{
#' \item \code{fit}, predictions (as for \code{se.fit = FALSE});
#' \item \code{se.fit}, estimated standard errors.
#' }
#'
#' @importFrom stats qbinom
#'
#' @examples
#' data(test1data)
#' attr(test1data$y, "formula")
#'
#' ##* Two-stage logistic model fit with unknown maximum success probability
#' # A multiplicative intercept included
#' msbfit <- msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                    lambda.formula = ~ 1,
#'                    weights = Total,
#'                    data = test1data)
#'
#' summary (msbfit)
#'
#' #* Expected binomial responses for a few observations
#' predict(msbfit, test1data[1:10,])
#'
#' #* Expected success probability and standard errors
#' predict(msbfit, test1data[1:10,], type = "prob", se.fit = TRUE)
#'
#'
#' @exportS3Method stats::predict
predict.msbm <- function (object,
                          newdata = NULL, input.me, stage.me,
                          type = c("response", "prob", "link", "quantile"),
                          at = 0.5, se.fit = FALSE, fisher.matrix = TRUE,
                          na.action = na.pass, ...) {
  type <- match.arg(type)
  na.act <- object$na.action
  object$na.action <- NULL

  se.fit <- as.logical(se.fit)[1]
  fisher.matrix <- as.logical(fisher.matrix)[1]
  if (is.null(newdata) & (missing(input.me) &
                          missing(stage.me))) {
    switch(type,
           prob = {
             pred <- object$fitted.values
           },
           response = {
             pred <- object$fitted.values * object$weights
           },
           link = {
             pred <- object$linear.predictors
           },
           quantile = {
             stopifnot(is.numeric(at))
             stopifnot(all(at >= 0))
             stopifnot(all(at <= 1))
             pred <- sapply(at, FUN = function (pv) {
               stats::qbinom (p = pv, size = object$weights, prob = object$fitted.values,
                              lower.tail = TRUE, log.p = FALSE)
             })
             if (length(at) == 1)
               pred <- c(pred)
           })


    if (!is.null(na.act))
      pred <- napredict(na.act, pred)

    if (se.fit & !identical(type, "quantile")) {
      smry <- summary(object,
                      fisher.matrix = fisher.matrix,
                      se.values = TRUE, frame = FALSE)

      sepred <- switch(type,
                       prob = smry$se.mu,
                       response = smry$se.mu * object$weights,
                       NA)

      if (identical(type, "link")) {
        frame <- object$frame %||% model.frame(object)
        eval(get.msbm.dims())
        endpj <- cumsum(pj + intercepts)
        startpj <- c(0, endpj[-q]) + 1
        nnn <- NROW(frame$input.matrix)

        ffun <- function(j) {
          dict <- frame$stage.dictionary[[j]]

          if (attr(dict, 'only.offset')) {
            return(rep(0, length.out = nnn))
          }

          covj <- smry$cov.scaled[startpj[j]:endpj[j], startpj[j]:endpj[j]]
          if (pj[j] > 0) {
            X <- frame$input.matrix[, dict, drop = FALSE]

            if (attr(dict, 'intercept')) {
              X <- cbind(1, X)
            }
            sej <- sqrt(rowSums ((X %*% covj) * X))
          }
          else {
            sej <- rep(sqrt(covj[1]), length.out = nnn)
          }

          return(sej)
        }

        if (q == 1) {
          sepred <- cbind(ffun(1))
        }
        else {
          sepred <- sapply(1:q, FUN = ffun)
        }
        colnames(sepred) <- names(frame$stage.dictionary)

        if (d > 0) {
          cov_alpha <- smry$cov.scaled[(endpj[q] + 1):(endpj[q] + d), (endpj[q] + 1):(endpj[q] + d)]

          if (attr(frame$alpha.dictionary, 'intercept')) {
            if (d == 1) {
              se.eta.alpha <- sqrt(cov_alpha)
            }
            else {
              Xalpha <- cbind(1, frame$input.matrix[, frame$alpha.dictionary, drop = FALSE])
              se.eta.alpha <- sqrt(rowSums((Xalpha %*% cov_alpha) * Xalpha))
            }
          }
          else {
            Xalpha <- frame$input.matrix[, frame$alpha.dictionary, drop = FALSE]
            se.eta.alpha <- sqrt(rowSums((Xalpha %*% cov_alpha) * Xalpha))
          }
          sepred <- cbind(sepred, alpha = se.eta.alpha)
        }
        else {
          sepred <- cbind(sepred, alpha = 0)
        }

        if (r > 0) {
          cov_lambda <- smry$cov.scaled[(endpj[q] + d + 1):(endpj[q] + d + r), (endpj[q] + d + 1):(endpj[q] + d + r)]

          if (attr(frame$lambda.dictionary, 'intercept')) {
            if (r == 1) {
              se.eta.lambda <- sqrt(cov_lambda)
            }
            else {
              Xlambda <- cbind(1, frame$input.matrix[, frame$lambda.dictionary, drop = FALSE])
              se.eta.lambda <- sqrt(rowSums((Xlambda %*% cov_lambda) * Xlambda))
            }
          }
          else {
            Xlambda <- frame$input.matrix[, frame$lambda.dictionary, drop = FALSE]
            se.eta.lambda <- sqrt(rowSums((Xlambda %*% cov_lambda) * Xlambda))
          }
          sepred <- cbind(sepred, lambda = se.eta.lambda)
        }
        else {
          sepred <- cbind(sepred, lambda = 0)
        }

      }

      if (!is.null(na.act))
        sepred <- napredict(na.act, sepred)

      pred <- list(fit = pred, se.fit = sepred)
    }
    else if (se.fit & identical(type, "quantile")) {
      pred <- list(fit = pred, se.fit = NA)
    }
  }
  else {
    if (is.null(newdata)) {
      newdata <- object$data
    }

    {
      mf <- object$call
      m <- match(x = c("formula", "input.me", "stage.me",
                       "alpha.formula", "lambda.formula",
                       "data", "weights", "sample.weights", "subset",
                       "na.action", "drop.unused.levels", "frames"),
                 table = names(mf),
                 nomatch = 0L)
      mf <- mf[c(1L, m)]
      mf[[1L]] <- quote(msbreg::msbm.frame)
      mfcall <- mf
      mf$drop.unused.levels <- TRUE
      mf$verbose <- FALSE
      mf$data <- newdata
      mf$na.action <- na.action
      if (!missing(input.me)) {
        mf$input.me <- input.me
      }
      if (!missing(stage.me)) {
        mf$stage.me <- stage.me
      }
      frame <- catch.conditions({
        eval(mf, parent.frame())
      })$value

      if (any(class(frame) %in% c("simpleError", "error",
                                  "condition", "try-error"))) {
        frame <- eval(mf, environment(formula))
      }
    }

    mutated <- mutate.params (theta = object$coefficients,
                              frame = frame,
                              link = object$link,
                              information = se.fit,
                              observed = se.fit & !fisher.matrix)

    switch(type,
           response = {
             pred <- mutated$mu * frame$weights
           },
           link = {
             pred <- mutated$etas
           },
           prob = {
             pred <- mutated$mu
           },
           quantile = {
             stopifnot(is.numeric(at))
             stopifnot(all(at >= 0))
             stopifnot(all(at <= 1))
             pred <- sapply(at, FUN = function (pv) {
               stats::qbinom (p = pv, size = frame$weights, prob = mutated$mu,
                              lower.tail = TRUE, log.p = FALSE)
             })
             if (length(at) == 1)
               pred <- c(pred)
           })
    if (!is.null(na.act))
      pred <- napredict(na.act, pred)

    if (se.fit & !identical(type, "quantile")) {
      cov.scaled <- catch.conditions({
        if (fisher.matrix)
          solve(mutated$EFinfo)
        else
          solve(mutated$OFinfo)
      })$value

      if (any(class(cov.scaled) %in% c("simpleError", "error",
                                          "condition", "try-error"))) {
        sepred <- NA
      }
      else {

        switch(type,
               prob = {
                 dmu <- mutated$mu.theta
                 sepred <- rowSums ((dmu %*% cov.scaled) * dmu)
               },
               response = {
                 dmu <- mutated$mu.theta
                 sepred <- rowSums ((dmu %*% cov.scaled) * dmu) * frame$weights
               },
               {
                 sepred <- NA
               })

        if (identical(type, "link")) {
          eval(get.msbm.dims())
          endpj <- cumsum(pj + intercepts)
          startpj <- c(0, endpj[-q]) + 1
          nnn <- NROW(frame$input.matrix)

          ffun <- function(j) {
            dict <- frame$stage.dictionary[[j]]

            if (attr(dict, 'only.offset')) {
              return(rep(0, length.out = nnn))
            }

            covj <- cov.scaled[startpj[j]:endpj[j], startpj[j]:endpj[j]]
            if (pj[j] > 0) {
              X <- frame$input.matrix[, dict, drop = FALSE]

              if (attr(dict, 'intercept')) {
                X <- cbind(1, X)
              }
              sej <- sqrt(rowSums ((X %*% covj) * X))
            }
            else {
              sej <- rep(sqrt(covj[1]), length.out = nnn)
            }

            return(sej)
          }

          if (q == 1) {
            sepred <- cbind(ffun(1))
          }
          else {
            sepred <- sapply(1:q, FUN = ffun)
          }
          colnames(sepred) <- names(frame$stage.dictionary)

          if (d > 0) {
            cov_alpha <- cov.scaled[(endpj[q] + 1):(endpj[q] + d), (endpj[q] + 1):(endpj[q] + d)]

            if (attr(frame$alpha.dictionary, 'intercept')) {
              if (d == 1) {
                se.eta.alpha <- sqrt(cov_alpha)
              }
              else {
                Xalpha <- cbind(1, frame$input.matrix[, frame$alpha.dictionary, drop = FALSE])
                se.eta.alpha <- sqrt(rowSums((Xalpha %*% cov_alpha) * Xalpha))
              }
            }
            else {
              Xalpha <- frame$input.matrix[, frame$alpha.dictionary, drop = FALSE]
              se.eta.alpha <- sqrt(rowSums((Xalpha %*% cov_alpha) * Xalpha))
            }
            sepred <- cbind(sepred, alpha = se.eta.alpha)
          }
          else {
            sepred <- cbind(sepred, alpha = 0)
          }

          if (r > 0) {
            cov_lambda <- cov.scaled[(endpj[q] + d + 1):(endpj[q] + d + r), (endpj[q] + d + 1):(endpj[q] + d + r)]

            if (attr(frame$lambda.dictionary, 'intercept')) {
              if (r == 1) {
                se.eta.lambda <- sqrt(cov_lambda)
              }
              else {
                Xlambda <- cbind(1, frame$input.matrix[, frame$lambda.dictionary, drop = FALSE])
                se.eta.lambda <- sqrt(rowSums((Xlambda %*% cov_lambda) * Xlambda))
              }
            }
            else {
              Xlambda <- frame$input.matrix[, frame$lambda.dictionary, drop = FALSE]
              se.eta.lambda <- sqrt(rowSums((Xlambda %*% cov_lambda) * Xlambda))
            }
            sepred <- cbind(sepred, lambda = se.eta.lambda)
          }
          else {
            sepred <- cbind(sepred, lambda = 0)
          }
        }
      }

      if (!is.null(na.act))
        sepred <- napredict(na.act, sepred)

      pred <- list(fit = pred, se.fit = sepred)
    }
    else if (se.fit & identical(type, "quantile")) {
      pred <- list(fit = pred, se.fit = NA)
    }
  }

  return(pred)
}
