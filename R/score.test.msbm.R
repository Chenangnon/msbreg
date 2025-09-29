# score.test
#' Score Test for Fitted Model Objects
#'
#' Performs score tests comparing a fitted model to alternatives.
#'
#' @param object an object of a class with a \code{score.test} method.
#' Typically, a fit, for instance \code{"lm"}, \code{"glm"} and \code{"msbm"}
#' class objects.
#'
#' @param ... further arguments passed to or from other methods.
#' This may include one or more alternative models of the same class as \code{object};
#' and arguments such as \code{which} indicating which parameters to test,
#' \code{type} for using the expected (Fisher) or observed information matrix,
#' \code{trace} for verbose, and \code{ncores} for the number of cores for
#' parallel computing.
#'
#' @param which integer vector, which score types statistic of a
#' \code{score.test} object should be extracted?
#' The default selects score statistics based on the Fisher information matrix.
#' Do not change it unless you are sure about what you are doing.
#'
# @param which character or integer indicating which parameters are targeted.
# Defaults to \code{which = NULL}, corresponding to the asymptote parameter
# \eqn{\lambda} of the Multistage Binomial (MSB) model.
# No alternative is currently available.
#
# @param null.values numeric vector of null values of target parameters.
#'
#' @param scope defines the range of alternative models to be examined.
#' This is used only when \code{...} does not include any model which nests
#' \code{object}.
#' This should be either a single formula, or a list containing components
#' \code{upper} and \code{lower}, both formulae.
#' See the details for how to specify the formulae and how they are used.
#'
#' @param steps the maximum number of null models to be considered.
#' The default is 1000 (essentially as many as required).
#' It is typically used to stop the process early when \code{...} is empty
#' and \code{scope} is a big model.
#'
#' @param type character, type of score statistic to compute. The default
#' (\code{type = "expected"}) uses the Fisher information matrix.
#' The alternative (\code{type = "observed"}) uses the observed information
#' matrix.
#'
#' @param penalty logical, should the test account for the penalty when a
#' penalized fitting criterion was used for object?
#'
#' @param identity.link.lambda should the identity link function be considered
#' for the asymptote parameter \eqn{\lambda} in a MSB model?
#' The default (\code{identity.link.lambda = FALSE}) uses a \code{log} link.
#'
#' Note that the link used by fitter function is generally a binomial link
#' function which represents the \eqn{\lambda = 1} as \eqn{\infty} and
#' results into an identically zero score component for \eqn{\lambda} and
#' thereby a singular Fisher information matrix.
#' A different parametrization is thus required/better for \eqn{\lambda} when
#' performing any of the statistical trinity (score, likelihood ratio, and
#' Wald tests).
#'
#' @param side \code{NULL} \code{NA}, \code{-1}, or \code{1} indicating if
#' one-sided derivatives should be used when numerical derivative are required
#' (only relevant when \code{type = "observed"} and \code{penalty = TRUE}).
#'
#' @param trace integer indicating if details should be printed during
#' computations.
#'
#' @param ncores integer, number of threads for parallel processing.
#'
#' @details
#' This is a generic function to perform the Rao score test or modified versions
#' on a model fit object.
#'
#' The method for an \code{object} of class "\link{lm}" or "\link{glm}" proceeds
#' as follows.
#' If a unique model is provided, the default method tests for
#'
#'
#'
#' When
#'
#' The function \code{scoretable} generates a summary of one or many
#' \code{score.test} results.
#'
#' @return An object of class ‘score.test’ that inherits from the class ‘htest’,
#' a list with components:
#'
#'
#'
#' The \code{print} method for \code{"score.test"} class objects calls
#' \code{scoretable}.
#'
#' The routine \code{extractscores} returns a vector for a \code{score.test}
#' input \code{object} and a matrix when \code{object} is a \code{score.list}.
#'
#' @aliases score.test.msbm
#'
#' @seealso See \link[msbreg]{profile.msbm} for computing log-likelihood profiles
#' for MSB model fits.
#'
#'
# @export score.test.msbm
#' @exportS3Method msbreg::score.test
#'
score.test.msbm <- function (object, ...,  scope, steps = 1000,
                             type = c("expected", "observed"), penalty = TRUE,
                             identity.link.lambda = FALSE, side = NULL,
                             trace = object$control$trace, ncores = NULL) {
  out <- score.test.msbm.intern (object, ...,
                                 type = type, penalty = penalty, side = side,
                                 identity.link.lambda = identity.link.lambda,
                                 dname = deparse1(substitute(object)),
                                 trace = trace, ncores = ncores)
  out$call <- match.call()

  out <- structure(list(out),
                   class = "score.list")

  out
}

score_msbm <- function (object, ...,
                        scope,
                        lambda.scope,
                        alpha.scope,
                        kappa.scope,
                        steps = 1000,
                        type = c("expected", "observed"), penalty = TRUE,
                        identity.link.lambda = FALSE,
                        dname = deparse1(substitute(object)),
                        trace = FALSE, ncores = NULL) {
  alternatives <- list(...)
  if (length(alternatives)) {
    keep1 <- sapply(alternatives, FUN = is, class2 = "msbm.frame")
    keep2 <- sapply(alternatives, FUN = is, class2 = "msbm")
    keep <- keep1 | keep2
    if (!any(keep)) {
      alternatives <- list()
      if (trace) {
        warning("no 'msbm.frame' or 'msbm' ob object found in the '...' arguments: ignoring '...'")
      }
    }
    else {
      alternatives <- alternatives[keep]
      alt.frames <- keep1[keep]
    }
  }
  mdots <- length(alternatives) > 0

  if (!mdots & !all(c(missing(scope), missing(lambda.scope),
                                     missing(alpha.scope), missing(kappa.scope)))) {
    #* formula (stages)
    if (missing(scope)) {
      formula.lower <- formula.upper <- NULL
    }
    else {
      if (is.list(scope)) {
        formula.lower <- scope$lower
        formula.upper <- scope$upper
      }
      else {
        formula.lower <- object$call$formula
        formula.upper <- scope
      }
    }

    #* lambda.formula (asymptote)
    if (missing(lambda.scope)) {
      lambda.formula.lower <- lambda.formula.upper <- NULL
    }
    else {
      if (is.list(lambda.scope)) {
        lambda.formula.lower <- lambda.scope$lower
        lambda.formula.upper <- lambda.scope$upper
      }
      else {
        lambda.formula.lower <- object$call$lambda.formula %||% ~ 1
        lambda.formula.upper <- lambda.scope
      }
    }

    #* alpha.formula (minimum succes probability)
    if (missing(alpha.scope)) {
      alpha.formula.lower <- alpha.formula.upper <- NULL
    }
    else {
      if (is.list(alpha.scope)) {
        alpha.formula.lower <- alpha.scope$lower
        alpha.formula.upper <- alpha.scope$upper
      }
      else {
        alpha.formula.lower <- object$call$alpha.formula %||% ~ 0
        alpha.formula.upper <- alpha.scope
      }
    }

    #* kappa.formula (outcome correlation parameter)
    if (missing(kappa.scope)) {
      kappa.formula.lower <- kappa.formula.upper <- NULL
    }
    else {
      if (is.list(kappa.scope)) {
        kappa.formula.lower <- kappa.scope$lower
        kappa.formula.upper <- kappa.scope$upper
      }
      else {
        kappa.formula.lower <- object$call$kappa.formula %||% ~ 0
        kappa.formula.upper <- kappa.scope
      }
    }

    #* ADD code to generate alternative model frames based on the inpute scopes

  }

  # Case no alternative model is supplied
  if (!length(alternatives)) {
    out <- score.test.msbm.intern(object,
                                  type = type, penalty = penalty, dname = dname,
                                  trace = trace, ncores = ncores,
                                  identity.link.lambda = identity.link.lambda)


    out <- structure(list(out),
                     class = "score.list")

    return(out)
  }

  # Case the user provided alternative models
  if (mdots) {

  }

  # Case the user did not provide alternative models


}

#' @exportS3Method base::print
print.score.list <- function(x, ...) {
  if(length(x) == 1) {
    # Just print the 'htest' class result
    print(x[[1]])
  }

  # Use score table to summarize results from many score tests

}

# Basic score test for ONE fitted msbm model
score.test.msbm.intern <- function (object, ...,
                                    type = c("expected", "observed"), penalty = TRUE,
                                    side = NULL, dname = deparse1(substitute(object)),
                                    trace = object$control$trace, ncores = NULL,
                                    identity.link.lambda = FALSE)  {

  setNumThreads(ncores)
  type < match.arg(type, choices = c("expected", "observed"),
                   several.ok = TRUE)
  type <- type[1]
  penalty <- as.logical(penalty)[1] + 1
  model.call <- object$call

  #* Model frame
  frame <- model.frame(object)
  control <- object$control
  link <- object$link
  direct.lambda0 <- identity.link.lambda
  objective <- logLikvalue <- rep(NA, 2)

  if (length(frame$lambda.dictionary) == 0 & attr(frame$lambda.dictionary, "intercept") == 0) {
    #* Consider model with lambda = 1 as NULL model (if object$unit.asymptote is not null, use it for initialization)
    #* Compare with fitted model in object

    #* Specify the identity link for lambda, if requested
    direct.lambda <- direct.lambda0

    #* Null model
    NULLmodel <- object
    H1theta <- NA

    #* Full model call
    object <- object$call
    object$lambda.formula <- ~ 1
    object$start <- NULLmodel$coefficients
    # msbm.fit automatically adjust the start argument
    # It includes `(Intercept).lambda` = link$linkfun(.99) for the lambda parameter

    #* Update the msbm.frame object
    object$method <- "model.frame"
    frame <- try ({eval(object)}, silent = TRUE) # model.frame(object)
    if (is(frame, 'try-error')) {
      frame <- try ({eval(object, envir = environment(NULLmodel$formula))}, silent = TRUE)
    }
    if (is(frame, 'try-error')) {
      nparent <- sum(sapply(1:10, FUN = sys.parent) > 0)

      if (nparent > 0) {
        for (iparent in seq_len(nparent)) {
          if (is(frame, 'try-error'))
            frame <- try ({eval(object, envir = parent.frame(iparent))}, silent = TRUE)
          else
            break
        }
        if (is(frame, 'try-error'))
          stop(paste0("alternative model frame evaluation failled: ", frame))
      }
      else {
        stop(paste0("alternative model frame evaluation failled: ", frame))
      }
    }

    #* Model AIC
    AIC <- c(AIC(NULLmodel), NA)

  }
  else {
    #* Consider model with lambda = 1 as NULL model (if object$unit.asymptote is not null, use it for initialization)
    #* Compare with fitted model in object
    #*
    #* Problem: the information is singular!
    #*
    #* We circumvent the singularity issue using a log link function (with restriction) for lambda
    #* Consider the parametrization: lambda = exp(delta[1] + delta[-1] x V) with the general (sample dependent) restrictions delta[1] <= 0 delta[-1] x V <= 0.
    #* The test is now for H0: delta = 0.
    #*

    #* Specify the identity or log link for lambda
    if (length(frame$lambda.dictionary) > 0 | attr(frame$lambda.dictionary, "intercept") == 0) {
      direct.lambda <- FALSE
    }
    else {
      direct.lambda <- direct.lambda0
    }

    #* Null model call
    NULLcall <- object$call
    NULLcall$lambda.formula <- ~ 0

    # Initialization
    if (length(frame$lambda.dictionary) == 0 &
        attr(frame$lambda.dictionary, "intercept") == 1 &
        !is.null(object$unit.asymptote)) {
      NULLcall$start <- object$unit.asymptote$par
    }
    else {
      NULLcall$start <- object$coefficients[1:(frame$dims$p + frame$dims$d)]
    }

    # Fit the null model
    NULLmodel <- try ({eval(NULLcall)}, silent = TRUE)
    if (is(NULLmodel, 'try-error')) {
      NULLmodel <- try ({eval(NULLcall, envir = environment(object$formula))}, silent = TRUE)
    }
    if (is(NULLmodel, 'try-error')) {
      nparent <- sum(sapply(1:10, FUN = sys.parent) > 0)

      if (nparent > 0) {
        for (iparent in seq_len(nparent)) {
          if (is(NULLmodel, 'try-error'))
            NULLmodel <- try ({eval(NULLcall, envir = parent.frame(iparent))}, silent = TRUE)
          else
            break
        }
        if (is(NULLmodel, 'try-error'))
          stop(paste0("null model evaluation failled: ", NULLmodel))
      }
      else {
        stop(paste0("null model evaluation failled: ", NULLmodel))
      }
    }
    H1theta <- object$coefficients

    #* Model AIC
    AIC <- c(AIC(NULLmodel), AIC(object))
  }
  objective[1] <- NULLmodel$criterion
  logLikvalue[1] <- logLik(NULLmodel)[1]
  exp.lambda <- !direct.lambda

  #* Null parameter vector
  theta <- c(NULLmodel$coefficients, rep(NA, frame$dims$r))
  if (exp.lambda) {
    theta[-(1:(frame$dims$p + frame$dims$d))] <- rep(0, frame$dims$r)
  }
  else {
    theta[-(1:(frame$dims$p + frame$dims$d))] <- c(Inf, if (frame$dims$r > 1) rep(0, frame$dims$r - 1))
  }
  if (frame$dims$r == 1 & attr(frame$lambda.dictionary, "intercept")) {
    names(theta) <- c(names(NULLmodel$coefficients), "(Intercept).lambda")
  }

  #* Score vector and expected information matrix for the log-likelihood criterion
  mutatedtheta <- mutate.params (theta, frame = frame,
                                 link = link,
                                 information = TRUE,
                                 observed = identical(type, "observed"),
                                 jerk = !identical(control$criterion, "ML"),
                                 direct.lambda = direct.lambda,
                                 exp.lambda = exp.lambda)

  #* Save likelihood
  logLikvalue[2] <- mutatedtheta$loglike

  # Compute the score statistics
  covscores_e <- try({mp.inverse(mutatedtheta$EFinfo)}, silent = TRUE)
  if (is(covscores_e, 'try-error')) {
    S_e <- NA
    objective[2] <- NA
  }
  else {
    S_e <- sum(mutatedtheta$score * c(covscores_e %*% mutatedtheta$score))
    objective[2] <- - mutatedtheta$loglike -
      if(identical(control$criterion, "ML")) 0 else 0.5 * sum(log(attr(covscores_e, 'singular.values')))
  }

  if (is.null(mutatedtheta$OFinfo)) {
    S_o <- NA
  }
  else {
    covscores_o <- try({mp.inverse(mutatedtheta$OFinfo)}, silent = TRUE)
    if (is(covscores_o, 'try-error')) {
      S_o <- NA
    }
    else {
      S_o <- sum(mutatedtheta$score * c(covscores_o  %*% mutatedtheta$score))
    }
  }

  statisticlist <- list(logLik = c(S_e, S_o))
  names(statisticlist$logLik) <- c("expected", "observed")

  #* Score statistic based on fitting criterion if penalized-likelihood
  if (!identical(control$criterion, "ML")) {
    if (is(covscores_e, 'try-error')) {
      statisticlist$criterion <- c(NA, NA)
      names(statisticlist$criterion) <- c("expected", "observed")
    }
    else {
      penscore <- mutatedtheta$score + 0.5 * c(t(mutatedtheta$Einfo.theta) %*% c(covscores_e))
      S_e <- sum(penscore * c(covscores_e %*% penscore))

      if (identical(type, "observed")) {
        #* (penalized-)deviance function of the alternative model
        eval(toget_half_devfun())

        #* Score vector and observed information matrix for the fitting criterion
        if (!is.null(side)) {
          stopifnot(all(side %in% c(-1, 1, NA)))
          if (length(side) != frame$dims$npars) {
            side <- numeric(frame$dims$npars) + NA
            side[-(1:(frame$dims$p + frame$dims$d))] <- -1
          }
        }

        thetaNum <- theta
        if (!exp.lambda) {
          thetaNum[frame$dims$p + frame$dims$d+1] <- link$linkfun(1-control$epsilon)
        }

        if (!is.null(side)) {
          scfun <- function(x) c(numDerivjacobian (half_devfun, x = x, side = side))
          objhessian <- numDerivjacobian (scfun, x = thetaNum, side = side)
        }
        else {
          objhessian <- numDerivhessian (half_devfun, x = thetaNum)
        }

        covscores_o <- try({mp.inverse(objhessian)}, silent = TRUE)
        if (is(covscores_o, 'try-error')) {
          S_o <- NA
        }
        else
          S_o <- sum(penscore * c(covscores_o %*% penscore))

        statisticlist$criterion <- c(S_e, S_o)
        names(statisticlist$criterion) <- c("expected", "observed")
      }
      else {
        statisticlist$criterion <- c(S_e, NA)
        names(statisticlist$criterion) <- c("expected", "observed")
      }
    }
  }
  else {
    statisticlist$criterion <- statisticlist$logLik
  }

  np <- c(frame$dims$npars - frame$dims$r, frame$dims$npars)
  parameter <- c(df = frame$dims$r)
  AICc <- ifelse((frame$dims$nobs - np - 1) > 0, AIC + 2 * np * (np + 1)/(frame$dims$nobs - np - 1), NA)
  names(AIC) <- names(AICc) <- c("H0", "H1")

  ### ADD score and FIM to output

  out <- list(statistic = c(`Xbar-squared` = as.numeric(statisticlist[[penalty]][type])),
              parameter = parameter,
              p.value = 0.5 * pchisq(statisticlist[[penalty]][type], df = parameter, lower.tail = FALSE),
              method = "Score test for H0: lambda = 1 under MSB models",
              data.name = dname,
              H0 = "H0: lambda = 1",
              H0theta = theta,
              H0objective = c(H0.logLik = as.numeric(logLikvalue[1]),
                              H0.criterion = as.numeric(objective[1])),
              H1theta = H1theta,
              H1objective = c(H1.logLik = as.numeric(logLikvalue[2]),
                              H1.criterion = as.numeric(objective[2])),
              model.call = model.call,
              scores = statisticlist,
              AIC = AIC, AICc = AICc)

  class(out) <- c("score.test", "htest")

  return(out)
}
