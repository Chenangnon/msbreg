#
#' Binomial \code{glm} to \code{msbm}
#'
#' Function to convert a binomial \code{"glm"} class object into a
#' \code{"msbm"} class object.
#'
#' @param object an object inheriting from the class \code{"glm"},
#' and fitted with the \link[stats]{binomial} family (that is,
#' \code{object$family$family = "binomial"}).
#'
#' @param refit logical, should a new fit be performed to replace the
#' \code{coefficient} component (and related quantities) of the
#' returned value?
#'
#' @param criterion fitting criterion, only used if \code{refit = TRUE}
#' (see \link[msbreg]{msbm.control} and \link[msbreg]{msbreg}).
#'
#' @param method,start,control,y,frame arguments for \link[msbreg]{msbreg},
#' only used if \code{refit = TRUE}.
#'
#' @param ... further arguments to pass to or from other methods.
#'
#' Currently, \code{alpha.formula}, \code{lambda.formula}, \code{me.input}
#' and \code{stage.me} are passed to \link[msbreg]{msbm.frame} when
#' \code{refit = TRUE} or \code{method = "model.frame"}.
#'
#' @details
#' This function aims to provide a simple routine to allow application
#' of methods developed for the \code{"msbm"} class to \code{"glm"} class
#' objects.
#'
#' @return an "\link[msbreg]{msbm.frame}" class object if
#' \code{method = "model.frame"}. Otherwise, an object inheriting
#' from class \code{"glm2msbm"} which is essentially an \code{"msbm"}
#' class object, but some components of the returned object miss
#' some attributes present in regular \code{"msbm"} objects.
#'
#' @export glm.to.msbm
#' @importFrom stats residuals
#'
#' @examples
#' data("infert", package = "datasets")
#'
#' ## Logistic regression fit
#' GLMfit <- glm (case ~ spontaneous,
#'                data = infert,
#'                family = binomial())
#'
#' summary(GLMfit)
#'
#' ## MSB format
#' MSBfit <- glm.to.msbm (GLMfit)
#'
#' summary(MSBfit)
#'
glm.to.msbm <- function(object, refit = FALSE,
                        criterion = "ML", method = "msbm.fit",
                        start = object$coefficients,
                        control = object$control,
                        y = !is.null(object$y), frame = FALSE,
                        ...) {
  stopifnot(inherits(object, "glm"))
  stopifnot(identical(object$family$family, "binomial"))

  #* Initialize call and MSB model frame
  mcall <- mf <- object$call

  #* Transfer any free "offset" to the "formula" component
  formula <- object$formula
  offset <- object$offset
  if (!is.null(offset) & !is.null(object$call$offset)) {
    any.offset <- any(offset != 0)
    if (any.offset) {
      fcharacter <- as.character(formula[[3]])
      fcharacter <- paste0(fcharacter,
                           " + offset(",
                           object$call$offset, ")")
      formula[[3]] <- as.name(fcharacter)
      mcall$formula <- mf$formula <-
        formula <- as.formula(formula)
    }
  }

  #* Update the call
  mcall$family <- mcall$etastart <- mcall$mustart <-
    mcall$offset <- mcall$model <- mcall$x <-
    mcall$singular.ok <- mcall$contrasts <- NULL
  mcall$link <- object$family$link
  mcall[[1]] <- quote(msbreg)

  #* Build the MSB model frame
  m <- match(x = names(formals(msbreg)),
             table = names(mf),
             nomatch = 0L)
  mf <- mf[c(1L, m)]
  if (!is.character(method) && !is.function(method))
    stop("invalid 'method' argument")
  if (refit | identical(method, "model.frame")) {
    # etc <- names(mfexpand)
    mfexpand <- match.call(expand.dots = TRUE)
    mf$input.me <- mfexpand$input.me
    mf$stage.me <- mfexpand$stage.me
    mf$alpha.formula <- mfexpand$alpha.formula
    if (is.null(mf$alpha.formula))
      mf$alpha.formula <- ~ 0
    mf$lambda.formula <- mfexpand$lambda.formula
    if (is.null(mf$lambda.formula))
      mf$lambda.formula <- ~ 0
  }
  else {
    mf$alpha.formula <- ~ 0
    mf$lambda.formula <- ~ 0
  }
  mf$drop.unused.levels <- TRUE
  mf$verbose <- FALSE
  mf[[1L]] <- quote(msbreg::msbm.frame)

  mf <- eval(mf, parent.frame())
  if (identical(method, "model.frame"))
    return(mf)

  #* Update mcall
  mcall$alpha.formula <- mf$call$alpha.formula
  mcall$lambda.formula <- mf$call$lambda.formula

  #* Get Binomial link function
  link <- object$family$link
  eval(get.linkfun())

  if (!refit) {
    glm.slots <- c("coefficients", "fitted.values", "rank", "linear.predictors",
                   "deviance", "aic", "null.deviance", "iter",
                   "df.residual", "df.null", "y", "converged", "boundary",
                   "formula", "data", "method")
    keep <- names(object) %in% glm.slots
    out <- object[keep]
    yname <- colnames(model.frame(object))[1]
    out$ylab <- c(yname = yname,
                  ycolname = get.ycolname (yname))

    out$weights <- object$prior.weights
    out$sample.weights <- 1

    # Correct "linear.predictors" to have three columns (add "alpha" and "lambda")
    out$linear.predictors <- cbind(stage.1 = object$linear.predictors,
                                   alpha = -Inf,
                                   lambda = Inf)

    out$control <- control[names(control) %in%
                         names(formals(msbm.control))]
    out$control$criterion <- "ML"

    # Additional slots of for a "msbm" class object
    add.slots <- c("alpha.values", "lambda.values", "null.rank", "link", "bic",
                   "hic", "sample.weights", "deviance.resid", "nobs", "y.resid",
                   "criterion", "fisher", "singular", "start", "IC", "me",
                   "me.offset", "me.sd", "dispersion", "dims", "frame")

    out$alpha.values <- 0
    out$lambda.values <- 1
    out$null.rank <- mf$dims$nobs - out$df.null
    out$link <- link
    out$hic <- out$aic + 2 * mf$dims$npars * (log(log(mf$dims$nobs)) - 1)
    out$bic <- out$aic +     mf$dims$npars * (-2 + logb(mf$dims$nobs))
    out$deviance.resid <- stats::residuals(object, type = "deviance")
    out$nobs <- length(out$deviance.resid)
    if (identical(class(attr(mf, "na.action")), "omit")) {
      out$y.resid <- numeric(out$nobs) + NA
      out$y.resid[!((1:out$nobs) %in% attr(mf, "na.action"))] <-
        mf$y - na.omit(out$fitted.values)
    }
    else
      out$y.resid <- mf$y - out$fitted.values
    out$criterion <- "ML"
    out$fisher <- t(object$R) %*% object$R
    QR <- qr(out$fisher)
    out$fisher <- structure(out$fisher,
                            qr = structure(QR,
                                           solve = catch.conditions({
                                             solve.qr(QR)
                                           })$value),
                            logLik = logLik (object))

    out$singular <- out$rank < NCOL(out$fisher)
    object$call$start <- NULL
    out$IC <- out$me <- out$me.offset <- FALSE
    out$dispersion <- 1
    out$dims <- mf$dims
    out$frame <- if (frame) mf
    out$call <- mcall

    class(out) <- c(class(object)[!(class(object) %in%
                                      c("glm", "lm"))], "glm2msbm")
    return(out)
  }
  else {
    control <- control[names(control) %in%
                         names(formals(msbm.control))]
    if(missing(criterion) & !is.null(control$criterion)) {
      criterion <- control$criterion
    }
    control$criterion <- criterion
    control <- do.call("msbm.control", control)
    nstart <- length(start)
    nparnames <- length(mf$parnames)
    if (!(nstart %in% c(0, nparnames))) {
      start <- if (nstart < nparnames)
        c(start, numeric(nparnames - start))
      else
        start[1:nparnames]
    }

    #* Call fitter function
    fit <- eval(call(if (is.function(method)) "method" else method,
                     frame = mf, y = mf$y, start = start,
                     link = link, control = control))

    #* Add the response and the model frame if requested
    fit$y <- if (y) fit$y %||% drop.attr(mf$y)
    fit$na.action <- attr(mf, "na.action")

    return(structure(c(fit,
                       list(call = mcall, formula = formula,
                            link = link, dispersion = 1, data = object$data,
                            control = control, method = method,
                            dims = mf$dims,
                            frame = if (frame) mf)),
                     class = c(fit$class, c("glm2msbm", "msbm"))))

  }
}
