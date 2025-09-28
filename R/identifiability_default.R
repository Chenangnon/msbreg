# Define a generic 'identifiability' function
#
#' Practical Identifiability of a Model
#'
#' \code{identifiability} (a generic function) takes in a model frame or a
#' fitted-model object and returns a logical indicating the identifiability
#' of model parameters (coefficients for all contrasts in the model). The
#' output generally has some \link{attributes} which are statistics for
#' diagnosing the source of potential non-identifiability. The default
#' method handle \code{"model.frame"} class objects
#' (\link[stats]{model.frame}), and objects of related classes such as
#' \code{"terms"}, \code{"formula"} and \code{"glm"}.
#'
#' @param object an object of a class with an \code{identifiability} method.
#'
#' @param ... further arguments passed to or from other methods.
#' For the default method, optional arguments for \link[stats]{model.frame}
#' such as \code{data}, \code{subset}, \code{na.action},
#'  \code{drop.unused.levels}, and \code{xlev}.
#'
#' @details
#' This is a generic function with a default method for \link[stats]{model.frame}
#' class objects. The default method takes a standard \code{"model.frame"}
#' object and check the rank of the corresponding design matrix, and the
#' separability of the response.
#'
#' The default method can summarily handle any \code{object} of a class
#' possessing a \code{model.frame} method that returns a \code{"model.frame"}
#' class output (when applied to the input \code{object}). This includes for
#' instance \link{formula}, \link{lm}, \link{glm} ... objects (note that
#' \link{lm} and \link{glm} objects have \code{$qr} list elements that
#' are sufficient for identifiability analysis in \code{lm}s, but not
#' sufficient for \code{glm}s), in addition to proper \code{"model.frame"} objects.
#'
#' @return a logical value: is the model represented by the first argument
#' identifiable? The returned value inherits from class \code{"id.analysis"}
#' which has simple \link{print} and \link{summary} methods.
#'
#' The logical output indicates practical identifiability (the model may be
#' structurally identifiable while the output is \code{FALSE}).
#'
#' For the default method, the returned value has the following
#' attributes (which may be missing or \code{NULL}):
#'
#' \item{\code{summary}}{ a named vector of information on the singularity and
#' separability of the design matrix, and details from both singularity
#' test and separation test (see \link[msbreg]{test.separation}):}
#'
#' \itemize{
#' \item \code{singular}: binary scalar, is the design matrix singular?
#' \item \code{rank.all}: numeric, rank of the design matrix, computed for
#' all observations;
#' \item \code{rank.successes}: numeric, rank of the design matrix, computed for
#' non-zero observations (\code{y > 0});
#' \item \code{rank.failures}: numeric, rank of the design matrix, computed for
#' zero observations (\code{y = 0});
#' \item \code{separable}: binary scalar, does the response in the model separate
#' any of the (group of) predictors included in the model?
#' \item \code{converged}: binary scalar, did the separation detection
#' algorithm converge?
#' \item \code{niter}: number of iterations used by the employed algorithm;
#' }
#'
#' \item{\code{aliased}}{ \code{NULL}, or a list of integer vectors
#' (one for each model component) indicating columns of the design
#' matrix \code{x} with aliased coefficients (these columns are
#' removed before any separation related calculation, leading
#' to a non-singular reduced design matrix, used for the subsequent
#' separation analysis);}
#' \item{\code{separated}}{ \code{NULL}, or a list of integer vectors
#' (one for each model component) indicating which observations
#' (if any) are separated by the response;}
#' \item{\code{resid}}{a named vector of additional information on the
#' achieved convergence threshold during separation analysis, i.e. the
#' absolute values of iterated weighted least-square residuals:}
#'
#' \itemize{
#' \item \code{min}: observed minimum of the absolute residuals;
#' \item \code{max}: observed maximum of the absolute residuals
#' (convergence is determined by \code{max < epsilon});
#' }
#'
#' \item{\code{epsilon}}{ the tolerance value used for the convergence
#' of the separation detection algorithm.}
#'
#' There may be other attributes depending on the method used: see
#' the appropriate documentation.
#'
#' Note that if the attribute \code{summary[1] = 1}, the \code{aliased}
#' attribute indicates columns that have linear dependence with
#' the other columns and have been removed to obtain a non-singular
#' reduced design matrix, before separation related calculations
#' were performed.
#' Likewise, if the attribute \code{summary[5] = 1}, the \code{separated}
#' attribute indicates observations that are separated by the response.
#'
#' @aliases identifiability.default
#' @export identifiability
#' @export identifiability.default
#'
#' @seealso See \link[msbreg]{identifiability.msbm} for methods
#' for multistage binomial model frame and fit.
#' See also \link[msbreg]{test.separation} for testing quasi/complete
#' separation in binary data regression.
#'
#' @examples
#' # Generate a dataset
#' # (binary response unrelated to the predictors)
#' set.seed(10)
#' df <- data.frame(y = rbinom(15, 1, 0.5),
#'                  # First continuous predictor
#'                  x1 = rnorm(15),
#'                  # Second continuous predictor
#'                  x2 = rnorm(15),
#'                  # Numeric, but discrete predictor
#'                  x3 = rpois(15, 5),
#'                  # Factor with three levels
#'                  f1 = rep(c('a', 'b', 'c'),
#'                           length.out = 15))
#'
#' # Test for identifiability given a model formula
#' library(msbreg)
#' identifiability (y ~ x1 + f1, data = df)
#'
#' # Test for identifiability in a standard binomial model
#' glm1 <- glm(y ~ x1 + f1, family = binomial(), data = df,
#'             na.action = na.omit)
#' identifiability (glm1, data = df)
#'
identifiability <- function(object, ...) {
  UseMethod("identifiability")
}

identifiability.default <- function(object, ...) {
  if (inherits(object, "msbm.frame")) {
    return(identifiability.msbm.frame (object, ...))
  }

  if (inherits(object, "msbm")) {
    return(identifiability.msbm (object, ...))
  }

  dotnames <- ...names()
  if ("data" %in% dotnames) {
    if (inherits(object, "model.frame")) {
      return(identifiability.model.frame (object, ...))
    }
    else {
      return(identifiability.model.frame (model.frame(object, ...), ...))
    }
  }
  else {
    mcall <- object$call
    if (!is.null(mcall)) {
      data <- mcall$data
      if (!is.null(data)) {
        Env <- environment(mcall)
        if (!is.null(Env)) {
          data <- catch.conditions({
            eval(mcall$data, envir = Env, enclos = parent.frame())
          })$value

          if (any(class(data) %in% c("simpleError", "error",
                                    "condition", "try-error"))) {
            data <- catch.conditions({
              eval(mcall$data, envir = parent.frame())
            })$value

            if (any(class(data) %in% c("simpleError", "error",
                                       "condition", "try-error"))) {
              data <- NULL
            }
          }
        }
        else {
          data <- catch.conditions({
            eval(mcall$data, envir = parent.frame())
          })$value

          if (any(class(data) %in% c("simpleError", "error",
                                     "condition", "try-error"))) {
            data <- NULL
          }
        }

        if (!is.null(data)) {
          if (inherits(object, "model.frame")) {
            return(identifiability.model.frame (object, data = data, ...))
          }
          else {
            return(identifiability.model.frame (model.frame(object, data = data, ...),
                                                data = data, ...))
          }
        }
        else {

          if (inherits(object, "model.frame")) {
            return(identifiability.model.frame (object, ...))
          }
          else {
            return(identifiability.model.frame (model.frame(object, ...), ...))
          }
        }
      }
      else {
        if (inherits(object, "model.frame")) {
          return(identifiability.model.frame (object, ...))
        }
        else {
          return(identifiability.model.frame (model.frame(object, ...), ...))
        }
      }
    }
    data <- as.data.frame(object)
    if (is.null(data)) {
        return(identifiability.model.frame (object,
                                            data = environment(formula(object)),
                                            ...))
      }
      else {
        return(identifiability.model.frame (model.frame(object, ...),
                                            data = data, ...))
      }
  }
}
setGeneric(name = "identifiability",
           def = identifiability.default)

#' @exportS3Method msbreg::identifiability
identifiability.formula <- function (object, data, ...) {
  mframe <- model.frame(object, data, ...)

  identifiability.model.frame (mframe, ...)

}

#' @exportS3Method msbreg::identifiability
identifiability.terms <- identifiability.formula

#' @exportS3Method msbreg::identifiability
identifiability.model.frame <- function (object, data, ...) {
  # Extract y and the design matrix
  eval(extract.model.frame())

  # Final result
  return(.identifiability.x.y (x = X, y = Y,
                               weights = ntrials,
                               offset = Offset,
                               method = "iterative_rectifier", ...))
}

#' @exportS3Method base::print
print.id.analysis <- function(x, ...) {
  print(as.logical(x), ...)
}

#' @exportS3Method base::summary
summary.id.analysis <- function(object, ...) {
  smry <- cbind(attr(object, "summary"))

  if (is.null(smry)) {
    out <- object[1]
    print(object, ...)
  }
  else {
    print(smry, ...)
    out <- smry
  }

  if(object)
    cat("     * Identifiable binary data model \n")
  else
    cat("     * Non identifiable binary data model \n")

  return(invisible(out))
}
