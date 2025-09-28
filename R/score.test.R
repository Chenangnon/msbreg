# score.test
#' Score Test for Fitted Model Objects
#'
#' Performs score tests comparing a fitted (null) model to alternatives.
#'
#' @param object an object of a class with a \code{score.test} method.
#' Typically, a fit, for instance \code{"lm"} or \code{"glm"} class objects.
#'
#' @param ... further arguments passed to or from other methods.
#' This may include one or more objects specifying alternative models.
#' These may be objects of class "\link{formula}"s, or (non evaluated)
#' \link{call}s, or depending on the method for the class of \code{object},
#' objects of the same class as \code{object}.
#'
#' @param scope defines the range of alternative models to be examined.
#' This is used only when \code{...} does not include any model which nests
#' \code{object}.
#' This should be either a single formula, or a list containing components
#' \code{upper} and \code{lower}, both formulae.
#' See the details for how to specify the formulae and how they are used.
#'
#' @param trace integer indicating if details should be printed during
#' computations. Larger values may give more detailed information.
#'
#' @param steps the maximum number of null models to be considered.
#' The default is 1000 (essentially as many as required).
#' It is typically used to stop the process early when \code{...} is empty
#' and \code{scope} is a big model.
#'
#' @param ncores integer, number of threads for parallel processing.
#'
#' @details
#' This is a generic function to perform the Rao score test (or modified versions)
#' on a model fit.
#' The score test provides a faster alternative to the Wald (obtained from
#' a call to \link[stats]{summary.glm}) and the likelihood ratio tests.
#' Indeed, only the null model need be fitted.
#' The argument \code{...} allows to precisely indicate the desired alternatives
#' to the model \code{object}.
#' If instead, a range of alternatives is desired, the argument \code{scope}
#' allows to delimit it.
#'
#' The details given below are for an \code{object} of class "\link{lm}" or "\link{glm}".
#' Note that all variables in all candidate models should be available from the
#' null model fit \code{object}, unless \code{scope} contains fitted models.
#'
#' If alternative models are specified through the argument \code{...},
#' the function only retains those models which contain \code{object} (the null model)
#' as a nested model. It then computes the score statistics for the null model
#' versus each of the alternatives.
#'
#' If no alternative model is given through the argument \code{...}, or none of the
#' provided alternatives is larger than \code{object}, the function considers
#' the argument \code{scope}.
#' The right-hand-side of the \code{upper} component of \code{scope} is always
#' considered as an alternative model, provided it is larger than \code{object}.
#' If \code{scope} is a single formula, it specifies the \code{upper} component,
#' and the \code{lower} model is \code{object}.
#' Any model smaller than \code{object} is replaced by \code{object}.
#' If both the \code{lower} and \code{upper} components indicated by \code{scope}
#' are identical to \code{object}, no computation is performed and a zero statistic
#' value is returned, comparing \code{object} with itself.
#'
# The function \code{scoretable} generates a summary of one or many
# \code{score.test} results.
#'
#' @return An object that inherits from the class ‘htests’,
#' a list with components:
#'
#' Write 'print.htests' to display a table of test results (seee 'stats:::print.htest')
#'
#' The \code{print} method for \code{"score.test"} class objects calls
#' \code{scoretable}.
#'
#' @aliases score.test.default
#' @export score.test
# @export score.test.default
#'
#' @seealso See \link[msbreg]{profile.msbm} for computing log-likelihood profiles
#' for MSB model fits.
#'
# @note
# Both function were inspired by the (non-generic) \link[secr]{score.test}
# and \link[secr]{score.table} of the package *secr* \insertCite{murray2025secr}{msbreg}.
#'
#' @seealso \link{score.test.msbm} for the \code{score.test} method for
#' fitted Multistage Binomial (MSB) models.
#'
score.test <- function(object, ...) {
  UseMethod("score.test")
}

score.test.default <- function(object, ...) {
  # Piece of cake if 'lm' or 'glm'
  if (inherits(object, "glm") | inherits(object, "lm")) {
    dname <- deparse1(substitute(object))
    out <- score.test.glm.intern (object, ..., dname = dname)
    out$call <- match.call()
    out$call[[1L]] <- quote(msbreg::score.test)

    out <- structure(list(out),
                     class = "score.list")

    return(out)
  }

  if (inherits(object, "msbm")) {
    dname <- deparse1(substitute(object))
    out <- score.test.msbm.intern (object, ..., dname = dname)
    out$call <- match.call()
    out$call[[1L]] <- quote(msbreg::score.test)

    out <- structure(list(out),
                     class = "score.list")

    return(out)
  }

  stop(paste0("There is no 'score.test' method for objects of class ",
              paste0(class(object), ".", collapse = ", ")))

}
setGeneric(name = "score.test",
           def = score.test.default)

#' @rdname score.test
#' @exportS3Method msbreg::score.test
score.test.lm <- function(object, ..., scope, steps = 1000, trace = FALSE, ncores = NULL) {
  dname <- deparse1(substitute(object))
  out <- score.test.glm.intern (object, ..., scope = scope, steps = steps,
                                 dname = dname, trace = trace, ncores = ncores)
  out$call <- match.call()

  out
}

#' @rdname score.test
#' @exportS3Method msbreg::score.test
score.test.glm <- score.test.lm

score.test.glm.intern <- function(object, ..., scope, steps = 1000,
                                  dname = deparse1(substitute(object)),
                                  trace = FALSE, ncores = NULL) {
  if(inherits(object, "lm"))
    scorefun <- function(formula) {
      coefficients <- object$coefficient

    }
  else
    scorefun <- function(formula) {
      coefficients <- object$coefficient

    }


}
