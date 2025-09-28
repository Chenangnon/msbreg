#' Extracting the Model Frame from a Fit
#'
#' \code{model.frame.msbm} is a method of the \code{S3} generic function
#' \link[stats]{model.frame} for class \code{"msbm"}. It returns an object of
#' class \code{"msbm.frame"} which is a list of \link{matrix} and dictionary
#' for extracting the model components and variables.
#'
#' @param formula an object of class \code{"msbm"} (as returned by \link{msbreg}).
#'
#' @param ... a mix of further optional arguments such as \code{alpha.formula},
#' \code{lambda.formula}, \code{me.offset}, \code{me.model}, \code{data},
#' \code{weights}, \code{sample.weights}, \code{subset}, \code{na.action},
#' \code{drop.unused.levels}, and \code{frames} to pass to \link{msbm.frame}
#' (when that call is required).
#'
#' @details
#' For an object of fitted-model class \code{"msbm"}, the method returns either
#' the saved model frame used when fitting the model (if any, when selected by
#' argument \code{frame = TRUE}) or pass the call used to fit the model on to
#' the function \link[msbreg]{msbm.frame}.
#'
# @name model.frame
#'
#' @examples
#' # Generate a dataset including a few missing values
#' # (binary response unrelated to the predictors)
#' set.seed(10)
#' dfex <- data.frame(y = c(NA, rbinom(14, 1, 0.5)),
#'                    # First continuous predictor
#'                    x1 = c(rnorm(7), NA, rnorm(7)),
#'                    # known standard deviation of x1
#'                    SDx1 = c(rexp(7, 1), NA, rexp(7, 1)),
#'                    # Second continuous predictor
#'                    x2 = c(rnorm(14), NA),
#'                    # known standard deviation of x2
#'                    SDx2 = c(rexp(14, 1), NA),
#'                    # Numeric, but discrete predictor
#'                    x3 = rpois(15, 5),
#'                    # Factor with three levels
#'                    f1 = rep(c('a', 'b', 'c'),
#'                             length.out = 15))
#'
#' # Build a MSB model frame with two stages
#' mframe <- msbm.frame (formula = y ~ x1 * f1| x2 * x3,
#'                       me.offset = list(x1 = ~ SDx1,
#'                                        x2 = ~ SDx2),
#'                       data = dfex, frames = TRUE)
#'
#'
#'
#' # Create a pseudo "msbm" class object (only has a $frame component)
#' pseudofit <- structure(list(frame = mframe),
#'                        class = 'msbm')
#'
#' # Extract a model frame from the pseudo-fit
#' model.frame (pseudofit)
#'
#' # Create a pseudo "msbm" class object (has $call and $data components)
#' pseudofit <- structure(list(call = mframe$call, data = dfex),
#'                        class = 'msbm')
#'
#' # Extract a model frame from the pseudo-fit
#' model.frame (pseudofit)
#'
#' @exportS3Method stats::model.frame
model.frame.msbm <- function(formula, ...) {
  # Extract the msbm.frame (if any) from the fit
  mf <- formula$frame

  if (inherits(mf, "msbm.frame")) {
    return(mf)
  }

  # Build the msbm.frame
  mf <- formula$call
  m <- match(x = c("formula", "input.me", "stage.me",
                   "alpha.formula", "lambda.formula",
                   "data", "weights", "sample.weights", "subset",
                   "na.action", "drop.unused.levels", "frames"),
             table = names(mf),
             nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(msbreg::msbm.frame)
  mcall <- mf
  mf$verbose <- FALSE
  if (is.null(mf$data))
    mf$data <- formula$data

  mframe <- catch.conditions({
    eval(mf, envir = formula$data, enclos = parent.frame())
  })$value

  if (any(class(mframe) %in% c("simpleError", "error",
                               "condition", "try-error"))) {
    mframe <- catch.conditions({
      eval(mf, envir = parent.frame())
    })$value

    if (any(class(mframe) %in% c("simpleError", "error",
                                 "condition", "try-error"))) {
      mframe <- catch.conditions({
        eval(mf, environment(formula$formula))
      })$value
    }

    if (any(class(mframe) %in% c("simpleError", "error",
                                 "condition", "try-error"))) {
      mframe <- catch.conditions({
        eval(mf, environment(formula))
      })$value
    }

    if (any(class(mframe) %in% c("simpleError", "error",
                                 "condition", "try-error"))) {
      mframe <- catch.conditions({
        eval(mf, environment(formula$formula))
      })$value
    }

    if (any(class(mframe) %in% c("simpleError", "error",
                                 "condition", "try-error"))) {
      stop(mframe)
    }
  }

  mframe$call <- mcall

  return(mframe)
}
