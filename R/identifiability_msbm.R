# @export identifiability.msbm
# @export identifiability.msbm.frame

#' Practical Identifiability of a Model
#'
#' These are methods of the generic \link[msbreg]{identifiability} for
#' Multistage Binomial (MSB) model frames (class \code{"msbm.frame"} from
#' \link[msbreg]{msbm.frame}) and model fits (class \code{"msbm"} from
#' \link[msbreg]{msbreg}). They return a logical indicating the identifiability
#' of model parameters (coefficients for all contrasts in the model).
#'
#' @param object an object of class \code{"msbm"} (model fit) or
#' \code{"msbm.frame"} (model frame).
#'
#' @param continuous.levels numeric scalar, minimum number of measured/observed
#' levels required for a numerical predictor in the \code{"msbm.frame"} to be
#' considered as a continuous variable during identifiability analysis.
#' The default value is \code{continuous.levels = 4}.
#'
#' For an \code{object} of class \code{"msbm"}, the argument
#' \code{continuous.levels}is used only if \code{full = TRUE}.
#'
#' @param full logical, should the identifiability of the model frame be
#' checked? The default \code{FALSE} means that for an object of class
#' \code{"msbm"}, only the rank of the Fisher information matrix is checked
#' and identifiability is declared if the information matrix is of full rank.
#' Setting \code{full = TRUE} means that the identifiability of the model frame
#' is also checked, and the result included in the output.
#'
#' @param ... further arguments passed to or from other methods.
#' For the default method, optional arguments for \link[stats]{model.frame}
#' such as \code{data}, \code{subset}, \code{na.action},
#'  \code{drop.unused.levels}, and \code{xlev}.
#'
#' @details
#' The \code{identifiability} method for \code{"msbm.frame"} class objects
#' assesses each model component separately (i.e. each individual linear
#' predictor), and then checks cross-identifiability (non-interchangeability
#' of the linear predictors).
#'
#' For an \code{object} of class \code{"msbm"}, the method
#' checks that the Fisher information matrix of the fit is of full rank.
#' If \code{full = TRUE}, the method also obtains a \code{"msbm.frame"}
#' from \code{object} (by extracting a \code{$frame} component if available,
#' building one otherwise), and then uses the method for \code{"msbm.frame"}
#' class objects: identifiability is declared if both processes indicate
#' identifiability.
#'
#' @return a logical value: is the model represented by the first argument
#' identifiable? The returned value inherits from class \code{"id.analysis"}
#' which has simple \link{print} and \link{summary} methods.
#'
#' The logical output indicates practical identifiability (the model may be
#' structurally identifiable while the output is \code{FALSE}).
#'
#' For an \code{object} of class \code{"msbm.frame"} the method returns
#' the attributes \code{separated} and \code{aliased} (see the default
#' \link[msbreg]{identifiability} method for details),
#' each being either \code{NULL}, or a list with vector elements
#' corresponding to different model components.
#' In addition, the following attributes are returned:
#' \describe{
#' \item{\code{summary}}{ a matrix giving the collection of individual
#' vector \code{summary} attributes (as defined above for the default method)
#' for different model components (in columns);}
#' \item{\code{continuous}}{ a two-row matrix giving the number of continuous
#' predictors (first row) only found in each model stage (columns), and
#' the presence of a continuous offset (second row) in each stage;}
#' \item{\code{resid}}{ a matrix giving the above defined \code{resid}
#' attribute for each model component (in rows);}
#' \item{\code{epsilon}}{ scalar as defined above.}
#' }
#'
#' For an \code{object} of class \code{"msbm"} the above attributes are
#' by default \code{NULL} (not returned), but are available when
#' \code{full = TRUE}. The output also has attributes \code{rank}
#' (the rank of the Fisher information matrix) and \code{npars} (the
#' number of model parameters).
#'
#'
#' @export identifiability
# @export identifiability.msbm.frame
#'
#' @seealso See \link[msbreg]{test.separation} for testing quasi/complete
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
#' # Build a MSB model frame with three stages
#' library(msbreg)
#' mframe <- msbm.frame (formula = y ~ offset(x1) | f1 +  x2 | x3,
#'                       data = df)
#'
#' # Test for identifiability in the MSB model frame
#' id <- identifiability (mframe)
#' id
#' summary(id)
#'
#' # Same test with the requirement: a variable "is continuous"
#' # if it has at least 15 distinct observed values
#' identifiability (mframe, data = df, continuous.levels = 15)
#'
#' # Test for identifiability in the MSB model fit
#' mfit <- msbreg (formula = y ~ offset(x1) | f1 +  x2 | x3,
#'                 data = df)
#'
#' identifiability (mfit)
#'

# @rdname identifiability
#' @exportS3Method msbreg::identifiability
#'
#' @import stats
identifiability.msbm <- function (object, full = FALSE, continuous.levels = 4, ...) {

  if (full) {
    id <- identifiability.msbm.frame(model.frame (object),
                                     continuous.levels = continuous.levels, ...)
  }
  else {
    id <- TRUE
  }

  if (is.null(object$rank) | is.na(object$rank)) {
    mrank <- matrix.rank (fisher.test(object, ...))
  }
  else {
    mrank <- object$rank
  }

  id[1] <- id[1] & (mrank == object$dims$npars)
  attr(id, "npars") <- object$dims$npars
  attr(id, "rank") <- mrank

  class(id) <- 'id.analysis'

  id
}

#setMethod ("identifiability",
#           signature = "msbm",
#           definition = identifiability.msbm)
