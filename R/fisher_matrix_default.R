# fisher.matrix
#' Fisher Information Matrix
#'
#' Compute the Fisher information matrix of a model (represented by
#' a fit, or a \code{model.frame} like object) given a parameter vector.
#'
#' @param object an object of a class with a \code{fisher.matrix} method.
#' Typically, a fit (for instance \code{"lm"} and \code{"glm"} class
#' objects) or a \code{model.frame} like object.
#'
#' @param sigma a numeric scalar (positive real), optional
#' dispersion parameter in the Generalized Linear Model framework
#' (e.g. residual standard deviation in Linear Models).
#'
#' When \code{sigma} is missing or \code{NULL} (the default), it is
#' inferred from the model frame if the model frame includes
#' a response, otherwise, \code{sigma = 1} is used.
#'
#' @param coefficients a numeric vector of all MSB model parameters.
#'
#' @param link optional, a \code{"link"} object (as returned by
#' \link[msbreg]{link}); defaults to \code{logit()}
#' (see \link[msbreg]{logit}).
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' This is a generic function to extract the Fisher information
#' matrix from a model fit or a model frame (given a vector of
#' model parameters when required).
#'
#' The default method can summarily handle any \code{object} of a class
#' possessing a \code{model.frame} method that returns a \code{"model.frame"}
#' class output (when applied to the input \code{object}). This includes for
#' instance \link{formula}, \link{terms}, ... objects, in addition to
#' proper \link{model.frame} objects. The returned information matrix
#' corresponds to a Linear Model based on the corresponding
#' \code{"model.frame"}. The default method also handles
#' objects of classes \link{lm} and \link{glm}.
#'
#' For the default method, the \code{...} arguments can include \code{data},
#' \code{subset}, \code{na.action}, \code{drop.unused.levels}, and \code{xlev};
#' all passed to \link[stats]{model.frame} when supplied.
#'
#' Objects of class \code{"msbm"} (as built by \link[msbreg]{msbreg})
#' have a dedicated \code{fisher.matrix} method  which takes a Multistage
#' Binomial (MSB) model fit and extracts the Fisher information matrix.
#' When called on an \code{object} of class \code{"msbm.frame"},
#' \code{fisher.matrix} takes the additional arguments:
#'
#' \itemize{
#' \item \code{coefficients} - *required*, a numeric vector of all MSB model
#' parameters; and
#'
#' \item \code{link} - *optional*, a \code{"link"} object (as returned by
#' \link[msbreg]{link}); defaults to \code{logit()}
#' (see \link[msbreg]{logit}).
#' }
#'
#' @return a \eqn{p \times p} matrix where \eqn{p} is the number of
#' model parameters, typically the number of columns of the matrix
#' of contrasts from a model frame.
#'
#' The output inherits from class \code{"fisher.matrix"} (which
#' encompasses symmetric non-negative definite matrices) and has
#' attribute:
#'
#' \item{\code{sigma}}{ numeric scalar, the used \code{sigma} argument.}
#'
#' For an \code{object} of class \code{"msbm.frame"} with a non \code{NULL}
#' response component (\code{attr(object, "response") == 1}), the output
#' has an additional attribute:
#'
# \item{\code{logLik}}{ numeric vector of individual log-likelihood
# contributions of each data point in the MSB model frame (so that
# \code{sum($logLik)} is the log-likelihood function evaluated at
# the used vector of model parameters);}
#' \item{\code{score}}{ numeric matrix with \eqn{p} columns (corresponding
#' to model parameters) and \eqn{n} rows (corresponding to data points in
#' the MSB model frame); the \eqn{i}th row of \code{score} gives the
#' contribution of the \eqn{i}th data point to the model score (first
#' derivative of the log-likelihood function with respect to all model
#' parameters); this means that \code{colSums(score)} gives the model
#' score vector.}
#'
#' There may be other attributes depending on the method used: see
#' the appropriate documentation.
#'
#' There is a simple \code{print} method for \code{"fisher.matrix"} class objects.
#'
#' @aliases fisher.matrix.default
#' @export fisher.matrix
# @export fisher.matrix.default
# @export fisher.matrix.msbm.frame
#'
#' @seealso See \link[msbreg]{identifiability} for testing identifiability
#' of regression models.
#'
#' @examples
#' ##* Generate a dataset
#' set.seed(10)
#' df <- data.frame(y = rbinom(15, 1, 0.5),
#'                  x1 = rnorm(15),
#'                  x2 = rpois(15, 5))
#'
#' ##* Fisher information matrix from a model frame (linear model assumed)
#' library(msbreg)
#' mf <- model.frame(y ~ x1 + x2, data = df)
#' fisher.matrix (mf) # dispersion estimated using the response 'y'
#'
#' fisher.matrix (mf, sigma = 1) # dispersion set to one
#'
#' mf <- model.frame( ~ x1 + x2, data = df)
#' fisher.matrix (mf) # sigma = 1 if 'y' is missing (and sigma is not specified)
#'
#' ##* Fisher information matrix for a linear model fit
#' lm1 <- lm (y ~ x1 + x2, data = df)
#' fisher.matrix (lm1)
#'
#' ##* Fisher information matrix for a generalized linear model fit
#' glm1 <- glm (y ~ x1 + x2, family = binomial(), data = df)
#' fisher.matrix (glm1)
#'
#' ##* Fisher information matrix for a MSB model
#' data(test1data)
#'
#' #* Build a MSB model frame
#' msbframe <- msbm.frame (formula = y/Total ~ x1 + offset(off1) | x2,
#'                         # minimum probability set to zero
#'                         alpha.formula = ~ 0,
#'                         # maximum probability is constant, allowed to be < 1
#'                         lambda.formula = ~ 1,
#'                         weights = Total,
#'                         data = test1data, frames = TRUE)
#'
#' msbframe
#'
#' #* Model parameter vector
#' theta0 <- attr(test1data$y,"theta")
#' theta0
#'
#' #* Fisher Information Matrix
#' FImat <- fisher.matrix (msbframe, coefficients = theta0, link = "logit")
#'
#' #* Print first two parameters
#' print(FImat[1:2, 1:2], digits = 5)
#'
#' #* Strip out parameter names to print them all
#' drop.attr(FImat, radical = "P")
#'
fisher.matrix <- function(object, ...) {
  UseMethod("fisher.matrix")
}

fisher.matrix.default <- function(object, ...) {
  # Piece of cake if 'lm' or 'glm'
  if (inherits(object, "glm")) {
    return(fisher.matrix.glm (object, ...))
  }
  if (inherits(object, "lm")) {
    return(fisher.matrix.lm (object, ...))
  }
  if (inherits(object, "msbm")) {
    return(fisher.matrix.msbm (object, ...))
  }
  if (inherits(object, "msbm.frame")) {
    return(fisher.matrix.msbm.frame (object, ...))
  }

  # Added HERE because "fisher.matrix" was not being taken as a generic
  # Got this warning: In .S3methods(generic.function, class, envir, all.names = all.names,  :
  # function 'fisher.matrix' appears not to be S3 generic; found functions that look like S3 methods
  # See sloop::s3_methods_generic("fisher.matrix")
  # Pretty sure this is related to "msbm.frame",
  # Warning during installation: no definition for class “msbm.frame”
  #  if (inherits(object, "msbm.frame")) {
  #    return(fisher.matrix.msbm.frame (object, ...))
  #  }

  dotnames <- catch.conditions({
    ...names() # Seems absent in R/4.0.2
  })$value
  if (any(class(dotnames) %in% c("simpleError", "error",
                                 "condition", "try-error"))) {
    dotnames <- names(list(...))
    dotok <- FALSE
  }

  if ("data" %in% dotnames) {
    if (inherits(object, "model.frame")) {
      return(fisher.matrix.model.frame (object, ...))
    }
    else {
      return(fisher.matrix.model.frame (model.frame(object, ...),
                                         ...))
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
            return(fisher.matrix.model.frame (object, data = data, ...))
          }
          else {
            return(fisher.matrix.model.frame (model.frame(object, data = data, ...),
                                              data = data, ...))
          }
        }
        else {

          if (inherits(object, "model.frame")) {
            return(fisher.matrix.model.frame (object, ...))
          }
          else {
            return(fisher.matrix.model.frame (model.frame(object, ...),
                                              ...))
          }
        }
      }
      else {
        if (inherits(object, "model.frame")) {
          return(fisher.matrix.model.frame (object, ...))
        }
        else {
          return(fisher.matrix.model.frame (model.frame(object, ...),
                                             ...))
        }
      }
    }
    data <- as.data.frame(object)
      if (is.null(data)) {
        return(fisher.matrix.model.frame (model.frame(object, ...),
                                          data = environment(formula(object)),
                                          ...))
      }
      else {
        return(fisher.matrix.model.frame (model.frame(object, ...),
                                          data = data, ...))
      }

  }

}
setGeneric(name = "fisher.matrix",
           def = fisher.matrix.default)

# @rdname fisher.matrix
#' @exportS3Method msbreg::fisher.matrix
fisher.matrix.formula <- function (object, sigma = NULL, data, ...) {
  mframe <- model.frame(object, data, ...)

  fisher.matrix.model.frame (mframe, sigma = sigma, data = data, ...)
}
#.S3method("fisher.matrix", "formula", "fisher.matrix.formula")

# @rdname fisher.matrix
#' @exportS3Method msbreg::fisher.matrix
fisher.matrix.terms <- fisher.matrix.formula

# @rdname fisher.matrix
#' @exportS3Method msbreg::fisher.matrix
fisher.matrix.model.frame <- function (object, sigma = NULL, data, ...) {
  # Extract y and the design matrix
  eval(extract.model.frame())

  # Final result
  return(.fisher.matrix.x.y (x = X, y = Y,
                             sigma = sigma,
                             weights = ntrials,
                             offset = Offset, ...))
}

#' @exportS3Method base::print
print.fisher.matrix <- function(x, n = Inf,
                                digits = getOption("digits") - 3, ...) {
  cat("\n *** Fisher Information Matrix *** \n")

  out <- as.data.frame(as.matrix(x))
  if (is.finite(n)) {
    n <- min(n, NCOL(out))
    out <- out[1:n, 1:n]
  }

  print(out, digits = digits, ...)

  if (!is.null(attr(x, "sigma"))) {
    cat(" *Dispersion (sigma): ",
        round(attr(x, "sigma"), digits, ...),
        "\n")
  }
}
