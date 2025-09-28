#' Drop Attributes
#'
#' Remove Attributes from Objects. This is a generic with
#' a default method which remove the \code{'names'} attribute
#' from a vector input (the \code{'names'} attribute is set
#' to \code{NULL}) and \code{'dimnames'} from a matrix input.
#' For all other class object, the default method strips out
#' all \link{attributes} of the input \code{x}, except for the
#' \link{class} of \code{x}.
#'
#' @param x any \code{R} object.
#'
#' @param radical character, optional prefix to be used to build
#' a simplified \link{dim} attribute for \code{x}.
#'
#' @param ... further arguments passed to or from other methods.
#' Currently not used.
#'
#' @details
#' Methods can be written for \code{drop.attr} to remove
#' (from an object) certain attributes that are not (or no longer)
#' useful in a given context. This can be useful to remove attributes
#' before saving an object for a specific use for which, some attributes
#' are not of any use. For an object with a large amount of information
#' saved as attributes, this will make the object smaller in size and
#' may reduce the disk space required for saving.
#' The function may also be useful for displaying an object without
#' attributes that may be distracting or take too much space.
#'
#' Note that if mandatory information for a class of objects
#' is stored as an attribute rather than a slot, a call
#' to \code{drop.attr} would remove such information
#' and make the object an invalid member of the class.
#' Also note that the behavior of many objects when passed
#' through many generic functions depend on their
#' attributes, hence removing attributes will most likely
#' change the result of such generic functions, or result
#' into errors (methods defined for a class may require
#' some attributes during processing).
#'
#' @return For the default method,
#' an object identical to the input \code{x},
#' but with all attributes set to \code{NULL},
#' except for the \link{class} of \code{x}.
#'
#' For the \code{'fisher.matrix'} method,
#' a \code{dim} attribute is additionally conserved
#' when \code{radical} is not \code{NULL}.
#'
#' @aliases drop.attr
#' @export drop.attr
#'
#' @seealso \link{attributes} for details on object attributes.
#'
#' @examples
#' # Create a named numeric vector
#' Vec <- c(`(Intercept)` = 1, x1 = 3, x2 = 5)
#' Vec
#'
#' # Remove the 'names' attributes and display
#' drop.attr(Vec)
#'
drop.attr <- function (x, ...) {
  UseMethod("drop.attr")
}

drop.attr.default <- function(x, ...) {
  # If 'x' is a vector, just remove the 'names' attribute, if any
  if (inherits(x, "vector")) {
    names(x) <- NULL
    return(x)
  }

  if (inherits(x, "fisher.matrix")) {
    return(drop.attr.fisher.matrix(x, ...))
  }

  if (inherits(x, "matrix")) {
    dimnames(x) <- list(NULL, NULL)
    return(x)
  }

  # Otherwise, strip out all attributes
  classx <- class(x)
  attributes(x) <- NULL
  classx -> class(x)
  return(x)
}
setGeneric(name = "drop.attr",
           def = drop.attr.default)

#' @rdname drop.attr
#' @exportS3Method msbreg::drop.attr
drop.attr.fisher.matrix <- function(x, radical = NULL, ...) {
  if (!is.null(radical)) {
    stopifnot(is.character(radical))
    nn <- length(radical)
    if (nn == 1)
      nm <- paste0(radical, 1:NCOL(x))
    else if (nn == NCOL(x))
      nm <- radical
    dimnames(x) <- list(nm, nm)
  }
  else {
    dimnames(x) <- list(NULL, NULL)
  }

  return(x)

}
### setMethod commented out to prevent the warning:
### Undocumented S4 methods:
###      generic 'drop.attr' and siglist 'fisher.matrix'
#setMethod("drop.attr",
#          signature(x = "fisher.matrix"),
#          definition = drop.attr.fisher.matrix)
