#' @rdname fisher.matrix
#' @exportS3Method msbreg::fisher.matrix
fisher.matrix.msbm <- function (object,
                                coefficients = object$coefficients,
                                sigma = NULL, ...) {
  # Default sigma if required
  eval(toget.msbm.sigma())

  # If the Fisher matrix is already there
  okay <- !is.null(object$fisher)
  if (!missing(coefficients)) {
    okay <- okay & all(object$coefficients == coefficients)
  }

  if (okay) {
    out <- object$fisher[1:object$dims$npars, 1:object$dims$npars]
    if (object$dispersion != 1) {
      out <- out * object$dispersion^2
    }
    out <- out / sigma^2

    return(structure(out,
                     sigma = sigma,
                     score = attr(object$fisher, "score.i"),
                     class = c("fisher.matrix",
                               class(out))))
  }

  mutation <- mutate.params (coefficients,
                             frame = model.frame.msbm (object),
                             link = link(object),
                             information = TRUE,
                             observed = FALSE)

  out <- mutation$EFinfo / (sigma^2)

  return(structure(out,
                   sigma = sigma,
                   score = mutation$score.i,
                   class = c("fisher.matrix",
                             class(out))))

}
### setMethod commented out to prevent the warning:
### Undocumented S4 methods:
###      generic 'fisher.matrix' and siglist 'msbm'
#setMethod ("fisher.matrix",
#           signature = "msbm",
#           definition = fisher.matrix.msbm)
