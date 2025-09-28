
#' @rdname fisher.matrix
#' @exportS3Method msbreg::fisher.matrix
fisher.matrix.msbm.frame <- function (object, coefficients,
                                      link = "logit",
                                      sigma = NULL, ...) {

  # Get Binomial link function
  eval(get.linkfun())

  # Default sigma if required
  eval(toget.msbm.sigma())

  mutation <- mutate.params (theta = coefficients,
                             frame = object,
                             link = link,
                             score = attr(object, "response") == 1,
                             information = TRUE,
                             observed = FALSE)

  out <- mutation$EFinfo / (sigma^2)

  return(structure(out,
                   sigma = sigma,
#                   logLik = if (attr(object, "response")) mutation$loglike.i,
                   score = if (attr(object, "response")) mutation$score.i,
                   class = c("fisher.matrix",
                             class(out))))

}

### setMethod commented out to prevent the warning:
### Undocumented S4 methods:
###      generic 'fisher.matrix' and siglist 'msbm.frame'
#setMethod ("fisher.matrix",
#           signature = "msbm.frame",
#           definition = fisher.matrix.msbm.frame)

toget.msbm.sigma <- function() {
  expression({
    if (!missing(sigma)) {
      if (!is.null(sigma)) {
        stopifnot(is.numeric(sigma))
        stopifnot(sigma > 0)

        if (sigma != 1)
          warning(paste0("non-unit dispersion supplied for the binomial family"))
      }
      else {
        sigma <- 1
      }
    }
    else {
      sigma <- 1
    }
  })
}
