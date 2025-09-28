
#'
# @param x description
#'
# @exportS3Method print link
#'
#' @exportS3Method base::print
print.link <- function(x, ...) {
  cat("\nFamily:", x$family)
  cat("\nLink function:", x$link, "\n")
  if (!identical(x$mefamilyname, 'none')) {
    mefamilyname <- switch(x$mefamilyname,
                           norm   = "Gaussian",
                           logis  = "Logistic",
                           bridge = "Bridge",
                           lstab = "Log-Positive-Stable",
                           conjugate = "Conjugate",
                           x$mefamilyname)
    cat("Measurement error family:", mefamilyname, "\n\n")
  }
  invisible(x)
}
