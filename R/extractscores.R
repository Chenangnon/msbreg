
#' @rdname score.test.msbm
#' @export extractscores
extractscores <- function(object, which = c(1,3,5,6), ...) {
  out <- sapply(object,
                FUN = function(x) {
                  c(unlist(x$scores), x$H1objective)[which]
                  })

  nm <- rownames(out)
  if (!is.null(nm)) {
    nm <- paste0("sc.", nm)
  }
  if (NCOL(out) == 1) {
    out <- c(out)
    names(out) <- nm
  }
  else {
    out <- t(out)
    colnames(out) <- nm
  }

  out
}
