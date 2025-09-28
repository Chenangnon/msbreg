#' @exportS3Method base::print
print.msbm.frame <- function (x, n = 6L, ...) {
  cat("\n *** Multistage binomial model frame *** \n")
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")
  eval(print_frame_msbm())
}

print_frame_msbm <- function() {
  expression({
    cat("\n Matrix of stages (rows) and terms (columns)  \n")
    print(cbind(x$tables$stage.terms, x$tables$stage.offsets), ...)

    nrows <- min(n, x$dims$nobs)
    out.main <- x$input.matrix[1:nrows, , drop = FALSE]
    #out.main <- cbind(x$y, x$input.matrix)[1:nrows, , drop = FALSE]
    #colnames(out.main)[1] <- x$auxy$yname
    attr(out.main, 'na.action') <- attr(x$input.matrix, 'na.action')
    cat("\n Head of the matrix of contrasts (", NROW(x$input.matrix), "rows * ",
        NCOL(x$input.matrix), "columns )", "\n")
    print(out.main, ...)

    if (x$me) {
      cat("\n Matrix of measurement errors (rows) and contrasts (columns)  \n")
      print(x$tables$me.contrasts, ...)

      cat("\n Head of the matrix of standard deviations (", NROW(x$sd.input.matrix), "rows * ",
          NCOL(x$sd.input.matrix), "columns )", "\n")
      out.main <- x$sd.input.matrix[1:nrows, , drop = FALSE]
      attr(out.main, 'na.action') <- attr(x$sd.input.matrix, 'na.action')
      print(out.main, ...)

    }
  })
}

#setMethod ("print",
#           signature = "msbm.frame",
#           definition = print.msbm.frame)

#' @exportS3Method base::summary
summary.msbm.frame <- function(object, n, ...) {
  print.msbm.frame (object, n = n, ...)
}
