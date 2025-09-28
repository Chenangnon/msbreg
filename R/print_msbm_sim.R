#'
# @exportS3Method print msbm.sim
#' @exportS3Method base::print
print.msbm.sim <- function (x, digits = getOption("digits") - 3,
                            index = 1L:ncol(x$estimates),
                            ...) {
  cl <- x$call
  estimates <- matrix(x$estimates[, index], nrow = nrow(x$estimates))
  estimates[is.infinite(estimates)] <- NaN

  if ("(Intercept).lambda" %in% colnames(x$estimates)) {
    il <- which("(Intercept).lambda" == colnames(x$estimates))
    Infs <- is.infinite(estimates[,il])
    if (any(Infs)) {
      estimates[Infs, il] <- 36
      il <- which("se.(Intercept).lambda" == colnames(x$estimates))
      zeros <- estimates[,il] == 0 & Infs
      estimates[zeros, il] <- 1e-8
    }
  }


  allNA <- apply(estimates, MARGIN = 2L, FUN = function(t) all(is.na(t)))
  ind1 <- index[allNA]
  index <- index[!allNA]
  estimates <- matrix(estimates[, !allNA], nrow = nrow(estimates))
  rnames <- x$stat.names %||% paste("t", index, "*", sep = "")
  if (length(index) == 0L)
    out.table <- NULL
  else {
    out.table <- cbind(apply(estimates, MARGIN = 2L,
                             FUN = mean, na.rm = TRUE),
                       apply(estimates, MARGIN = 2L,
                             FUN = sd, na.rm = TRUE),
                       apply(estimates, MARGIN = 2L,
                             FUN = min, na.rm = TRUE),
                       apply(estimates, MARGIN = 2L,
                             FUN = median, na.rm = TRUE),
                       apply(estimates, MARGIN = 2L,
                             FUN = max, na.rm = TRUE),
                       apply(estimates, MARGIN = 2L,
                             FUN = IQR, na.rm = TRUE))
    dimnames(out.table) <- list(rnames,
                                c(" Mean", " SD", " Min",
                                  " Median", " Max", " IQR"))
  }

  cat("\nMonte Carlo Simulations \n\n")

  cat("\nCall:\n")
  dput(cl, control = NULL)
  cat("\n\nSummary statistics :\n")
  if (!is.null(out.table))
    print(out.table, digits = digits)
  if (length(ind1) > 0L)
    for (j in ind1) cat(paste("WARNING: All values of ",
                              rnames[j], "* are NA\n", sep = ""))
  invisible(x)
}

#setMethod ("print",
#           signature = "msbm.sim",
#           definition = print.msbm)
