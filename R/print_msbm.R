#' Summarizing a Multistage Binomial Model Fits
#'
#' Print a quick summary of a multistage binomial model fit.
#' These are \link{methods} for classes \code{"msbm"} and
#' \code{"summary.msbm"}.
#'
# @param object an object of class "msbm",
# usually, the result of a call to \link[msbreg]{msbreg}.
#
#' @param x an object of class \code{"msbm"} or \code{"summary.msbm"},
#' typically, a result of a call to \link[msbreg]{msbreg} or
#' \link[msbreg]{summary.msbm}.
#'
#' @param digits,symbolic.cor,signif.stars,show.residuals Same
#' as for \link[stats]{summary.glm}:
#' \itemize{
#' \item \code{digits} the number of significant digits to use when printing;
#'
#' \item \code{symbolic.cor} logical, if \code{TRUE}, print the correlations
#' in a symbolic form (see \link[stats]{symnum}) rather than as numbers;
#'
#' \item \code{signif.stars} logical, if \code{TRUE}, *significance stars*
#' are printed for each coefficient;
#'
#' \item \code{show.residuals} logical, if \code{TRUE} then a summary of the
#' deviance residuals is printed at the head of the output.;
#' }
#'
#' @param me logical, if \code{TRUE}, the measurement error component
#' of the model is summarized and printed.
#'
#' @param frame logical, if \code{TRUE}, the model frame (component
#' \code{$frame} of the input) is returned/printed.
#'
#' @param ... further arguments passed to or from other methods.
#
# @usage
# ## S3 method for class 'msbm'
# print (x, digits = max(3L, getOption("digits") - 3L),
#        me = FALSE, frame = FALSE, ...)
#
# ## S3 method for class 'msbm'
# summary (object, fisher.matrix = TRUE, se.values = FALSE,
#          IC = NULL, correlation = FALSE, symbolic.cor = FALSE,
#          rsquared.method = 'KL', cor.method = "pearson", frame = FALSE, ...)
#
# ## S3 method for class 'summary.msbm'
# print (x, digits = max(3L, getOption("digits") - 3L),
#        symbolic.cor = x$symbolic.cor,
#        signif.stars = getOption("show.signif.stars"),
#        show.residuals = FALSE, frame = FALSE,
#        me = FALSE, ...)
#
#' @details
#' These functions mimic the similar \code{print} methods for
#' objects of classes "\link[stats]{glm}" and "\link[stats]{summary.glm}".
#'
#' @aliases print.summary.msbm
#'
#' @return
#' The \code{print} methods invisibly returns \code{x}, after printing the main
#' features of \code{x}.
#'
# @exportS3Method print msbm
# @export print.msbm
#'
#' @examples
#' # Example 1
#' set.seed(10)
#' data("test1data")
#'
#' MSBfit = msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                  alpha.formula = ~ 0,
#'                  lambda.formula = ~ 0,
#'                  weights = Total,
#'                  data = test1data)
#' MSBfit # Print coefficient estimates (descriptive)
#'
#' summary (MSBfit) # Print a more detailed summary of the fit (inference)

#' @exportS3Method base::print
print.msbm <- function (x,
                        digits = max(3L, getOption("digits") - 3L),
                        me = FALSE,
                        frame = FALSE,
                        ...) {
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients")
    cat(":\n")
    print.default(format(x$coefficients, digits = digits),
                  print.gap = 2, quote = FALSE)
  }
  else cat("No coefficients\n\n")

  alpha.lambda.values <- x$alpha.values * x$lambda.values
  if (length(unique(alpha.lambda.values)) == 1 & length(unique(x$lambda.values)) == 1) {
    cat("\nSuccess probability bounds: [")
    cat(format(alpha.lambda.values[1], digits = digits))
    cat(", ")
    cat(format(x$lambda.values[1], digits = digits))
    cat("]")
  }
  else {
    cat("\nMinimum success probability: ")
    if (length(unique(alpha.lambda.values)) == 1) {
      cat(format(c(Minimum = alpha.lambda.values[1]), digits = digits))
    }
    else {
      cat("\n ")
      print.default(format(summary(alpha.lambda.values),
                           digits = digits),
                    print.gap = 2, quote = FALSE)
    }
    cat("\nMaximum success probability: ")
    if (length(unique(x$lambda.values)) == 1)
      cat(format(c(Maximum = x$lambda.values[1]), digits = digits))
    else {
      cat("\n ")
      print.default(format(summary(x$lambda.values), digits = digits),
                    print.gap = 2, quote = FALSE)
    }
  }
  cat("\n ")

  if (me & (x$me | x$me.offset)) {
    cat("\nStandard deviation(s) of measurement error(s): ")
    cat("\n ")
    print.default(format(summary(x$me.sd), digits = digits),
                  print.gap = 2, quote = FALSE)
  }

  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
      x$df.residual, "Residual\n")
  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("Null Deviance:\t   ",
      format(signif(x$null.deviance, digits)),
      "\nResidual Deviance: ",
      format(signif(x$deviance, digits)),
      "\tAIC:", format(signif(x$aic, digits)))
  cat("\n")

  if(frame & !is.null(x$frame)) {
    x0 <- x
    x <- x$frame
    n <- 6L
    eval(print_frame_msbm())
    invisible(x0)
  }
  else
    invisible(x)
}
# print_frame_msbm
#setMethod ("print",
#           signature = "msbm",
#           definition = print.msbm)

# @exportS3Method print summary.msbm
#' @rdname print.msbm
#' @exportS3Method base::print
print.summary.msbm <- function (x, digits = max(3L, getOption("digits") - 3L),
                                symbolic.cor = x$symbolic.cor,
                                signif.stars = getOption("show.signif.stars"),
                                show.residuals = FALSE,
                                frame = FALSE,
                                me = FALSE, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")

  if (show.residuals) {
    cat("Deviance Residuals: \n")
    x$deviance.resid <- x$deviance.resid * sign(x$y.resid)
    if (x$df.residual > 5) {
      x$deviance.resid <- setNames(quantile(x$deviance.resid,
                                            na.rm = TRUE),
                                   c("Min", "1Q", "Median", "3Q", "Max"))
    }
    xx <- zapsmall(x$deviance.resid, digits + 1L)
    print.default(xx, digits = digits, na.print = "", print.gap = 2L)
  }

  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
     overparametrized <- x$overparametrized %||% FALSE

    df <- if ("df" %in% names(x))
      x[["df"]]
    else NULL
    if (any(is.na(df))) {
      cat("\nCoefficients:\n")
    }
    else if (!is.null(df) && (nsingular <- df[3L] - df[1L])) {
      if (overparametrized) {
        cat("\nCoefficients: (", nsingular, " with zero standard error(s) because of overparametrization)\n",
            sep = "")
      }
      else {
        cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
          sep = "")
      }
    }
    else cat("\nCoefficients:\n")
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn,
                                                               colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)

    cat("\n(Dispersion parameter for the binomial family taken to be ",
        format(x$dispersion), ")\n")
  }

  if (length(unique(x$alpha.lambda)) == 1 & length(unique(x$lambda)) == 1) {
    cat("\nSuccess probability bounds: [")
    cat(format(x$alpha.lambda[1], digits = digits))
    cat(", ")
    cat(format(x$lambda[1], digits = digits))
    cat("]")
  }
  else {
    cat("\nMinimum success probability: ")
    if (length(unique(x$alpha.lambda)) == 1) {
      cat(format(c(Minimum = x$alpha.lambda[1]), digits = digits))
    }
    else {
      cat("\n ")
      print.default(format(summary(x$alpha.lambda),
                           digits = digits),
                    print.gap = 2, quote = FALSE)
    }
    cat("\nMaximum success probability: ")
    if (length(unique(x$lambda)) == 1)
      cat(format(c(Maximum = x$lambda[1]), digits = digits))
    else {
      cat("\n ")
      print.default(format(summary(x$lambda), digits = digits),
                    print.gap = 2, quote = FALSE)
    }
  }

  if (me & (x$me | x$me.offset)) {
    cat("\n\nStandard deviation(s) of measurement error(s): ")
    cat("\n ")
    print.default(format(summary(x$me.sd), digits = digits),
                  print.gap = 2, quote = FALSE)
  }

  cat("\n\n")

  cat(apply(cbind(paste(format(c(" Null", " Residual"),
                               justify = "right"), "deviance:"),
                  format(unlist(x[c("null.deviance", "deviance")]),
                         digits = max(5L, digits + 1L)),
                  " on ", format(unlist(x[c("df.null", "df.residual")])),
                  " degrees of freedom\n"),
            1L, paste, collapse = " "), sep = "")

  if (!anyNA(unlist(x[c("null.deviance", "deviance")]))) {
    LRstat <- diff(unlist(x[c("deviance", "null.deviance")]))
    LRdf <- diff(unlist(x[c("df.residual", "df.null")]))
    cat("      LR-statistic:", formatC(LRstat, digits = digits),
        " on ", LRdf, " DF,  p-value:",
        format.pval(pchisq(LRstat, df = LRdf, lower.tail = FALSE), digits = digits))
    cat("\n")
  }

  if (!is.null(x$r.squared)) {
    r.squared <- x$r.squared[1]
    adj.r.squared <- 1 - (1 - r.squared) * x[["df.null"]] / x[["df.residual"]]

    cat("Multiple R-squared:", formatC(r.squared, digits = digits))
    cat(",  Adjusted R-squared:", formatC(adj.r.squared, digits = digits))
    cat("\n")

  }

  if (nzchar(mess <- naprint(x$na.action)))
    cat("  (", mess, ")\n", sep = "")
  cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),
      "\n\n", sep = "")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      }
      else {
        correl <- format(round(correl, 2L), nsmall = 2L,
                         digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
      cat("\n")
    }
  }

  if (identical(x$criterion, "ML")) {
    cat("Fit: Maximum Likelihood \n")
  }
  else if (identical(x$criterion, "MLJ")) {
    cat("Fit: Maximum A Posteriori (with Jeffreys non-informative prior) \n")
  }
  else if (identical(x$criterion, "MLPJ")) {
    cat("Fit: Maximum A Posteriori (with Jeffreys prior for slopes) \n")
  }
  else if (identical(x$criterion, "MLJIC")) {
    cat("Fit: Maximum A Posteriori (with Jeffrey prior and intercept-correction) \n")
  }
  cat(paste0("(Converged?: ", x$converged, ")\n"))

  if(frame & !is.null(x$frame)) {
    x0 <- x
    x <- x$frame
    n <- 6L
    eval(print_frame_msbm())

    invisible(x0)
  }
  else
    invisible(x)
}

#setMethod ("print",
#           signature = "summary.msbm",
#           definition = print.summary.msbm)
