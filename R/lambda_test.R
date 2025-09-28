#' Likelihood Ratio Test
#'
#' Likelihood ratio test for \eqn{\lambda = \lambda_0} in multistage
#' binomial models.
#'
#' @param object description
#'
#' @param lambda scalar in \code{(0, 1]}.
#'
#' @param alternative When \code{lambda = 1}, the function always uses
#' \code{alternative = 'less'}.
#'
#' @param simulate.p.value description
#'
#' @param B description
#'
#' @export lambda.test
#'
lambda.test <- function (object, lambda = 1, # Only H0: lambda = 1 is implemented
                         alternative = "less", # Not implemented
                         simulate.p.value = FALSE,  # Not implemented
                         B = 2000) {
  if(!inherits(object, "msbm") & inherits(object, "glm")) {
    object <- glm.to.msbm(object)
  }

  mf <- object$frame <- model.frame(object)
  r <- mf$dims$r
  intercept <- attr(mf$lambda.dictionary, "intercept")
  method <- object$method
  link <- object$link
  control <- object$control

  if (lambda == 1) {
    if (r == 0) {
      attr(mf$lambda.dictionary, "intercept") <- mf$dims$r <- 1
      mf$parnames <- c(mf$parnames, "(Intercept).lambda")
      object1 <- eval(call(if (is.function(method)) "method" else method,
                           frame = mf, y = mf$y, start = c(object$coefficients, 0),
                           link = link, control = control))

      Dev0 <- object$deviance
      Dev1 <- object1$deviance
      df <- 1
      H0 <- "H0: lambda = 1"
      H1 <- "H1: 0 < lambda < 1"
    }
    else if (r == 1) {
      if (intercept) {
        attr(mf$lambda.dictionary, "intercept") <- mf$dims$r <- 0
        H0 <- "H0: lambda = 1"
        H1 <- "H1: 0 < lambda < 1"
      }
      else {
        mf$dims$r <- 0
        attrL <- attributes(mf$lambda.dictionary)
        attrL$names <- NULL
        mf$lambda.dictionary <- numeric(0L)
        attributes(mf$lambda.dictionary) <- attrL
        H0 <- "H0: lambda_i = 1"
        H1 <- "H1: 0 < lambda_i < 1"
      }

      mf$dims$npars <- mf$dims$npars - 1
      mf$parnames <- mf$parnames[-length(mf$parnames)]
      mf$parindex$intercepts <- mf$parindex$intercepts[-length(mf$parindex$intercepts)]
      if (object$IC) {
        mf$parindex$eta.intercepts <- mf$parindex$eta.intercepts[-length(mf$parindex$eta.intercepts)]
      }

      object0 <- eval(call(if (is.function(method)) "method" else method,
                           frame = mf, y = mf$y, start = object$coefficients[-length(object$coefficients)],
                           link = link, control = control))

      Dev0 <- object0$deviance
      Dev1 <- object$deviance
      df <- 1
    }
    else {
      attr(mf$lambda.dictionary, "intercept") <- mf$dims$r <- 0
      H0 <- "H0: lambda_i = 1"
      H1 <- "H1: 0 < lambda_i < 1"

      attrL <- attributes(mf$lambda.dictionary)
      attrL$names <- NULL
      mf$lambda.dictionary <- numeric(0L)
      attributes(mf$lambda.dictionary) <- attrL

      object0 <- eval(call(if (is.function(method)) "method" else method,
                           frame = mf, y = mf$y, start = object$coefficients[1:sum(unlist(object$frame$dims[c("p", "d")]))],
                           link = link, control = control))


      Dev0 <- object0$deviance
      Dev1 <- object$deviance
      df <- r
    }
  }
  else {
    stopifnot(lambda > 0, lambda < 1)


  }

  DNAME <- paste0(c(deparse(object$call), "\n"), collapse = "")
  METHOD <- "Likelihood ratio (LR) test"
  STATISTIC <- max(Dev0 - Dev1, 0)
  PARAMETER <- c(df, object$dims$nobs)
  names(STATISTIC) <- "LR"
  names(PARAMETER) <- c("df", "nobs")

  if (PARAMETER[1] == 1)
    PVAL <- 0.5 * pchisq(STATISTIC, PARAMETER[1], lower.tail = FALSE)
  else {
    # Could be different from chibar(0, df)!
    # Check the number of mixed chisq distributions given df, and their weights
    PVAL <- 0.5 * pchisq(STATISTIC, PARAMETER[1], lower.tail = FALSE)
  }

  structure(list(statistic = STATISTIC,
                 parameter = PARAMETER,
                 H0 = H0,
                 H1 = H1,
                 p.value = PVAL,
                 method = METHOD,
                 data.name = DNAME),
            class = "htest")
}
