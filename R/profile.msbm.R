#'
#' Method for Profiling \code{msbm} Objects
#'
#' EXPERIMENTAL FUNCTION.
#' Investigates the profile log-likelihood function
#' for a fitted model of class "msbm".
#'
#' @param fitted the original fitted model object.
#'
#' @param which the original model parameters which should be profiled.
#' This can be a numeric or character vector.
#' By default, all parameters are profiled.
#'
#' @param alpha real in \code{(0, 1)}, highest significance level allowed for the profile z-statistics.
#'
#' @param maxsteps integer, maximum number of points to be used for profiling each parameter.
#'
#' @param del positive real, suggested change on the scale of the profile t-statistics.
#' Default value chosen to allow profiling at about 10 parameter values.
#'
#' @param trace logical, should the progress of profiling be reported?
#'
#' @param test character, profile Likelihood Ratio test or Rao Score test.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' The profile z-statistic is defined either as (case \code{test = "LRT"}) the
#' square root of change in deviance with an appropriate sign, or
#' (case \code{test = "Rao"}) as the similarly signed square root of the
#' Rao Score test statistic.
#' The latter is defined as the squared gradient of the profile log likelihood
#' divided by the profile Fisher information, but more conveniently calculated
#' via the deviance of a Gaussian GLM fitted to the residuals of the profiled model.
#'
#' @return
#' A list of classes "profile.msbm" and "profile" with an element for each
#' parameter being profiled. The elements are data-frames with two variables:
#'
#' \item{\code{par.vals}}{ a matrix of parameter values for each fitted model;}
#'
#' \item{\code{tau} or \code{z}}{ the profile t or z-statistics (the name
#' depends on whether there is an estimated dispersion parameter).}
#'
#' @author
#' Chenangnon Tovissode built the function on the model of
#' \link[stats]{profile.glm}.
#'
#' @exportS3Method stats::profile
profile.msbm <- function (fitted, which = 1:npar, alpha = 0.01, maxsteps = 10,
                          del = zmax/5, trace = FALSE, test = c("LRT", "Rao"), ...) {

  stop("not yet implemented")

  test <- match.arg(test)
  Pnames <- names(B0 <- coef(fitted))
  nonA <- !is.na(B0)
  pv0 <- t(as.matrix(B0))
  p <- npar <- length(Pnames)
  if (is.character(which))
    which <- match(which, Pnames)
  summ <- summary(fitted)
  std.err <- summ$coefficients[, "Std. Error", drop = FALSE]
  mf <- model.frame(fitted)
  Y <- model.response(mf)
  n <- NROW(Y)
  O <- model.offset(mf)
  if (!length(O))
    O <- rep(0, n)
  W <- model.weights(mf)
  if (length(W) == 0L)
    W <- rep(1, n)
  OriginalDeviance <- deviance(fitted)
  DispersionParameter <- summ$dispersion
  X <- model.matrix(fitted)
  fam <- family(fitted)
  switch(fam$family, binomial = , poisson = , `Negative Binomial` = {
    zmax <- sqrt(qchisq(1 - alpha, 1))
    profName <- "z"
  }, gaussian = , quasi = , inverse.gaussian = , quasibinomial = ,
  quasipoisson = , {
    zmax <- sqrt(qf(1 - alpha, 1, n - p))
    profName <- "tau"
  })
  prof <- vector("list", length = length(which))
  names(prof) <- Pnames[which]
  for (i in which) {
    if (!nonA[i])
      next
    zi <- 0
    pvi <- pv0
    a <- nonA
    a[i] <- FALSE
    Xi <- X[, a, drop = FALSE]
    pi <- Pnames[i]
    for (sgn in c(-1, 1)) {
      if (trace)
        message("\nParameter: ", pi, " ", c("down", "up")[(sgn +
                                                             1)/2 + 1])
      step <- 0
      z <- 0
      LP <- X[, nonA, drop = FALSE] %*% B0[nonA] + O
      while ((step <- step + 1) < maxsteps && abs(z) <
             zmax) {
        bi <- B0[i] + sgn * step * del * std.err[Pnames[i],
                                                 1]
        o <- O + X[, i] * bi
        fm <- glm.fit(x = Xi, y = Y, weights = W, etastart = LP,
                      offset = o, family = fam, control = fitted$control)
        LP <- Xi %*% fm$coefficients + o
        ri <- pv0
        ri[, names(coef(fm))] <- coef(fm)
        ri[, pi] <- bi
        pvi <- rbind(pvi, ri)
        zz <- (fm$deviance - OriginalDeviance)/DispersionParameter
        if (zz > -0.001)
          zz <- max(zz, 0)
        else stop("profiling has found a better solution, so original fit had not converged")
        if (test == "Rao") {
          r <- fm$residuals
          w <- fm$weights
          fml <- glm.fit(x = X, y = r, weights = w, control = fitted$control,
                         intercept = FALSE)
          zz <- (fml$null.deviance - fml$deviance)/DispersionParameter
          zz <- max(zz, 0)
        }
        z <- sgn * sqrt(zz)
        zi <- c(zi, z)
      }
    }
    si <- order(zi)
    prof[[pi]] <- structure(data.frame(zi[si]), names = profName)
    prof[[pi]]$par.vals <- pvi[si, , drop = FALSE]
  }
  val <- structure(prof, original.fit = fitted, summary = summ)
  class(val) <- c("profile.glm", "profile")
  attr(val, "test") <- test
  val
}
