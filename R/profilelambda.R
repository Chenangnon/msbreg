#' Profile log-likelihood for a MSB model asymptote
#'
#' Computes and optionally plots profile log-likelihoods for the
#' asymptote \eqn{\lambda} of a Multistage Binomial model.
#'
#' @param object an object of class \code{"msbm"} typically, a result of a call
#' to \link[msbreg]{msbreg}.
#'
#' @param lambda vector of values of \eqn{\lambda}. The default value is
#' taken around the estimate of \eqn{\lambda} in \code{object}.
#'
#' @param step small positive real specifying the steps of \eqn{\lambda}.
#' Only used when \code{lambda} is \code{NULL}. Defaults to 0.01.
#'
#' @param eps small positive real, tolerance for \eqn{\lambda = 0} and
#' \eqn{\lambda = 1}. Defaults to 1e-04.
#'
#' @param neg.criterion logical, should the negative of the criterion
#' used to fit \code{object} be considered instead of the log-likelihood?
#' This is relevant only when the MSB model was fitted using a penalized
#' likelihood (the default for \link{msbreg}).
#'
#' @param plotit logical, should the result be plotted? Defaults to \code{FALSE}.
#'
#' @param interp logical, should spline interpolation be used when plotting?
#' Defaults to \code{TRUE} if plotting with \code{lambda} of length less than 200.
#'
#' @param xlab,ylab labels of the axes when \code{plotit = TRUE}.
#' Defaults to \code{'lambda'} and \code{'log-Likelihood'} respectively.
#'
#' @param ci.level level of the confidence interval to be shown on the plot
#' when \code{plotit = TRUE}.
#'
#' @param boundary \code{NULL}, or logical, should the \eqn{\bar{\chi}^2} distribution
#' having a point mass with probability 0.5 at zero be used for the confidence
#' interval? The alternative is the common \eqn{\chi^2} distribution.
#' The default depends on the estimate of \eqn{\lambda} in \code{object} and
#' the \code{lambda} value with the highest profile log-likelihood value.
#'
#' @param ... additional parameters to be passed to or from other methods.
#'
#' @details
#' This is a routine to compute the log-likelihood profile for a constant
#' asymptote \eqn{\lambda} in a MSB model (See \link[msbreg]{msbreg-package}
#' for an introduction to the MSB model). It is possibly time demanding as
#' it repeatedly fits the model for different \eqn{\lambda} values in \code{lambda}.
#'
#' The input argument \code{boundary} should be left as the default \code{NULL},
#' unless you know what you are doing (i.e. how this affects the returned/plotted
#' confidence interval). Using the wrong \code{boundary} value will make the
#' interval too conservative (\code{boundary = FALSE} by error) or too liberal
#' (\code{boundary = TRUE} by error).
#'
#' If \code{plotit = TRUE}, the function plots the log-likelihood function profiled
#' for \eqn{\lambda} and indicates a confidence interval (of level \code{ci.level},
#' defaults to \code{95%}) around the \eqn{\lambda} values corresponding to the maximum
#' observed log-likelihood value. If \code{interp = TRUE}, spline interpolation
#' is used to give a smoother plot.
#'
#' @return
#' A object of class \code{'profile.lambda'} which is a vector (of the same
#' length as the unique values in \code{lambda}) with the attributes:
#'
#' \item{\code{max.index}}{ index of the \eqn{\lambda} value corresponding to
#' the maximum profiled log-likelihood value (or the negative criterion value)
#' in \code{lambda};}
#'
#' \item{\code{pci}}{ profiled confidence interval for \eqn{\lambda};}
#'
#' \item{\code{lambda}}{ the used vector \code{lambda};}
#'
#' \item{\code{boundary}}{ the used \code{boundary} value (logical);}
#'
#' \item{\code{ML}}{ the log-likelihood value corresponding to an
#' unrestricted estimate of the MSB model in \code{object};}
#'
#' \item{\code{L1}}{ the log-likelihood value corresponding to \eqn{\lambda = 1}.}
#'
#' When \code{plotit = TRUE}, the output is invisibly returned (not automatically
#' printed as usually when an output is not assigned).
#'
#' The \code{'profile.lambda'} class has a \link[base]{plot} method
#' for visualizing the calculated log-likelihood profile.
#'
#' @aliases profilelambda
#' @export profilelambda
#'
#' @importFrom stats spline
#' @importFrom graphics segments
#' @importFrom graphics abline
#' @importFrom grDevices dev.hold
#' @importFrom graphics text
#' @importFrom stats qchisq
#'
#' @examples
#' ## Two-stage logistic model
#' data(test1data)
#' attr(test1data$y, "formula")
#'
#' # Fit and summarize the model
#' MSBres <- msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                    lambda.formula = ~ 1, weights = Total,
#'                    data = test1data)
#'
#' summary (MSBres)
#'
#' # Compute the profile log-likelihood and
#' # plot with the default 95% confidence interval
#' PL <- profilelambda(MSBres, plotit = TRUE,
#'                     ylab = 'Penalized log-likelihood')
#'
#' # Show a 90% confidence interval
#' plot(PL, ci.level = .9,
#'      ylab = 'Penalized log-likelihood')
#'
profilelambda <- function(object, lambda = NULL, step = 0.01, eps = 1e-04,
                          neg.criterion = TRUE, plotit = FALSE, interp = NULL,
                          xlab = expression(lambda), ylab = "Log-Likelihood",
                          ci.level = 0.95, boundary = NULL, ...) {

  # object must be of class 'msbm'
  stopifnot(inherits(object, 'msbm'))
  stopifnot(eps > 0)
  stopifnot(eps <= 0.05)
  stopifnot(ci.level > 0.5, ci.level < 1)
  if(!is.null(boundary)) {
    boundary <- boundary[1]
    stopifnot(is.logical(boundary))
  }

  # Include model frame
  object$frame <- model.frame(object)

  # MAP estimate
  pd <- object$frame$dim$p+object$frame$dim$d
  object$call$start <- object$coefficients[1:pd]
  Nullfit <- update(object, lambda.formula. = ~ 0, evaluate = TRUE)
  if (object$frame$dim$r > 0) {
    object$call$start <- object$coefficients[1:(pd+1)]
    Fullfit <- update(object, lambda.formula. = ~ 1, evaluate = TRUE)
  }
  else {
    object$call$start <- c(object$coefficients[1:pd], object$link$linkfun(min(object$lambda.values[1], .99)))
    Fullfit <- update(object, lambda.formula. = ~ 1, evaluate = TRUE)
  }

  # Lambda values
  step <- max(1e-8, min(step[1], 0.1))
  mapL <- Fullfit$lambda.values[1]
  if (is.null(lambda)) {
    if (mapL <= 0.05) {
      lambda <- seq(from = eps, to = 0.1, by = step)
    }
    else if (mapL <= 0.1) {
      lambda <- seq(from = 0.01, to = 0.15, by = step)
    }
    else if (mapL <= 0.15) {
      lambda <- seq(from = 0.05, to = 0.2, by = step)
    }
    else if (mapL <= 0.2) {
      lambda <- seq(from = 0.1, to = 0.3, by = step)
    }
    else if (mapL >= 0.95) {
      lambda <- c(seq(from = 0.9, to = 1, by = step), 1)
    }
    else if (mapL > 0.9) {
      lambda <- seq(from = 0.85, to = 0.99, by = step)
    }
    else if (mapL > 0.8) {
      lambda <- seq(from = 0.75, to = 0.95, by = step)
    }
    else {
      lambda <- seq(from = mapL - 0.05, to = mapL + 0.05, by = step)
    }
  }
  else {
    lambda <- lambda[lambda > eps & lambda <= 1]
  }
  lambda <- as.vector(sort(unique(lambda)))
  m <- length(lambda)
  stopifnot(m > 0)

  # Function to extract profile likelihood
  extractLL <- if(neg.criterion)
    function(object) - object$criterion[1]
  else
    function(object) logLik(object)[1]

  # Isolate the case lambda = 1
  if (1 %in% lambda)
    m_0 <- m - 1
  else
    m_0 <- m

  # Loop over lambda values
  loglik <- numeric(m)
  for (k in 1L:m_0) {
    object$call$lambda.formula <- as.formula(paste0("~ offset(qlogis(", lambda[k], ")) + 0"))
    object$call$start <- object$coefficients[1:pd]
    MAPfit_lv <- eval(object$call)

    loglik[k] <- extractLL(MAPfit_lv)
  }

  if(m_0 < m) {
    loglik[m] <- extractLL(Nullfit)
  }

  ml <- which.max(loglik)[1L]
  lambda.maxLik <- lambda[ml]
  if(is.null(boundary)) {
    boundary <- as.vector((lambda.maxLik >= 1 - eps) | (mapL >= 1 - eps))
  }

  Lmax <- loglik[ml]
  if (boundary) {
    lim <- Lmax - qchisq(2 * ci.level - 1, df = 1, lower.tail = TRUE)/2
  }
  else {
    lim <- Lmax - qchisq(ci.level, df = 1, lower.tail = TRUE)/2
  }
  pci <- lambda[range(which(loglik > lim))]
  attr(pci, 'ci.level') <- ci.level

  out <- structure(loglik,
                   lambda = lambda,
                   max.index = ml,
                   pci = pci,
                   boundary = boundary,
                   ML = extractLL(Fullfit)[1],
                   L1 = extractLL(Nullfit),
                   class = 'profile.lambda')

  if (plotit) {

    plot.profile.lambda (out, interp = interp, xlab = xlab, ylab = ylab,
                         ci.level = ci.level, ...)
  }
  else {

    return(out)
  }
}

#' @exportS3Method base::plot
plot.profile.lambda <- function (x, interp = NULL, ci.level = 0.95,
                                 xlab = expression(lambda),
                                 ylab = "Log-Likelihood", ...) {
  out <- x; rm(x)
  stopifnot(ci.level > 0.5, ci.level < 1)

  m <- length(out)
  if (is.null(interp))
    interp <- m < 200
  else
    interp <- interp[1] == 1

  lambda <- attr(out, 'lambda')
  loglik <- as.vector(out)

  if (interp) {
    sp <- stats::spline(lambda, loglik, n = 200)
    lambda <- sp$x
    loglik <- sp$y
    m <- length(lambda)
  }

  # Plotting code adapted from MASS::boxcox
  mx <- which.max(loglik)[1L]
  Lmax <- loglik[mx]
  if (attr(out, 'boundary')) {
    lim <- Lmax - qchisq(2 * ci.level - 1, df = 1, lower.tail = TRUE)/2
  }
  else {
    lim <- Lmax - qchisq(ci.level, df = 1, lower.tail = TRUE)/2
  }
  dev.hold()
  on.exit(dev.flush())
  plot(lambda, loglik, xlab = xlab, ylab = ylab, type = "l",
       ylim = range(loglik, lim))
  plims <- par("usr")
  abline(h = lim, lty = 3)
  y0 <- plims[3L]
  scal <- (1/10 * (plims[4L] - y0))/par("pin")[2L]
  scx <- (1/10 * (plims[2L] - plims[1L]))/par("pin")[1L]
  text(lambda[1L] + scx, lim + scal, paste0(" ", 100 * ci.level, "%"), xpd = TRUE)
  la <- lambda[mx]
  if (mx > 1 && mx < m)
    segments(la, y0, la, Lmax, lty = 3)
  ind <- range((1L:m)[loglik > lim])

  if (loglik[1L] < lim) {
    i <- ind[1L]
    x <- lambda[i - 1] + ((lim - loglik[i - 1]) * (lambda[i] -  lambda[i - 1]))/(loglik[i] - loglik[i - 1])
    segments(x, y0, x, lim, lty = 3)
  }

  if (loglik[m] < lim) {
    i <- ind[2L] + 1
    x <- lambda[i - 1] + ((lim - loglik[i - 1]) * (lambda[i] - lambda[i - 1]))/(loglik[i] - loglik[i - 1])
    segments(x, y0, x, lim, lty = 3)
  }

  invisible(out)
}
