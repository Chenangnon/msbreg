# https://bookdown.dongzhuoer.com/hadley/r-pkgs/man#multiple-man

#' Draw Model Plots
#'
#' Draws the curve or perspective plot of a function/model over a region.
#' \code{courbe} is a generic function.
#'
#' @param object an \code{R} object of a class with a \code{courbe} method.
#' The default method works on a function of one variable (of class
#' "\link{function}") to be drawn, or a numeric vector (of class
#' "\link{numeric}") of length from 2 to 5 (representing a univariate
#' or a bivariate multistage binomial (MSB) model).
#'
#' @param ... additional argument to \code{courbe}.
#' The default method for an \code{object} of class "\link{function}"
#' takes any argument of \link[graphics]{curve} or \link[graphics]{persp}.
#' For object of class \code{"numeric"}, \code{linkinv} and \code{alpha} (
#' \eqn{\alpha} model component) are used.
#'
#' @param from,to,n,add,type,xname,yname,xlab,ylab,zlab,log,xlim,ylim,zlim Plotting
#' arguments. See \link[graphics]{curve}.
#'
#' @param jitter logical, should the plot be jittered?
#' Defaults to \code{TRUE}. This is only used when object is of class
#' \code{msbm} and has only one covariate (2-D plot).
#'
#' @param jitter.range optional numeric scalar, range of the jitter/random
#' noise (half \code{jitter.range} below and half \code{jitter.range} above).
#' Only used for 2-D plot when \code{jitter = TRUE}.
#'
#' @param seed an integer to seed the random generator of \code{R} before using
#' the function \link[stats]{runif} to add noise to the binary response.
#' Only used for 2-D plot when \code{jitter = TRUE}.
#'
# @param linkinv function, inverse link function of a MSB model.
#'
# @param alpha scalar in \eqn{[0, 1)}, minimum success probability
# of the response variable in a MSB model.
#'
#' @param main,sub,theta,phi,expand,ltheta,shade,r,d,scale,border,box,axes,nticks,ticktype Plotting
#' arguments. See \link[graphics]{persp}.
#'
#' @param col,lty,lwd,cex plotting arguments. See \link[graphics]{curve} or
#' \link[graphics]{persp}.
#'
# @usage
# courbe (object, ...)
#
# # Default S3 method
# courbe (object, from = NULL, to = NULL, n = 101,
#         add = FALSE, type = "l", xname = "x", yname = "y",
#         xlab = "x", ylab = NULL, zlab = NULL, log = NULL,
#         xlim = NULL, ylim = NULL, zlim = NULL,
#         linkinv = stats::plogis, alpha = 0, main = " ",
#         col = "black", lty = 1, lwd = 1, cex = 0.1,
#         theta = 50, phi = 10, expand = 0.8,
#         ltheta = 0, shade = 0.75, r = sqrt(3),
#         d = 1, scale = TRUE, border = NULL,
#         box = TRUE, axes = TRUE, nticks = 5,
#         ticktype = "detailed", ...)
#
#'
#' @details
#' The default \code{courbe} method is an extension
#' of the function \link[graphics]{curve} to handle low dimensional Multistage
#' Binomial (MSB) model curves. The default method indeed calls
#' \link[graphics]{curve} when \code{object} is a \code{"function"}
#' Note that, unlike \link{curve}, a first argument of class
#' \code{"expression"} or \code{"call"} to \code{courbe} will generally
#' result into an error, quoting the \code{"expression"} or \code{"call"}
#' using \link{quote} will however work.
#'
#' For a \code{"numeric"} class \code{object} that can be interpreted
#' as the parameters (\code{coefficients}) of a univariate (the length of
#' \code{object} is two or three) or a bivariate (the length of
#' \code{object} is four or five) MSB model, the default method
#' creates a function returning the response success probability
#' of the MSB model. For a univariate MSB model, the success probability
#' function is passed to \link[graphics]{curve}. For a bivariate MSB model,
#' the success probability function is passed to \link[graphics]{persp}.
#'
#' *MSB model from an object of class \code{"numeric"}*
#'
#' When \code{object} is a numeric vector of length two or three,
#' the first two elements of \code{object} are respectively the intercept
#' and the slope of the univariate (one stage) MSB model. If a third
#' element is present in \code{object}, it gives the maximum success
#' probability after applying an inverse link function (default is the
#' \link[msbreg]{logit} link). If otherwise absent, the maximum success
#' probability is set to one. Note that the third element of \code{object}
#' can be \link{Inf}, which corresponds to a maximum success probability
#' at one.
#'
#' When \code{object} is a numeric vector of length four or five,
#' a function returning the response success probability in a bivariate
#' MSB model is created and passed to \link[graphics]{persp}.
#' If the length of \code{object} is four, the fifth element is
#' taken to be \link{Inf} which corresponds to a maximum success
#' probability of one.
#'
#' For the default method, the arguments \code{...} can be any formal
#' argument accepted by \link[graphics]{curve} or \link[graphics]{persp}.
#'
#' The method for class \code{"msbm"} also works for MSB models with only
#' one or two predictors (all numeric). The difference with passing
#' a numeric vector to the default method is essentially the inclusion
#' of the observed data into the plot.
#'
#' Methods can be added for special objects.
#' The function \code{curves} is a simple alias of \code{courbe},
#' and methods should be written for \code{courbe} (and not for \code{curves}).
#'
#' @return Invisibly returns the plotted object.
#'
#' For the default method:
#'
#' If \code{object} is a \code{"function"} or a quoted \code{"expression"}
#' or \code{"call"}, the plotted object is a list with components \code{x}
#' and \code{y} of the points that were drawn, returned invisibly.
#'
#' If \code{object} is a numeric vector of length two or three,
#' the plotted object is also a list with components \code{x}
#' and \code{y} of the points that were drawn, returned invisibly.
#' The same format of output is returned for a univariate \code{"msbm"}
#' class object.
#'
#' If \code{object} is a numeric vector of length four or five, \code{courbe}
#' invisibly returns a \eqn{4 \times 4} *viewing transformation matrix*
#' as returned by \link[graphics]{persp}.
#' The same format of output is returned for a bivariate \code{"msbm"}
#' class object.
#'
#' @importFrom graphics curve
#' @importFrom graphics persp
#' @importFrom utils tail
#' @export courbe
#' @export curves
# @export courbe.default
#'
#' @aliases courbe
#' @aliases curves
#' @aliases courbe.default
#'
#' @examples
#' ##** Equivalent to 'curve'
#' courbe (sin, -2*pi, 2*pi, xname = "t")
#'
#' ##** Univariate MSB model with a maximum success probability at 1
#' courbe (c(3, -2), -1, 5, xname = "x")
#'
#' ##** Univariate MSB model with a maximum success probability at 0.75
#' courbe (c(3, -2, qlogis(0.75)), -1, 5, xname = "x", ylim = c(0, 1))
#'
#' ##** Two stage binomial model with a maximum success probability at 0.9
#' courbe (c(3, -2, 3, -2, qlogis(0.9)), -1, 5,)
#'
#' ##** MSB fit on the Infertility data
#' data ("infert", package = "datasets")
#'
#' ## Logistic regression fit to the infert data
#' GLMres <- glm (case ~ spontaneous,
#'                data = infert,
#'                family = binomial())
#'
#' GLMres
#'
#' ## MSB model fit with a maximum success probability in [0, 1)
#' require(msbreg)
#' MSBres <- msbreg (case ~ spontaneous,
#'                   lambda.formula = ~ 1,
#'                   data = infert)
#'
#' MSBres
#'
#' #* Plot the two model fits with the observed data
#' courbe(MSBres, col = 'blue', xlim = c(-2, 4))
#' Tb <- table(infert$case, infert$spontaneous)
#' points(0:2, Tb[2,]/colSums(Tb))
#' courbe(glm.to.msbm(GLMres), xlim = c(-2, 4), col = 'red', add = TRUE)
#' legend(x = -2, y = 0.9, legend = c("data", "GLM", "MSB"),
#'        col = c("black", "red", "blue"),
#'        pch = c("o", ".", "."), lty = c(0, 1, 1))
#'
#' ## Bivariate MSB model
#' MSBres <- msbreg (case ~ spontaneous | age,
#'                   lambda.formula = ~ 1,
#'                   data = infert)
#'
#' courbe(MSBres, theta = 0, col = "lightblue", phi = 30,
#'        expand = 1, ltheta = 120)
#'
courbe <- function(object, ...) {
  UseMethod("courbe")
}

courbe.default <- function (object, ...) {

  if (inherits(object, "msbm") | inherits(object, "glm2msbm")) {
    return(courbe.msmb(object, ...))
  }

  courbe_default_core (object, ...,
                       sexpr = substitute(object),
                       envCall = parent.frame())

}

setGeneric(name = "courbe",
           def = courbe.default)

curves <- courbe

courbe_default_core <- function(object,
                                from = NULL, to = NULL, n = 101,
                                add = FALSE, type = "l", xname = "x", yname = "y",
                                xlab = NULL, ylab = NULL, zlab = NULL, log = NULL,
                                xlim = NULL, ylim = NULL, zlim = NULL,
                                linkinv = stats::plogis, alpha = 0, main = " ",
                                col = NULL, lty = 1, lwd = 1, cex = 0.1,
                                theta = 50, phi = 10, expand = 0.8,
                                ltheta = 0, shade = 0.75, r = sqrt(3),
                                d = 1, scale = TRUE, border = NULL,
                                box = TRUE, axes = TRUE, nticks = 5,
                                ticktype = "detailed", ...,
                                sexpr = substitute(object),
                                envCall = parent.frame(n=2)) {

  if (!is.numeric(object)) {
    defaultcurve <- FALSE
    if (is.name(sexpr)) {
      func <- object

      if (is.null(ylab)) {
        ylab <- paste0(deparse(sexpr), "(", xname, ")")
      }

      defaultcurve <- TRUE
    }
    else {
      if (is.call (sexpr)) {
        func <- function(x) {
          lenv <- envCall
          lenv$x <- x
          assign(xname, value = x, envir = lenv)
          out <- catch.conditions({
            eval(sexpr, envir = lenv)
          })$value

          if (any(class(out) %in% c("simpleError", "error",
                                    "condition", "try-error"))) {
            out <- eval(object, envir = lenv)
          }

          if (is.call(out)) {
            out <- eval(out, envir = lenv)
            attr(out, 'deparse') <- TRUE
          }
          else {
            attr(out, 'deparse') <- FALSE
          }

          out
        }

        if (is.null(ylab)) {
          if (attr(func((from+to)/2), 'deparse')) {
            ylab <- deparse(eval(sexpr, envir = envCall))
          }
          else {
            ylab <- deparse(sexpr)
          }
        }

        defaultcurve <- TRUE
      }
    }

    if (defaultcurve) {
      if (is.null(xlab))
        xlab <- xname
      if (is.null(col))
        col <- "black"
      return(graphics::curve(expr = func, from = from, to = to,
                             n = n, add = add, type = type, xname = xname,
                             xlab = xlab, ylab = ylab, log = log, xlim = xlim,
                             col = col, lty = lty, lwd = lwd, cex = cex, ...))
    }
  }
  courbe_default_add (object,
                      from = from, to = to, n = n,
                      add = add, type = type, xname = xname, yname = yname,
                      xlab = xlab, ylab = ylab, zlab = zlab, log = log,
                      xlim = xlim, ylim = ylim, zlim = zlim,
                      linkinv = linkinv, alpha = alpha, main = main,
                      col = col, lty = lty, lwd = lwd, cex = cex,
                      theta = theta, phi = phi, expand = expand,
                      ltheta = ltheta, shade = shade, r = r,
                      d = d, scale = scale, border = border,
                      box = box, axes = axes, nticks = nticks,
                      ticktype = ticktype, ...)
}

courbe_default_add <- function(object, from = NULL, to = NULL, n = 101,
                               add = FALSE, type = "l", xname = "x", yname = "y",
                               xlab = NULL, ylab = NULL, zlab = NULL, log = NULL,
                               xlim = NULL, ylim = NULL, zlim = NULL,
                               linkinv = stats::plogis, alpha = 0, main = " ",
                               col = NULL, lty = 1, lwd = 1, cex = 0.1,
                               theta = 50, phi = 10, expand = 0.8,
                               ltheta = 0, shade = 0.75, r = sqrt(3),
                               d = 1, scale = TRUE, border = NULL,
                               box = TRUE, axes = TRUE, nticks = 5,
                               ticktype = "detailed", ...) {
  stopifnot(is.numeric(object))
  beta <- object
  p <- length(beta)
  stopifnot(p > 0)
  if (p == 1) {
    beta <- c(0, beta)
    p <- 2
  }
  else if (p > 5) {
    warning(paste0("In courbe (object, ...): numeric 'object' of length ", p,
                   ": only the first five elements of 'object' will be used."),
            call. = FALSE)
    beta <- beta[1:5]
    p <- 5
  }

  if (is.null(from))
    from <- -3
  if (is.null(to))
    to <- 3
  if (is.null(xlim))
    xlim <- c(from, to)
  else {
    xlim <- rep(xlim, length.out = 2)
  }
  if (is.null(xlab))
    xlab <- xname

  # lambda
  if (p %in% c(2, 4)) {
    lambda <- 1
    p <- p / 2
  }
  else {
    lambda <- linkinv (beta[p])
    beta <- beta[-p]
    p <- (p - 1) / 2
  }
  if (is.null(col))
    col <- if (p == 1) "black" else "lightblue"

  # Two dim plot if one predictor
  if (p == 1) {
    if (is.null(ylab))
      ylab <- paste0("Prob(Y = 1|", xname, ")")

    xfun <- function(x) {
      lambda * linkinv (beta[1] + beta[2] * x)
    }
    res <- graphics::curve(xfun,
                           from = xlim[1], to = xlim[2], ylim = ylim,
                           col = col, lty = lty, lwd = lwd, cex = cex,
                           xname = xname, xlab = xlab, ylab = ylab,
                           add = add, ...)

    return(invisible(res))
  }

  stopifnot(alpha >= 0, alpha < 1)
  if (is.null(ylim))
    ylim <- xlim
  else {
    ylim <- rep(ylim, length.out = 2)
  }
  if (is.null(zlim))
    zlim <- c(0, 1)
  else {
    zlim <- rep(zlim, length.out = 2)
  }

  if (is.null(ylab))
    ylab <- yname

  if (is.null(zlab))
    zlab <- paste0("Prob(Y = 1|", xname, ",", yname, ")")

  xvalues <- seq(xlim[1], xlim[2], length.out = n)
  yvalues <- seq(ylim[1], ylim[2], length.out = n)
  xyfun <- function(x, y) {
    linkinv (beta[1] + beta[2] * x) *
      linkinv (beta[3] + beta[4] * y)
  }

  zvalues <- outer (xvalues, yvalues, FUN = xyfun)
  zvalues <- lambda * (alpha + (1 - alpha) * zvalues)

  # op <- par(bg = "white")
  res <- graphics::persp (x = xvalues, y = yvalues, z = zvalues,
                          theta = theta, phi = phi, expand = expand,
                          ltheta = ltheta, shade = shade, main = main,
                          xlab = xlab, ylab = ylab, zlab = zlab,
                          zlim = zlim, r = r, d = d, scale = scale,
                          border = border, box = box, axes = axes,
                          nticks = nticks, ticktype = ticktype,
                          col = col)

  return(invisible(res))
}
