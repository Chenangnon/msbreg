#'
#' @importFrom graphics persp
#' @rdname courbe
#' @exportS3Method msbreg::courbe
courbe.msmb <- function (object, from = NULL, to = NULL, n = 101,
                         add = FALSE, type = "l", xname = NULL, yname = NULL,
                         xlab = NULL, ylab = NULL, zlab = NULL, log = NULL,
                         xlim = NULL, ylim = NULL, zlim = NULL,
                         main = " ", sub = NULL,
                         col = NULL, lty = 1, lwd = 1, cex = 0.1,
                         jitter = TRUE, jitter.range = .2, seed = NULL,
                         theta = 50, phi = 10, expand = 0.8,
                         ltheta = 0, shade = 0.75, r = sqrt(3),
                         d = 1, scale = TRUE, border = NULL,
                         box = TRUE, axes = TRUE, nticks = 5,
                         ticktype = "detailed", ...) {

  dims <- object$dims
  if (dims$q == 1 & all(dims$pj == 1)) {
    Xval <- with(object$data,
                 eval(object$call$formula[3][[1]]))
    if (missing(xlim)) {
      xlim <- range(Xval)
    }

    if (is.null(xlab))
      xlab <- deparse(object$call$formula[3][[1]])
    if (is.null(xname))
      xname <- xlab
    if (is.null(ylab)) {
      ylab <- paste0("Prob(",
                     object$ylab[2],
                     " = 1|", xlab, ")")
    }

    if (jitter) {
      stopifnot(is.numeric(jitter.range))
      eval(setsave.RNGstate())
      Jitter <- runif(NROW(object$data), min = -abs(jitter.range)/2, max = abs(jitter.range)/2)
    }
    else
      Jitter <- 0

    if (add) {
      points(x = Xval,
             y = with(object$data, eval(object$call$formula[2][[1]])) + Jitter,
             cex = cex, xlab = xlab, ylab = ylab)
    }
    else {
      plot(Xval,
           with(object$data, eval(object$call$formula[2][[1]])) + Jitter,
           xlim = xlim,
           cex = cex, xlab = xlab, ylab = ylab)
    }
    Lo <- object$alpha.values
    L <- object$lambda.values
    beta <- object$coefficients
    msb1func <- function(x) {
      L[1] * (Lo[1] + (1 - Lo[1]) * object$link$linkinv (beta[1] + beta[2] * x))
    }

    if (missing(col)) {
      col <- "black"
    }

    res <- graphics::curve(msb1func,
                           from = xlim[1], to = xlim[2], ylim = ylim,
                           col = col, lty = lty, lwd = lwd, xname = xname,
                           xlab = xlab, ylab = ylab, add = TRUE, n = n, ...)

    return(invisible(res))

  }

  if ((dims$q == 2 & all(dims$pj == 1)) | (dims$q == 1 & all(dims$pj == 2))) {

    X1val <- with(object$data, eval(object$call$formula[3][[1]][[2]]))
    X2val <- with(object$data, eval(object$call$formula[3][[1]][[3]]))

    if (missing(xlim)) {
      x1lim <- range(X1val)
      x2lim <- range(X2val)
    }
    else {
      if (is.list(xlim)) {
        x1lim <- xlim[[1]]
        x2lim <- xlim[[2]]
      }
      else {
        x1lim <- x2lim <- xlim
      }
    }

    if (is.null(xlab)) {
      x1lab <- deparse(object$call$formula[3][[1]][[2]])
      x2lab <- deparse(object$call$formula[3][[1]][[3]])
    }
    else {
      if (is.list(xlab)) {
        x1lab <- xlab[[1]]
        if (length(xlab) > 1)
          x2lab <- xlab[[2]]
        else
          x2lab <- x1lab
      }
      else {
        x1lab <- x2lab <- xlab
      }
    }
    if (is.null(ylab))
      ylab <- deparse(object$call$formula[2][[1]])

    if (jitter) {
      stopifnot(is.numeric(jitter.range))
      eval(setsave.RNGstate())
      Jitter <- runif(NROW(object$data), min = -abs(jitter.range)/2, max = abs(jitter.range)/2)
    }
    else
      Jitter <- 0

    x1values <- seq(x1lim[1], x1lim[2], length.out = n)
    x2values <- seq(x2lim[1], x2lim[2], length.out = n)

    Lo <- object$alpha.values
    L <- object$lambda.values

    object$frame <- object$frame %||% model.frame (object)
    if (dims$q == 1) {
      beta <- object$coefficients
      if (!object$frame$intercepts[1]) {
        beta <- c(0, beta)
      }

      msbfunc <- function(x, y) {
        object$link$linkinv (beta[1] + beta[2] * x + beta[3] * y)
      }
    }
    else {
      beta <- object$coefficients
      if (sum(object$frame$intercepts) == 0) {
        beta <- c(0, beta[1], 0, beta[2])
      }
      else {
        if (!object$frame$intercepts[1]) {
          beta <- c(0, beta[1:3])
        }
        if (!object$frame$intercepts[1]) {
          beta <- c(beta[1:2], 0, beta[3])
        }
      }

      msbfunc <- function(x, y) {
        object$link$linkinv (beta[1] + beta[2] * x) *
          object$link$linkinv (beta[3] + beta[4] * y)
      }
    }

    zvalues <- outer(x1values, x2values, FUN = msbfunc)
    zvalues <- L[1] * (Lo[1] + (1 - Lo[1]) * zvalues)

    if (missing(zlim)) {
      zlim <- c(0, 1)
    }
    else if (is.null(zlim)) {
      zlim <- c(0, 1)
    }

    # op <- par(bg = "white")
    if (missing(col)) {
      col <- "white"
    }

    if (is.null(zlab)) {
      zlab <- paste0("Prob(",
                     object$frame$auxy$ycolname,
                     " = 1|", x1lab, ", ", x2lab, ")")
    }

    res <- persp(x1values, x2values, zvalues, theta = theta, phi = phi,
                 expand = expand, col = col, ltheta = ltheta, shade = shade,
                 main = main, sub = sub, xlab = x1lab, ylab = x2lab, zlab = zlab, zlim = zlim,
                 r = r, d = d, scale = scale, border = border, box = box,
                 axes = axes, nticks = nticks, ticktype = ticktype, ...)

    return(invisible(res))

  }

  ###############################################################
  ###############################################################


  has.covariate <- dims$pj > 0
  nb.stg.has.covariate <- sum(has.covariate)
  stopifnot(nb.stg.has.covariate <= 2) # At most two stages with covariates
  if (nb.stg.has.covariate == 2) {
    stopifnot(all(dims$pj <= 1)) #
  }



  warning("Not implemented for MSB model of this structure: using only the 'coefficients' component of 'object'.")

  return(courbe_default_add(object$coefficients,
                            from = from, to = to, n = n,
                            add = add, type = type, xname = xname, yname = yname,
                            xlab = xlab, ylab = ylab, zlab = zlab, log = log,
                            xlim = xlim, ylim = ylim, zlim = zlim,
                            linkinv = object$link$linkinv, alpha = object$alpha.values[1],
                            main = main, col = col, lty = lty, lwd = lwd, cex = cex,
                            theta = theta, phi = phi, expand = expand,
                            ltheta = ltheta, shade = shade, r = r,
                            d = d, scale = scale, border = border,
                            box = box, axes = axes, nticks = nticks,
                            ticktype = ticktype, ...))
}
