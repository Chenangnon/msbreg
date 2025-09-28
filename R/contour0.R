
contour0.msbm <- function (x,
                           xlim = NULL, ylim = NULL, zlim,
                           xlab = NULL, ylab = NULL,
                           zlab = NULL,
                           panel = "panel.contourplot",
                           default.prepanel = "prepanel.default.levelplot",
                           cuts = 10,
                           nlevels = 200,
                           labels = TRUE,
                           contour = TRUE,
                           pretty = TRUE,
                           region = TRUE,
                           colorkey = list(title = zlab),
                           ...) {
  object <- x; rm(x)
  dims <- object$dims

  if ((dims$q == 2 & all(dims$pj == 1)) | (dims$q == 1 & all(dims$pj == 2))) {

    if (missing(xlim)) {
      X1val <- with(object$data, eval(object$call$formula[3][[1]][[2]]))
      X2val <- with(object$data, eval(object$call$formula[3][[1]][[3]]))

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

    # zlab
    if (is.null(zlab)) {
      yname <- deparse(object$call$formula[2][[1]])
      yname <- get.ycolname(yname)

      if (all(object$prior.weights == 1)) {
        zlab <- paste0("Prob(", yname, " = 1|", x1lab, ", ", x2lab, ")")
        if (nchar(zlab) > 25) {
          zlab <- paste0("Prob(", yname, " = 1)")
        }
      }
      else {
        zlab <- paste0("Mean(", yname, "|", x1lab, ", ", x2lab, ")")
        if (nchar(zlab) > 25) {
          zlab <- paste0("Mean(", yname, ")")
        }
      }
    }

    # Plot data
    x1values <- seq(x1lim[1], x1lim[2], length.out = nlevels)
    x2values <- seq(x2lim[1], x2lim[2], length.out = nlevels)
    data <- expand.grid(x = x1values, y = x2values)
    colnames(data) <- c(deparse(object$call$formula[3][[1]][[2]]),
                        deparse(object$call$formula[3][[1]][[3]]))
    Lo <- object$alpha.values
    L <- object$lambda.values

    object$frame <- object$frame %||% model.frame (object)
    if (dims$q == 1) {
      beta <- object$coefficients
      if (!object$frame$intercepts[1]) {
        beta <- c(0, beta)
      }

      msbfunc <- function(x) {
        y <- x[2]
        x <- x[1]
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

      msbfunc <- function(x) {
        y <- x[2]
        x <- x[1]
        object$link$linkinv (beta[1] + beta[2] * x) *
          object$link$linkinv (beta[3] + beta[4] * y)
      }
    }

    data$yname <- apply(data, MARGIN = 1, FUN = msbfunc)
    data$yname <- L[1] * (Lo[1] + (1 - Lo[1]) * data$yname)
    plotformula <- paste0("yname ~ ",
                          deparse(object$call$formula[3][[1]][[2]]), " + ",
                          deparse(object$call$formula[3][[1]][[3]]))
    plotformula <- as.formula(plotformula)

    if (missing(zlim)) {
      zlim <- c(0, 1)
    }
    else if (is.null(zlim)) {
      zlim <- c(0, 1)
    }

    # Plot
    out <- lattice::contourplot(plotformula, data = data,
                                panel = panel,
                                default.prepanel = default.prepanel,
                                xlim = x1lim, xlab = x1lab,
                                ylim = x2lim, ylab = x2lab,
                                zlim = zlim,
                                colorkey = colorkey,
                                cuts = cuts,
                                labels = labels,
                                region = region,
                                contour = contour,
                                pretty = pretty,
                                ...)

    print(out)

    return(invisible(out))
  }

  stop("not implemented for MSB model of this structure.")
}

