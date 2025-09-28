
#' @importFrom lattice contourplot
#' @importFrom lattice lattice.options
# @importFrom stringi stri_length
#'
#' @rdname contour.msbm
#' @exportS3Method graphics::contour
contour.glm <- function (x,
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

  if (inherits(object, "glm"))
    p <- length(attr(object$terms,"dataClasses")) - 1
  else
    p <- object$p

  if (p != 2)
    stop("'contour.glm' only handles bivariate model fits")

  # Predictor ranges
  if (!length(xlim)) {
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
  stopifnot(is.numeric(x1lim), is.numeric(x2lim),
            length(x1lim) == 2, length(x2lim) == 2)

  # Predictor labels
  if (!length(xlab)) {
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

    if (identical(object$family$family, "binomial") & all(object$prior.weights == 1)) {
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
  x <- seq(from = x1lim[1], to = x1lim[2], length.out = nlevels)
  y <- seq(from = x2lim[1], to = x2lim[2], length.out = nlevels)

  zclass <- attr(object$zmaxpterms,"dataClasses")
  if (length(zclass)) {
    ### ???
  }

  data <- expand.grid(x=x, y=y)
  colnames(data) <- c(deparse(object$call$formula[3][[1]][[2]]),
                      deparse(object$call$formula[3][[1]][[3]]))
  data$yname <- predict(object, newdata = data, type = "response")
  plotformula <- paste0("yname ~ ",
                        deparse(object$call$formula[3][[1]][[2]]), " + ",
                        deparse(object$call$formula[3][[1]][[2]]))
  plotformula <- as.formula(plotformula)

  # zlim
  if(missing(zlim))
    zlim <- range(data$yname, finite = TRUE)

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
