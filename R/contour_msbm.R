
#' Display Fitted Model Contours
#'
#' Create color contour plots for fitted model objects.
#'
#' @param x an object of class (or inheriting from class) "msbm".
#'
#' @param xlim,ylim,zlim,xlab,ylab,zlab,xname,yname labels for the contours.
#'
#' @param panel,default.prepanel See \link[lattice]{contourplot}.
#'
#' @param cuts number of contour levels desired, i.e. the number of levels the
#' range of z would be divided into..
#'
#' @param nlevels number of function evaluation points per covariate.
#'
#' @param labels,contour,pretty,region logicals, see \link[lattice]{contourplot}.
#'
#' @param colorkey A logical flag specifying whether a colorkey is to be drawn
#' alongside the plot, or a list describing the colorkey. See
#' \link[lattice]{contourplot} for a description of the components of a list
#' \code{colorkey}.
#'
#' @param ... Further arguments passed \link[lattice]{contourplot}.
#'
#' @details
#' The function is a \link{contour} method for \code{"msbm"} objects.
#' It uses \link[lattice]{contourplot} from package \code{lattice}.
#'
#' @aliases contours
#'
#' @importFrom lattice contourplot
#' @importFrom lattice lattice.getOption
#' @importFrom lattice lattice.options
#' @importFrom grDevices contourLines
#' @importFrom graphics points
#' @return An object of class \code{"trellis"} (returned invisibly).
#' See \link[lattice]{contourplot} for a description.
#'
#' @exportS3Method graphics::contour
contour.msbm <- function (x,
                          xname, yname,
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

  object$frame <- object$frame %||% model.frame(object)
  object$frame$frames <- object$frame$frames %||% update(object$frame, frames = TRUE)$frames
  mainframe <- object$frame$frames$maim

  newdata <- expand.grid()
  newdata$zname <- predict.msbm(object, newdata = newdata)

}
