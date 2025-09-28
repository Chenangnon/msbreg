
drop.NAs.from.x.y.weights.offset <- function() {
  expression({
    # Remove NAs from x, y, weights, and offset (if any)
    nobs <- NROW(x)
    rnames <- names(y)
    if (is.null(rnames)) {
      rnames <- 1:nobs
      if (length(y)) {
        names(y) <- rnames
      }
    }

    if (is.null(weights))
      weights <- 1

    if (is.null(offset))
      offset <- 0

    yx <- na.omit(cbind(y, weights, offset, x))
    dnaaction <- attr(yx, 'na.action')
    if (length(dnaaction)) {
      if (length(y)) {
        y <- yx[,1]
        weights <- yx[, 2]
        offset <- yx[, 3]
        x <- yx[, -(1:3), drop = FALSE]
      }
      else {
        weights <- yx[, 1]
        offset <- yx[, 2]
        x <- yx[, -(1:2), drop = FALSE]
      }
    }
    else {
      if (length(y)) {
      weights <- yx[, 2]
      offset <- yx[, 3]
      }
      else {
        weights <- yx[, 1]
        offset <- yx[, 2]
      }
    }
    rm(yx)
  })
}
