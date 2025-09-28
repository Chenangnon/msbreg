
extract.model.frame <- function() {
  expression({

    dnaaction <- attr(object, "na.action")

    Y <- stats::model.response (object, type = "any")

    if (!missing(data))
      X <- stats::model.matrix (object, data = data, ...)
    else
      X <- stats::model.matrix (object, ...)
    if (length(dnaaction))
      X <- X[-dnaaction,]

    if (length(dim(Y)) == 1L) {
      nm <- rownames(Y)
      dim(Y) <- NULL
      if (!is.null(nm))
        names(Y) <- nm
    }

    nobs <- NROW(X)
    stopifnot(nobs > 0)
    weights <- as.vector(model.weights(object))
    if (!is.null(weights)) {
      if (!is.numeric(weights))
        stop("'weights' must be a numeric vector")
      if (any(weights < 0))
        stop("negative weights not allowed")
      if (length(weights) != nobs)
        stop(gettextf("number of weights is %d should equal %d (number of observations)",
                      length(weights), nobs), domain = NA)
    }
    else {
      weights <- 1
    }

    if (NROW(Y) > 0) {
      if (NCOL(Y) == 1) {
        if (is.factor(Y))
          Y <- Y != levels(Y)[1L]
        Y[weights == 0] <- 0
        if (any(na.omit(Y) < 0 | na.omit(Y) > 1))
          stop("y values must be 0 <= y <= 1")
        my <- Y * weights
        Y <- round(my)
        if (any(abs(na.omit(my - Y)) > 0.001))
          warning(gettextf("non-integer #successes in a %s glm!",
                           "binomial"), domain = NA)
        ntrials <- weights
      }
      else if (NCOL(Y) == 2) {
        if (any(abs(Y - round(Y)) > 0.001))
          warning(gettextf("non-integer counts in a %s glm!",
                           "binomial"), domain = NA)
        n <- (Y1 <- Y[, 1L]) + Y[, 2L]
        Y <- Y1
        if (any(n0 <- n == 0))
          Y[n0] <- 0
        Y <- Y * weights
        ntrials <- weights * n
      }
      else stop(gettextf("for the '%s' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is #successes and col 2 is #failures",
                         "binomial"), domain = NA)

    if (length(ntrials) > 1)
      names(ntrials) <- names(Y)
    }
    else {
      ntrials <- 1
    }

    if (NROW(X) != nobs) { # Never TRUE ! (since nobs is now defined as NROW(X))
      if (!is.null(rownames(object))) {
        keep <- rownames(X) %in% rownames(object)
        X <- X[keep, , drop = FALSE]
        rownames(X) <- rownames(object)
      }
      else {
        if (length(dnaaction)) {
          notkeep <- (1:NROW(X)) %in% dnaaction
          X <- X[!notkeep, , drop = FALSE]
        }
      }
    }

    Offset <- stats::model.offset(object)
    if (is.null(Offset)) {
      Offset <- 0
    }
    else {

      if (length(dnaaction) & length(Offset) > 1)
        Offset <- Offset[-dnaaction]

      if (length(Offset) != nobs)
        stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                      length(Offset), nobs), domain = NA)
    }
  })
}
