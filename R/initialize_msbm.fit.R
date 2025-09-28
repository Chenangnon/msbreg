
# @export initialize.msbm.fit
#'
#' @importFrom stats lm.wfit

##* We could replace 'stats::lm.wfit' by 'stats::glm.fit'
#*  But we are not using 'stats::glm.fit' for two reasons:
#*   (1) we would need rational calculus algorithm to make
#*       use of the argument 'sample.weights' in 'stats::glm.fit'.
#*
#*       For instance, the package 'fractional' can be used to turn
#*       each 'sample.weights' into an approximate rational fraction
#*       ('fractional::fractional').
#*       Then, finding a common multiple for rational 'sample.weights'
#*       will help combine 'weights' with 'sample.weights' during the
#*       initialization.
#*
#*       See for e.g. 'as.numeric(fractional::fractional(w, sync = TRUE)) /
#*                      sum(fractional::fractional(w, sync = TRUE))'
#*       for a vector of weight 'w'.
#*
#*       Meanwhile, making use of the argument 'sample.weights' is a very
#*       simple task using 'stats::lm.wfit'.
#*
#*   (2) In every cases, 'stats::glm.fit' would take more time!
#*       We could reduce the number of iterations so that the difference
#*       in time is small, but then, 'stats::glm.fit' not converging does
#*       not warrant any result better than 'stats::lm.wfit'.
#*
##* We could also replace 'stats::lm.wfit' by 'betareg::betareg'
#*  But, is the dependence on the package 'betareg' worth it?
#*  This would also take significantly more time, unless only a
#*  few iterations are allowed (in that case, is it worth anything at all?).
#*
#*
initialize.msbm.fit <- function (frame, link = logit(),
                                 control = msbm.control(),
                                 bin.aug = 0.5) {

  #* Check that frame is an 'msbm.frame' object
  stopifnot(is.msbm.frame(frame))

  #* Set dimensions from mbm.frame: q, pj (for j = 1:q), d, r
  eval(get.msbm.dims())

  # Set alpha to an offset
  if (d == 0) {
    if (attr(frame$alpha.dictionary, 'offset')) {
      eta.alpha <- frame$offset.matrix[, attr(frame$alpha.dictionary,
                                              'offset.index')]
      alpha <- link$linkinv (eta.alpha)
    }
    else {
      eta.alpha <- -Inf
      alpha <- rep(0, length.out = nobs)
    }
  }
  else {
    if (attr(frame$alpha.dictionary, 'offset')) {
      eta.alpha <- frame$offset.matrix[, attr(frame$alpha.dictionary,
                                              'offset.index')]
      alpha <- link$linkinv (eta.alpha)
    }
    else
      alpha <- rep(0, length.out = nobs)
  }

  # Set lambda if constant (no covariate)
  if (r == 0) {
    if (attr(frame$lambda.dictionary, 'offset')) {
      eta.lambda <- frame$offset.matrix[, attr(frame$lambda.dictionary,
                                               'offset.index')]
      lambda <- link$linkinv (eta.lambda)
    }
    else {
      eta.lambda <- Inf
      lambda <- rep(1, length.out = nobs)
    }
  }
  else {
    if (attr(frame$lambda.dictionary, 'offset')) {
      eta.lambda <- frame$offset.matrix[, attr(frame$lambda.dictionary,
                                               'offset.index')]
      lambda <- link$linkinv (eta.lambda)
    }
    else
      lambda <- rep(1, length.out = nobs)
  }

  # Any conflict between alpha and lambda?
  if ((d == 0 | r == 0) & any(lambda <= alpha)) {
    stop("'lambda <= alpha': conflicting offsets in 'alpha' and 'lambda' components")
  }

  # only.offset
  only.offset <- sapply(frame$stage.dictionary,
                        FUN = function(dict) {
                          attr(dict, 'only.offset')
                        })

  # 'mustart', 'muhat' and 'logity'
  mustart <- (frame$y + bin.aug)/(frame$weights + 2*bin.aug)
  if (d == 0) {
    mustart <- pmax(mustart, alpha)
  }
  if (r == 0) {
    mustart <- pmin(mustart, lambda)
  }
  correctal <- mustart / lambda <= alpha
  if (any(correctal)) {
    if (d > 0) {
      alpha[correctal] <- pmax(0, pmin(alpha[correctal], mustart[correctal] / lambda[correctal] - control$epsilon))
    }
    else if (r > 0) {
      lambda[correctal] <- pmax(0, pmin(lambda[correctal], mustart[correctal] / alpha[correctal] - control$epsilon))
    }
    correctal <- mustart / lambda <= alpha
    if (any(correctal)) {
      mustart <- pmax(mustart, lambda * alpha)
    }
  }

  mustart -> muhat
  Pstart <- (mustart / lambda - alpha) / (1 - alpha)
  Phat <- Pstart
  logity <- link$linkfun (Phat^(1/sum(!only.offset)))
  #* Code used for 'glm.fit'
  #* muhat <- frame$y / frame$weights
  #* muhat[muhat == 1] <- 1 - control$epsilon
  #* muhat[muhat == 0] <- control$epsilon
  #* Phat <- (muhat / lambda - alpha) / (1 - alpha)
  #* logity <- link$linkfun (Phat^(1/sum(!only.offset)))

  #** Initialize betas (assume that alpha/lambda are 0/1 or set to offsets)
  #* Step 1.1: Roll over stages to initialize each betaj and get tildemuj
  ffun <- function(j) {
    dict <- frame$stage.dictionary[[j]]

    offsetj <- if (attr(dict, 'offset'))
      frame$offset.matrix[, attr(dict, 'offset.index')]
    else 0

    if (attr(dict, 'only.offset')) {
      muj <- rep(link$linkinv (offsetj),
                length.out = nobs)
      return(list(betaj = NULL, tildemuj = muj))
    }

    if (pj[j] > 0) {
      X <- frame$input.matrix[, dict, drop = FALSE]

      if (attr(dict, 'intercept')) {
        X <- cbind(1, X)
      }
    }
    else {
      X <- matrix(1, nrow = nobs, ncol = 1)
    }

    auxreg <- catch.conditions({
      stats::lm.wfit(x = X, y = logity,
                     w = rep(frame$sample.weights, length.out = nobs),
                     offset = rep(offsetj, length.out = nobs))

    })$value

    if (any(class(auxreg) %in% c("simpleError", "error", "condition", "try-error"))) {
      betaj <- rep(0, pj[j] + attr(dict, 'intercept'))

      muj <- rep(link$linkinv (offsetj),
                 length.out = nobs)

    }
    else {
      betaj <- auxreg$coefficients
      if (any(is.na(betaj))) {
        betaj[is.na(betaj)] <- 0
      }

      muj <- link$linkinv (auxreg$fitted.values)
    }

    return(list(betaj = betaj, tildemuj = muj))
  }
  beta.init <- lapply(1:q, FUN = ffun)
  tildemu <- sapply(beta.init,
                    FUN = function(x) {
                      x$tildemuj
                    })

  #* Step 1.2: Get bi and barmuj
  barmuj <- tildemu
  if (any(only.offset)) {
    bi <- (log(Phat) - rowSums(log(tildemu[,only.offset, drop = FALSE]))) /
      rowSums(log(tildemu[,!only.offset, drop = FALSE]))
  }
  else {
    bi <- log(Phat) / rowSums(log(tildemu[,!only.offset, drop = FALSE]))
  }
  barmuj[,!only.offset] <- tildemu[,!only.offset]^bi

  #* Step 1.3: Get hatmuj and betas
  ffun <- function(j) {
    dict <- frame$stage.dictionary[[j]]

    offsetj <- if (attr(dict, 'offset'))
      frame$offset.matrix[, attr(dict, 'offset.index')]
    else 0

    if (attr(dict, 'only.offset')) {
      muj <- rep(link$linkinv (offsetj),
                 length.out = nobs)
      return(list(betaj = NULL, hatmuj = muj))
    }

    if (pj[j] > 0) {
      X <- frame$input.matrix[, dict, drop = FALSE]

      if (attr(dict, 'intercept')) {
        X <- cbind(1, X)
        colnames(X)[1] <- paste0("(Intercept).", j)
      }
    }
    else {
      X <- matrix(1, nrow = nobs, ncol = 1)
      colnames(X) <- paste0("(Intercept).", j)
    }

    baretaj <- link$linkfun (barmuj[,j]) # Could we factor in measurement errors?
    #                                    # Would be too complex (in the best case, link$linkfun (barmuj[,j]) has a denominator which depends on beta_j)
    #                                    # But we could easily include ME offsets (if any!) since the denominator does not depend on beta_j
    auxreg <- catch.conditions({
      stats::lm.wfit(x = X, y = baretaj,
                     w = rep(frame$sample.weights, length.out = nobs),
                     offset = rep(offsetj, length.out = nobs))

    })$value

    if (any(class(auxreg) %in% c("simpleError", "error", "condition", "try-error"))) {
      betaj <- rep(0, pj[j] + attr(dict, 'intercept'))

      muj <- rep(link$linkinv (offsetj),
                 length.out = nobs)

    }
    else {
      betaj <- auxreg$coefficients
      if (any(is.na(betaj))) {
        betaj[is.na(betaj)] <- 0
      }

      muj <- link$linkinv (auxreg$fitted.values)
    }

    return(list(betaj = betaj, hatmuj = muj))
  }
  betas <- lapply(1:q, FUN = ffun)
  hatmuj <- sapply(betas,
                    FUN = function(x) {
                      x$hatmuj
                    })
  betas <- lapply(betas,
                  FUN = function(x) {
                    x$betaj
                  })

  #** Initialize lambda (assume that alpha = 0 or set to offset)
  Phat <- exp(rowSums(log(hatmuj)))
  if (r > 0) {
    L.init <- pmin(muhat / (alpha + (1 - alpha) * Phat), 1-control$epsilon)
    etaL.init <- link$linkfun (L.init)
    dict <- frame$lambda.dictionary

    offset.L <- if (attr(dict, 'offset'))
      frame$offset.matrix[, attr(dict, 'offset.index')]
    else 0

    if (r - attr(dict, 'intercept') > 0) {
      X <- frame$input.matrix[, dict, drop = FALSE]

      if (attr(dict, 'intercept')) {
        X <- cbind(1, X)
        colnames(X)[1] <- "(Intercept).lambda"
      }
    }
    else {
      X <- matrix(1, nrow = nobs, ncol = 1)
      colnames(X) <- "(Intercept).lambda"
    }

    auxreg <- catch.conditions({
      stats::lm.wfit(x = X, y = etaL.init,
                     w = rep(frame$sample.weights, length.out = nobs),
                     offset = offset.L)

    })$value


    if (any(class(auxreg) %in% c("simpleError", "error", "condition", "try-error"))) {
      delta <- rep(0, r + attr(dict, 'intercept'))

      hatL <- rep(link$linkinv (offset.L),
                 length.out = nobs)

    }
    else {
      delta <- auxreg$coefficients
      if (any(is.na(delta))) {
        delta[is.na(delta)] <- 0
      }

      hatL <- link$linkinv (auxreg$fitted.values)
    }
  }
  else {
    delta <- NULL
    hatL <- lambda
  }

  #** Initialize alpha
  #* Step 2.1
  if (d > 0) {
    A.init <- pmax(pmin((muhat / hatL - Phat) / (1 - Phat), 1 - control$epsilon), control$epsilon)
    etaA.init <- link$linkfun (A.init)
    dict <- frame$alpha.dictionary

    offset.A <- if (attr(dict, 'offset'))
      frame$offset.matrix[, attr(dict, 'offset.index')]
    else 0

    if (d - attr(dict, 'intercept') > 0) {
      X <- frame$input.matrix[, dict, drop = FALSE]

      if (attr(dict, 'intercept')) {
        X <- cbind(1, X)
        colnames(X)[1] <- "(Intercept).alpha"
      }
    }
    else {
      X <- matrix(1, nrow = nobs, ncol = 1)
      colnames(X) <- "(Intercept).alpha"
    }

    auxreg <- catch.conditions({
      stats::lm.wfit(x = X, y = etaA.init,
                     w = rep(frame$sample.weights, length.out = nobs),
                     offset = offset.A)

    })$value


    if (any(class(auxreg) %in% c("simpleError", "error", "condition", "try-error"))) {
      gamma <- rep(0, d + attr(dict, 'intercept'))

      hatA <- rep(link$linkinv (offset.A),
                  length.out = nobs)

    }
    else {
      gamma <- auxreg$coefficients
      if (any(is.na(gamma))) {
        gamma[is.na(gamma)] <- 0
      }

      hatA <- link$linkinv (auxreg$fitted.values)
    }
  }
  else {
    gamma <- NULL
    hatA <- alpha
  }

  structure(c(unlist(betas), gamma, delta),
            beta = betas,
            gamma = gamma,
            delta = delta,
            mujhat = hatmuj,
            alphahat = hatA,
            lambdahat = hatL)
}
