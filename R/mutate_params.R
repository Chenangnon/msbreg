#* Transform parameter vector into various quantities entering
#* the log-likelihood function or its derivarives (score) or
#* hessian matrix (negative information)
#*

# offset in alpha.matrix must be -Inf if alpha = 0 (set this manually
# after calling msbm.frame)
# Same for lambda = 1 (offset = Inf)
#
# Use pracma::newtonsys for solving the system: score = 0?
#
# Review the use of 'replicate': no need most of times
#
# Organization of vector and matrix (and related operations)
# have been designed to be optimal when the number of observations
# is larger than the number of MSB stages and the number of model parameters
mutate.params <- function (theta, frame, link = "logit",
                           score = FALSE, information = observed,
                           observed = FALSE, jerk = FALSE,
                           direct.lambda = FALSE,
                           exp.lambda = FALSE) {

  # Get Binomial link function
  eval(get.linkfun())

  # Compute basic quantities and the vector 'mu' of success probabilities
  eval(basic.mutate())

  # Re-compute the lambda component and mu if exp.lambda = TRUE
  if (exp.lambda) {
    if (r == 0) {
      if (attr(frame$lambda.dictionary, 'offset')) {
        loglambda <- eta.lambda
      }
      else {
        eta.lambda <- 0
        loglambda <- 0
      }
    }
    else {
      loglambda <- eta.lambda
    }

    # Compute the vector 'mu' of success probabilities
    mu <- exp(loglambda + logalpha) + (1 - exp(logalpha)) * exp(loglambda + rowSumlogmus)
    validmu <- is.finite(mu) & (mu >= 0 & mu <= 1)
  }

  #* Compute log-likelihood values
  size <- frame$weights
  if (length(size) > 1) {
    size <- size[validmu]
  }
  sweights <- frame$sample.weights
  if (length(sweights) > 1) {
    sweights <- sweights[validmu]
  }
  loglike.i <- numeric(nobs) + NA

  NotNULLy <- !is.null(frame$y)
  if (NotNULLy) {
    loglike.i[validmu] <- stats::dbinom (frame$y[validmu], size = size,
                                         prob = mu[validmu], log = TRUE)
    loglike.i[validmu] <- loglike.i[validmu] * sweights
    names(loglike.i) <- names(frame$y)
  }

  out <- list(loglike = sum(loglike.i[validmu]),
              loglike.i = loglike.i,
              mu = mu, validmu = validmu, logmus = logmus,
              etas = cbind(etas, alpha = eta.alpha, lambda = eta.lambda), SDs = SDs,
              logalpha = logalpha, loglambda = loglambda,
              betas = betas, gamma = gamma, delta = delta)

  if (!score & !information & !jerk) {
    return(out)
  }

  #* Compute score vector
  # Derivative of the inverse link function
  mu.eta <- link$mu.eta

  # hdot
  sigma_o2 <- link$sigma^2
  if (length(SDs) == 1) {
    Invfactors <- rep(1, nobs)
    Invfactorj <- function(j, drop = TRUE) {
      return(if (drop) Invfactors else cbind(Invfactors))
    }
  }
  else {
    Invfactors <- 1 / sqrt(1 + (SDs^2)/sigma_o2)
    Invfactorj <- function(j, drop = TRUE) {
      Invfactors[, j, drop = drop]
    }
  }
  hdots <- mu.eta (etas * Invfactors)

  #SB
  SBs <- NULL
  if (frame$me) {
    ffun <- function(j) {
      dict <- frame$stage.dictionary[[j]]

      if (!attr(dict, 'has.me')) {
        return(NULL)
      }
      betaj <- betas[[j]]

      if (attr(dict, 'intercept')) {
        betaj <- betaj[-1]
      }
      betaj <- betaj[attr(dict, 'contr.with.me')]
      SDj <- frame$sd.input.matrix[, attr(dict, 'me.index'), drop = FALSE]

      if (length(betaj) > 1) {
        SBj <- (SDj^2) * t(replicate(nobs, betaj))
      }
      else {
        SBj <- (SDj^2) * betaj[1]
      }

      return(SBj)
    }
    SBs <- lapply(1:q, FUN = ffun)
  }

  # Tau
  if (frame$me) {
    ffun <- function(j) {
      dict <- frame$stage.dictionary[[j]]

      if (attr(dict, 'only.offset')) {
        return(NULL)
      }

      if (pj[j] == 0) {
        # Only an intercept
        Xj <- matrix(1, nrow = nobs, ncol = 1) * Invfactorj(j, drop = FALSE)
        colnames(Xj) <- paste0('(Intercept).', j)
        return(Xj)
      }

      Xj <- frame$input.matrix[, dict, drop = FALSE]

      if (!attr(dict, 'has.me')) {
        if (attr(dict, 'intercept')) {
          cnm <- colnames(Xj)
          Xj <- cbind(1, Xj)
          colnames(Xj) <- c(paste0('(Intercept).', j), cnm)
        }

        if (!attr(dict, 'has.me.offset')) {
          return(Xj)
        }

        Tauj <- replicate( NCOL(Xj), Invfactorj(j) ) * Xj

        return(Tauj)
      }

      SBj <- SBs[[j]]

      n.mej <- length(attr(dict, 'contr.with.me'))
      XDiff <- replicate(n.mej, c(etas[,j] / (sigma_o2 + SDs[,j]^2)) ) * SBj

      Xj[, attr(dict, 'contr.with.me')] <-
        Xj[, attr(dict, 'contr.with.me'), drop = FALSE] - XDiff

      if (attr(dict, 'intercept')) {
        cnm <- colnames(Xj)
        Xj <- cbind(1, Xj)
        colnames(Xj) <- c(paste0('(Intercept).', j), cnm)
      }

      Tauj <- replicate( NCOL(Xj), Invfactorj(j) ) * Xj

      return(Tauj)
    }
  }
  else {
    ffun <- function(j) {
      dict <- frame$stage.dictionary[[j]]

      if (attr(dict, 'only.offset')) {
        return(NULL)
      }

      if (pj[j] == 0) {
        # Only an intercept
        Xj <- matrix(1, nrow = nobs, ncol = 1) * Invfactorj(j, drop = FALSE)
        colnames(Xj) <- paste0('(Intercept).', j)
        return(Xj)
      }

      Xj <- frame$input.matrix[, dict, drop = FALSE]

      if (attr(dict, 'intercept')) {
        cnm <- colnames(Xj)
        Xj <- cbind(1, Xj)
        colnames(Xj) <- c(paste0('(Intercept).', j), cnm)
      }

      if (!attr(dict, 'has.me.offset')) {
        return(Xj)
      }

      Tauj <- replicate(NCOL(Xj), Invfactorj(j)) * Xj
      return(Tauj)

    }
  }
  Taus <- lapply(1:q, FUN = ffun)

  # muj.betaj
  ffun <- function (j) {
    Tauj <- Taus[[j]]
    muj.betaj <- if (!is.null(Tauj)) replicate(NCOL(Tauj), c(hdots[,j])) * Tauj

    return(muj.betaj)
  }
  mu.betas <- lapply(1:q, FUN = ffun)

  # mu.beta
  ffun <- function (j) {
    muj.betaj <- mu.betas[[j]]

    if (is.null(muj.betaj)) {
      return(NULL)
    }

    mu_j <- (1 - exp(logalpha)) *
      exp(loglambda + rowSumlogmus - logmus[,j])

    mu.betaj <- replicate(NCOL(muj.betaj), mu_j) * muj.betaj

    return(mu.betaj)
  }
  mu.beta <- lapply(1:q, FUN = ffun)

  # mu.gamma
  if (d > 0) {
    hdot.alpha <- mu.eta (eta.alpha)
    c.alpha <- exp(loglambda) - exp(loglambda + rowSumlogmus)
    m.alpha <- replicate(d, c(c.alpha * hdot.alpha))
    if (attr(frame$alpha.dictionary, 'intercept')) {
      if (d == 1) {
        mu.gamma <- m.alpha
      }
      else {
        mu.gamma <- cbind(1, frame$input.matrix[, frame$alpha.dictionary, drop = FALSE]) * m.alpha
      }
      colnames(mu.gamma) <- c('(Intercept).alpha',
                                 colnames(frame$input.matrix[, frame$alpha.dictionary, drop = FALSE]))
    }
    else {
      mu.gamma <- frame$input.matrix[, frame$alpha.dictionary, drop = FALSE] * m.alpha
    }
  }
  else {
    mu.gamma <- NULL
  }

#  browser()

  # mu.delta
  if (r > 0) {
    if (direct.lambda) {
      hdot.lambda <- 1
    } else if (exp.lambda) {
      hdot.lambda <- exp(min(eta.lambda, 709))
    }
    else {
        hdot.lambda <- mu.eta (eta.lambda)
    }

    c.lambda <- exp(logalpha) + (1 - exp(logalpha)) * exp(rowSumlogmus)
    m.delta <- replicate(r, c(c.lambda * hdot.lambda))
    if (attr(frame$lambda.dictionary, 'intercept')) {
      if (r == 1) {
        mu.delta <- m.delta
      }
      else {
        mu.delta <- m.delta * cbind(1, frame$input.matrix[, frame$lambda.dictionary, drop = FALSE])
      }
      colnames(mu.delta) <- c('(Intercept).lambda',
                              colnames(frame$input.matrix[, frame$lambda.dictionary, drop = FALSE]))
    }
    else {
      mu.delta <- frame$input.matrix[, frame$lambda.dictionary, drop = FALSE] * m.delta
    }
  }
  else {
    mu.delta <- NULL
  }

  mu.theta <- cbind(do.call('cbind', mu.beta), mu.gamma, mu.delta)
  out$tau <- Taus
  out$mu.theta <- mu.theta
  out$c.alpha <- if (d > 0) c.alpha
  out$c.lambda <- if (r > 0) c.lambda

  if (NotNULLy) {
    c.score <- c(frame$sample.weights * (frame$y - frame$weights * mu) /
                   (mu * (1 - mu)))
    out$score.i <- mu.theta * replicate(NCOL(mu.theta), c.score)
    out$score <- colSums(out$score.i)

    rownames(out$score.i) <- names(frame$y)
    colnames(out$score.i) <- names(theta)
    if (is.null(colnames(out$score.i)))
      colnames(out$score.i) <- frame$parnames
    colnames(out$mu.theta) <- colnames(mu.theta) <- colnames(out$score.i)
    names(out$score) <- names(theta)
  }
  else if (score) {
    out$score <- rep(NA, length(theta))
    out$score.i <- matrix(NA, nrow = NROW(frame$input.matrix),
                          ncol = length(theta))
    names(out$score) <- rownames(out$score.i) <- rownames(frame$input.matrix)
    colnames(out$score.i) <- names(theta)
    warning("there is no response in the supplied model frame for score computation-returning NAs")
  }

  if (!information & !jerk) {
    return(out)
  }

  c.Einfo <- frame$sample.weights * frame$weights / (mu * (1 - mu))
  out$EFinfo <- t(out$mu.theta) %*% (replicate(NCOL(out$mu.theta), c.Einfo) * out$mu.theta)
  # Use Rfast::mat.mult(x, y) for speed? How reliable is 'Rfast' maintenance?

  dimnames(out$EFinfo) <- list(colnames(out$score.i),
                               colnames(out$score.i))
  if (!observed & !jerk) {
    return(out)
  }

  if (!NotNULLy) {
    warning("there is no response in the supplied model frame for observed information computation-returning NAs")
    out$OFinfo <- out$EFinfo + NA
    dimnames(out$OFinfo) <- dimnames(out$EFinfo)
    return(out)
  }

  #* Compute hessian matrix
  # Second derivative of the inverse link function
  mu.eta.eta <- link$mu.eta.eta

  # hddot
  if (length(SDs) == 1) {
    Numfactors <- rep(1 / (sigma_o2 + SDs^2), nobs)
    Numfactorj <- function(j, drop = TRUE) {
      return(if (drop) Numfactors else cbind(Numfactors))
    }
  }
  else {
    Numfactors <- 1 / (sigma_o2 + SDs^2)
    Numfactorj <- function(j, drop = TRUE) {
      Numfactors[, j, drop = drop]
    }
  }

  Inv2factors <- Invfactors * Numfactors
  if (length(SDs) == 1) {
    Inv2factorj <- function(j, drop = TRUE) {
      return(if (drop) Inv2factors else cbind(Inv2factors))
    }
  }
  else {
    Inv2factorj <- function(j, drop = TRUE) {
      Inv2factors[, j, drop = drop]
    }
  }
  hddots <- mu.eta.eta (etas * Invfactors)

  # dTau/dbeta (only components with MEs)
  Taudots <- NULL
  if (frame$me) {
    ffun <- function(j) {
      dict <- frame$stage.dictionary[[j]]

      if (!attr(dict, 'has.me')) {
        return(NULL)
      }

      Xj <- cbind(if (attr(dict, 'intercept')) 1,
                  frame$input.matrix[, dict, drop = FALSE])
      SDj <- frame$sd.input.matrix[, attr(dict, 'me.index'), drop = FALSE]
      ncj <- attr(dict, 'intercept') + pj[j]
      n.mej <- length(attr(dict, 'contr.with.me'))
      SBj <- SBs[[j]]
      #Tauj <- Taus[[j]]

      if (n.mej > 1) {
        SBBSj <- apply(SBj, MARGIN = 1,
                       FUN = function (x) {
                         c(outer(x, x))
                       })
        SBBSj <- t(SBBSj)

        outerSj <- apply(SDj^2, MARGIN = 1,
                         FUN = function (x) {
                           varmat <- diag(nrow = n.mej, ncol = n.mej)
                           diag(varmat) <- x
                           c(varmat)
                         })
        outerSj <- t(outerSj)
      }
      else {
        SBBSj <- SBj^2

        outerSj <- SDj^2
      }

      if (ncj > 1) {
        XjSBj <- apply(cbind(Xj, SBj),
                       MARGIN = 1,
                       FUN = function (x) {
                         xij <- x[1:ncj]
                         sij <- numeric(ncj)
                         sij[attr(dict, 'contr.with.me') +
                               attr(dict, 'intercept')] <- x[-(1:ncj)]
                         out <- outer(xij, sij)
                         c(out + t(out))
                       })
        XjSBj <- t(XjSBj)
      }
      else {
        # Only one covariate which also has ME
        XjSBj <- 2 * Xj * SBj
      }

      dTauj <- - replicate(ncj^2, Inv2factorj(j, drop = TRUE)) * XjSBj
      otherTerms <- 3 * replicate(n.mej^2, Numfactorj(j, drop = TRUE) ) * SBBSj - outerSj
      otherTerms <- replicate(n.mej^2, etas[,j, drop = TRUE] *
                                Inv2factorj(j, drop = TRUE)) * otherTerms

      # index to complete dtau for contrasts with MEs
      jpos <- attr(dict, 'contr.with.me') + attr(dict, 'intercept')
      me.pos <- matrix(0, nrow = ncj, ncol = ncj)
      me.pos[jpos, jpos] <- 1
      me.pos <- which(c(me.pos) > 0)

      dTauj[, me.pos] <- dTauj[, me.pos] + otherTerms

      return (dTauj)
    }
    Taudots <- lapply(1:q, FUN = ffun)
  }

  # muj.betaj.betaj (hessians of mu_ij)
  ffun <- function (j) {
    dict <- frame$stage.dictionary[[j]]
    if (attr(dict, 'only.offset')) {
      return(NULL)
    }
    Tauj <- Taus[[j]]
    n_c_tau <- NCOL(Tauj)

    if (n_c_tau > 1) {
      muj.betaj.betaj <- apply(Tauj, MARGIN = 1,
                               FUN = function(x) {
                                 c(outer(x, x))
                               })
      muj.betaj.betaj <- t(muj.betaj.betaj)
    }
    else {
      muj.betaj.betaj <- Tauj^2
    }

    muj.betaj.betaj <- replicate(NCOL(muj.betaj.betaj), hddots[,j, drop = TRUE]) * muj.betaj.betaj

    if (!attr(dict, 'has.me')) {
      return(muj.betaj.betaj)
    }

    Taudotj <- Taudots[[j]]
    Taudotj <- replicate(NCOL(Taudotj), hdots[,j, drop = TRUE]) * Taudotj

    muj.betaj.betaj <- muj.betaj.betaj + Taudotj

    return(muj.betaj.betaj)
  }
  mu.betas.betas <- lapply(1:q, FUN = ffun)

  # number of parameters per stage j (elements of betaj)
  n_j <- pj + intercepts

  # Indices for pairs of betaj components
  jk <- lapply(1:q, FUN = function(j) {
    cbind(j:q, j)
  })
  jk <- do.call('rbind', jk)

  # mu.beta.beta: hessian for each pair
  ffun <- function (l) {
    j <- jk[l,1]
    k <- jk[l,2]

    if (attr(frame$stage.dictionary[[j]], 'only.offset') |
        attr(frame$stage.dictionary[[k]], 'only.offset')) {
      return(NULL)
    }

    if (j == k) {
      # Diagonal block
      muj.betaj.betaj <- mu.betas.betas[[j]]
      mu_j <- (1 - exp(logalpha)) *
        exp(loglambda + rowSumlogmus - logmus[,j])

      mu.betaj.betaj <- replicate(NCOL(muj.betaj.betaj), mu_j) * muj.betaj.betaj

      return(mu.betaj.betaj)
    }
    else {
      # Off-diagonal blocks
      muj.betaj <- mu.betas[[j]]
      muk.betak <- mu.betas[[k]]
      mu_jk <- (1 - exp(logalpha)) *
        exp(loglambda + rowSumlogmus - logmus[,j] - logmus[,k])

      pjj <- NCOL(muj.betaj)
      mu.betajk <- apply(cbind(muj.betaj, muk.betak), MARGIN = 1,
                         FUN = function(x) {
                           c(outer(x[1:pjj], x[-(1:pjj)]))
                         })
      pjk <- pjj * NCOL(muk.betak)
      if (pjk > 1) {
        mu.betajk <- t(mu.betajk)
      }
      else {
        mu.betajk <- cbind(mu.betajk)
      }
      mu.betaj.betak <- replicate(pjk, mu_jk) * mu.betajk

      return(mu.betaj.betak)
    }
  }
  mu.beta.beta.list <- lapply(1:NROW(jk), FUN = ffun)

  # mu.gamma.gamma and related
  if (d > 0) {
    # mu.gamma.gamma
    hddot.alpha <- mu.eta.eta (eta.alpha)
    m.alpha <- replicate(d^2, c(c.alpha * hddot.alpha))
    if (attr(frame$alpha.dictionary, 'intercept')) {
      if (d == 1) {
        mu.gamma.gamma <- m.alpha
      }
      else {
        uut <- t(apply(frame$input.matrix[, frame$alpha.dictionary, drop = FALSE],
                       MARGIN = 1,
                       FUN = function (x) {
                         c(outer(c(1, x), c(1, x)))
                       }))
        mu.gamma.gamma <- m.alpha * uut
      }
    }
    else {
      if (d == 1) {
        uut <- frame$input.matrix[, frame$alpha.dictionary, drop = FALSE]^2
      }
      else {
        uut <- t(apply(frame$input.matrix[, frame$alpha.dictionary, drop = FALSE],
                       MARGIN = 1,
                       FUN = function (x) {
                         c(outer(x, x))
                       }))
      }
      mu.gamma.gamma <- m.alpha * uut
    }

    # mu.beta.gamma
    if (attr(frame$alpha.dictionary, 'intercept')) {
      if (d == 1)
        Xalpha <- 1
      else
        Xalpha <- cbind(1, frame$input.matrix[, frame$alpha.dictionary, drop = FALSE])
    }
    else {
      Xalpha <- frame$input.matrix[, frame$alpha.dictionary, drop = FALSE]
    }
    p_alpha <- NCOL(Xalpha)
    ffun <- function (j) {
      if (is.null(mu.betas[[j]]))
        return(NULL)
      muj.betaj.gamma <- apply(cbind(Xalpha, mu.betas[[j]]),
                               MARGIN = 1,
                               FUN = function(x) {
                                 c(outer(x[-(1:p_alpha)],
                                         x[1:p_alpha]))
                               })
      if (p_alpha > 1 | NCOL(mu.betas[[j]]) > 1) {
        muj.betaj.gamma <- t(muj.betaj.gamma)
      }
      else {
        muj.betaj.gamma <- cbind(muj.betaj.gamma)
      }

      mu_j <- - exp(loglambda + rowSumlogmus - logmus[,j]) * hdot.alpha
      mu.betaj.gamma <- replicate(NCOL(muj.betaj.gamma), mu_j) * muj.betaj.gamma

      return(mu.betaj.gamma)
    }
    mu.betaj.gamma.list <- lapply(1:q, FUN = ffun)

  }
  else {
    mu.betaj.gamma.list <- NULL
    mu.gamma.gamma <- NULL
  }

  # mu.delta.delta and related
  if (r > 0) {
    if (direct.lambda) {
      hddot.delta <- 0
    }
    else if (exp.lambda) {
      hddot.delta <- exp(min(eta.lambda, 709))
    }
    else {
      hddot.delta <- mu.eta.eta (eta.lambda)
    }
    m.delta <- replicate(r^2, c(c.lambda * hddot.delta))
    if (attr(frame$lambda.dictionary, 'intercept')) {
      if (r == 1) {
        mu.delta.delta <- m.delta
      }
      else {
        vvt <- t(apply(frame$input.matrix[, frame$lambda.dictionary, drop = FALSE],
                       MARGIN = 1,
                       FUN = function (x) {
                         c(outer(c(1, x), c(1, x)))
                       }))

        mu.delta.delta <- m.delta * vvt
      }
    }
    else {
      if (r == 1) {
        vvt <- frame$input.matrix[, frame$lambda.dictionary, drop = FALSE]^2
      }
      else {
        vvt <- t(apply(frame$input.matrix[, frame$lambda.dictionary, drop = FALSE], MARGIN = 1,
                     FUN = function (x) {
                       c(outer(x, x))
                     }))
      }
      mu.delta.delta <- m.delta * vvt
    }

    # mu.betaj.delta.list
    if (attr(frame$lambda.dictionary, 'intercept')) {
      if (r == 1)
        Xlambda <- 1
      else
        Xlambda <- cbind(1, frame$input.matrix[, frame$lambda.dictionary, drop = FALSE])
    }
    else {
      Xlambda <- frame$input.matrix[, frame$lambda.dictionary, drop = FALSE]
    }
    p_lambda <- NCOL(Xlambda)
    ffun <- function (j) {
      if (is.null(mu.betas[[j]]))
        return(NULL)
      muj.betaj.delta <- apply(cbind(Xlambda, mu.betas[[j]]),
                                 MARGIN = 1,
                                 FUN = function(x) {
                                   c(outer(x[-(1:p_lambda)],
                                           x[1:p_lambda]))
                                 })
      if (p_lambda > 1 | NCOL(mu.betas[[j]]) > 1) {
        muj.betaj.delta <- t(muj.betaj.delta)
      }
      else {
        muj.betaj.delta <- cbind(muj.betaj.delta)
      }

      mu_j <- (1 - exp(logalpha)) * exp(rowSumlogmus - logmus[,j]) * hdot.lambda

      mu.betaj.delta <- replicate(NCOL(muj.betaj.delta), mu_j) * muj.betaj.delta

      return(mu.betaj.delta)
    }
    mu.betaj.delta.list <- lapply(1:q, FUN = ffun)

    if (d > 0) {
      ffun <- function () {
        o.gamma.delta <- apply(cbind(Xlambda, Xalpha),
                                 MARGIN = 1,
                                 FUN = function(x) {
                                   c(outer(x[-(1:p_lambda)],
                                           x[1:p_lambda]))
                                 })
        if (p_alpha > 1 | p_lambda > 1) {
          o.gamma.delta <- t(o.gamma.delta)
        }
        else {
          o.gamma.delta <- matrix(c(o.gamma.delta), nrow = nobs, ncol = 1)
        }

        mu_j <- (1 - exp(rowSumlogmus)) * hdot.alpha * hdot.lambda

        omu.gamma.delta <- replicate(NCOL(o.gamma.delta), mu_j) * o.gamma.delta

        return(omu.gamma.delta)
      }
      mu.gamma.delta.list <- ffun()
    }
    else {
      mu.gamma.delta.list <- NULL
    }
  }
  else {
    mu.betaj.delta.list <- NULL
    mu.gamma.delta.list <- NULL
    mu.delta.delta <- NULL
  }

  # Build the total sample hessian
  ffun <- function (l) {
    j <- jk[l,1]
    k <- jk[l,2]

    if (is.null(mu.beta.beta.list[[l]])) {
      return(NULL)
    }

    tojk <- replicate(NCOL(mu.beta.beta.list[[l]]), c.score) *
      mu.beta.beta.list[[l]]

    hjk <- matrix(colSums(tojk),
                  nrow = n_j[j], ncol = n_j[k])
    return(hjk)
  }
  mu.beta2.list <- lapply(1:NROW(jk), FUN = ffun)

  ffun <- function (k) {
    if (attr(frame$stage.dictionary[[k]], 'only.offset')) {
      return(NULL)
    }

    # Build column-block k of the beta x beta matrix
    ffunk <- function(j) {
      if (attr(frame$stage.dictionary[[j]], 'only.offset')) {
        return(NULL)
      }

      # build row-block j of the column-block k
      if (k > j) {
        l <- which((jk[,1] == k) & (jk[,2] == j))
        bjk <- matrix(mu.beta2.list[[l]],
                      nrow = n_j[j], ncol = n_j[k],
                      byrow = TRUE)
      }
      else {
        # where the (j,k)-matrix is stored
        l <- which((jk[,1] == j) & (jk[,2] == k))
        bjk <- matrix(mu.beta2.list[[l]],
                      nrow = n_j[j], ncol = n_j[k],
                      byrow = FALSE)
      }
      return(bjk)
    }
    dk <- lapply(1:q, FUN = ffunk)
    do.call('rbind', dk)
  }
  mu.beta2 <- lapply(1:q, FUN = ffun)
  mu.beta2 <- do.call('cbind', mu.beta2)

  if (d > 0) {
    mu.gamma2 <- matrix(colSums(replicate(d^2, c.score) * mu.gamma.gamma),
                        nrow = d, ncol = d)

    ffun <- function(mu.betaj.gamma) {
      if (is.null(mu.betaj.gamma))
        return(NULL)
      matrix(colSums(replicate(NCOL(mu.betaj.gamma), c.score) * mu.betaj.gamma),
             ncol = d)
    }
    mu.beta.gamma <- lapply(mu.betaj.gamma.list, FUN = ffun)
    mu.beta.gamma <- do.call('rbind', mu.beta.gamma)
  }
  else {
    mu.gamma2 <- mu.beta.gamma <- NULL
  }

  if (r > 0) {
    mu.delta2 <- matrix(colSums(replicate(r^2, c.score) * mu.delta.delta),
                        nrow = r, ncol = r)

    ffun <- function(mu.betaj.delta) {
      if (is.null(mu.betaj.delta))
        return(NULL)
      matrix(colSums(replicate(NCOL(mu.betaj.delta), c.score) * mu.betaj.delta),
             ncol = r)
    }
    mu.beta.delta <- lapply(mu.betaj.delta.list, FUN = ffun)
    mu.beta.delta <- do.call('rbind', mu.beta.delta)

    if (d > 0) {
      mu.gamma.delta <- matrix(colSums(replicate(r * d, c.score) * mu.gamma.delta.list),
                               nrow = d, ncol = r)
    }
    else {
      mu.gamma.delta <- NULL
    }
  }
  else {
    mu.delta2 <- mu.beta.delta <- mu.gamma.delta <- NULL
  }

  OFinfo <- rbind(cbind(             mu.beta2,                        mu.beta.gamma,   mu.beta.delta),
                 cbind(if (d > 0) t(mu.beta.gamma),                  mu.gamma2,       mu.gamma.delta),
                 cbind(if (r > 0) t(mu.beta.delta), if (d * r > 0) t(mu.gamma.delta), mu.delta2))

  c.Oinfo <- frame$y / (mu ^ 2) + (frame$weights - frame$y) / ((1 - mu)^2)
  c.Oinfo <- frame$sample.weights * c.Oinfo
  OFinfo <- t(out$mu.theta) %*% (replicate(NCOL(out$mu.theta), c.Oinfo) * out$mu.theta) -
    OFinfo

  out$OFinfo <- OFinfo
  dimnames(out$OFinfo) <- dimnames(out$EFinfo)

  out$dtau <- Taudots
  out$dscore <- list(mu.beta.beta.list = mu.beta.beta.list,
                     mu.betaj.gamma.list = mu.betaj.gamma.list,
                     mu.betaj.delta.list = mu.betaj.delta.list,
                     mu.gamma.gamma = mu.gamma.gamma,
                     mu.delta.delta = mu.delta.delta,
                     mu.gamma.delta.list = mu.gamma.delta.list)

  if (!jerk) {
    return(out)
  }

  #* Compute expected jerk matrix (inflection, change in curvature/acceleration)
  # Get mu.theta.theta
  ffun <- function (i) {
    # Obtain the matrix ddot{\mu}_i
    # beta x beta block
    ffuni <- function (l) {
      j <- jk[l,1]
      k <- jk[l,2]

      if (is.null(mu.beta.beta.list[[l]])) {
        return(NULL)
      }

      hjk <- matrix(mu.beta.beta.list[[l]][i,],
                    nrow = n_j[j], ncol = n_j[k])
      return(hjk)
    }
    mu.beta2.list.i <- lapply(1:NROW(jk), FUN = ffuni)

    ffunl <- function (k) {
      if (attr(frame$stage.dictionary[[k]], 'only.offset')) {
        return(NULL)
      }

      # Build column-block k of the beta x beta matrix
      ffunk <- function(j) {
        if (attr(frame$stage.dictionary[[j]], 'only.offset')) {
          return(NULL)
        }

        # build row-block j of the column-block k
        if (k > j) {
          l <- which((jk[,1] == k) & (jk[,2] == j))
          bjk <- matrix(mu.beta2.list.i[[l]],
                        nrow = n_j[j], ncol = n_j[k],
                        byrow = TRUE)
        }
        else {
          # where the (j,k)-matrix is stored
          l <- which((jk[,1] == j) & (jk[,2] == k))
          bjk <- matrix(mu.beta2.list.i[[l]],
                        nrow = n_j[j], ncol = n_j[k],
                        byrow = FALSE)
        }
        return(bjk)
      }
      dk <- lapply(1:q, FUN = ffunk)
      do.call('rbind', dk)
    }
    mu.beta2.i <- lapply(1:q, FUN = ffunl)
    mu.beta2.i <- do.call('cbind', mu.beta2.i)

    # Alpha block
    if (d > 0) {
      mu.gamma2.i <- matrix(mu.gamma.gamma[i,], nrow = d, ncol = d)

      fffun <- function(mu.betaj.gamma) {
        if (is.null(mu.betaj.gamma))
          return(NULL)
        matrix(mu.betaj.gamma[i,], ncol = d)
      }
      mu.beta.gamma.i <- lapply(mu.betaj.gamma.list, FUN = fffun)
      mu.beta.gamma.i <- do.call('rbind', mu.beta.gamma.i)
    }
    else {
      mu.gamma2.i <- mu.beta.gamma.i <- NULL
    }

    # Lambda block
    if (r > 0) {
      mu.delta2.i <- matrix(mu.delta.delta[i,], nrow = r, ncol = r)

      fffun <- function(mu.betaj.delta) {
        if (is.null(mu.betaj.delta))
          return(NULL)
        matrix(mu.betaj.delta[i,], ncol = r)
      }
      mu.beta.delta.i <- lapply(mu.betaj.delta.list, FUN = fffun)
      mu.beta.delta.i <- do.call('rbind', mu.beta.delta.i)

      if (d > 0) {
        mu.gamma.delta.i <- matrix(mu.gamma.delta.list[i,],
                                   nrow = d, ncol = r)
      }
      else {
        mu.gamma.delta.i <- NULL
      }
    }
    else {
      mu.delta2.i <- mu.beta.delta.i <- mu.gamma.delta.i <- NULL
    }

    mu.theta.theta.i <- rbind(cbind(             mu.beta2.i,                        mu.beta.gamma.i,   mu.beta.delta.i),
                              cbind(if (d > 0) t(mu.beta.gamma.i),                  mu.gamma2.i,       mu.gamma.delta.i),
                              cbind(if (r > 0) t(mu.beta.delta.i), if (d * r > 0) t(mu.gamma.delta.i), mu.delta2.i))

    return(mu.theta.theta.i)

  }
  mu.theta.theta.list <- lapply(1:frame$dims$nobs, FUN = ffun)
  out$mu.theta.theta.list <- mu.theta.theta.list

  # Get dot{Einfo}
  vweights <- rep(frame$weights * frame$sample.weights, length.out = frame$dims$nobs)
  ffun <- function(i) {
    mu_i <- mu[i]
    mu.theta.i <- mu.theta[i,]
    mu.theta.theta.i <- mu.theta.theta.list[[i]]

    outi <- if (frame$dims$npars == 1) 1 else matrixcalcK.matrix(r = frame$dims$npars, c = frame$dims$npars)
    outi <- outi + diag(frame$dims$npars^2)
    outi <- c(outi %*% dirprod(mu.theta.theta.i, mu.theta.i))
    outi <- outi + ((2*mu_i - 1) / (mu_i * (1 - mu_i))) * c(dirprod(mu.theta.i, mu.theta.i, t(mu.theta.i)))
    outi <- (vweights[i] / (mu_i * (1 - mu_i))) * outi

    return(outi)
  }
  Einfo.theta <- sapply(1:frame$dims$nobs, FUN = ffun)
  if (frame$dims$npars == 1)
    out$Einfo.theta <- sum(Einfo.theta)
  else
    out$Einfo.theta <- matrix(rowSums(Einfo.theta), nrow = frame$dims$npars^2,
                              ncol = frame$dims$npars, byrow = FALSE)

  return(out)
}






