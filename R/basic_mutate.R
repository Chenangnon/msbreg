
basic.mutate <- function () {
  expression({
    #* Check that frame is an 'msbm.frame' object
    stopifnot(is.msbm.frame(frame))

    #* Check parameter dimension
    stopifnot(length(theta) == length(frame$parnames))

    #* Set dimensions from mbm.frame: q, pj (for j = 1:q), d, r
    eval(get.msbm.dims())

    #* Extract sub-parameter vectors
    # Regression parameters
    endpj <- cumsum(pj + intercepts)
    startpj <- c(0, endpj[-q]) + 1
    betas <- lapply(1:q, FUN = function(j) {
      theta[startpj[j]:endpj[j]]
    })

    # 'alpha' parameters
    gamma <- if (d > 0) theta[(endpj[q] + 1):(endpj[q] + d)]

    # 'lambda' parameters
    delta <- if (r > 0) theta[(endpj[q] + d + 1):(endpj[q] + d + r)]

    #* Roll over design matrices and compute linear predictors
    ffun <- function(j) {
      dict <- frame$stage.dictionary[[j]]

      offsetj <- if (attr(dict, 'offset'))
        frame$offset.matrix[, attr(dict, 'offset.index')]
      else 0

      if (attr(dict, 'only.offset')) {
        return(offsetj)
      }
      betaj <- betas[[j]]

      if (pj[j] > 0) {
        X <- frame$input.matrix[, dict, drop = FALSE]

        if (attr(dict, 'intercept')) {
          etaj <- betaj[1] + c(X %*% betaj[-1])
        }
        else {
          etaj <- c(X %*% betaj)
        }
      }
      else {
        etaj <- rep(betaj, length.out = nobs)
      }

      etaj <- etaj + offsetj

      return(etaj)
    }

    if (q == 1) {
      etas <- cbind(ffun(1))
    }
    else {
      etas <- sapply(1:q, FUN = ffun)
    }
    names(betas) <- colnames(etas) <- names(frame$stage.dictionary)

    #* Standard deviations
    if (frame$me | frame$me.offset) {
      ffun <- function(j) {
        dict <- frame$stage.dictionary[[j]]

        if (!attr(dict, 'has.me')) {

          if (!attr(dict, 'has.me.offset')) {
            return(numeric(nobs))
          }

          sd.offj <- frame$sd.offset.matrix[, attr(dict, 'me.offset.index'),
                                            drop = TRUE]

          return(sd.offj)
        }

        betaj <- betas[[j]]
        if (attr(dict, 'intercept')) {
          betaj <- betaj[-1]
        }
        betaj <- betaj[attr(dict, 'contr.with.me')]
        sdX <- frame$sd.input.matrix[, attr(dict, 'me.index'), drop = FALSE]
        varj <- c((sdX^2) %*% (betaj^2))

        if (attr(dict, 'has.me.offset')) {
          varj <- varj +
            frame$sd.offset.matrix[, attr(dict, 'me.offset.index'), drop = TRUE]^2
        }

        return(sqrt(varj))
      }

      if (q == 1) {
        SDs <- cbind(ffun(1))
      }
      else {
        SDs <- sapply(1:q, FUN = ffun)
      }
      colnames(SDs) <- names(frame$stage.dictionary)
    }
    else {
      SDs <- 0
    }

    #* Inverse link functions
    linkinv <- link$linkinv
    melinkinv <- link$melinkinv

    #* Compute the alpha component
    if (d == 0) {
      if (attr(frame$alpha.dictionary, 'offset')) {
        eta.alpha <- frame$offset.matrix[, attr(frame$alpha.dictionary,
                                                'offset.index')]
        logalpha <- linkinv(eta.alpha, log.p = TRUE)
      }
      else {
        eta.alpha <- -Inf
        logalpha <- -Inf
      }
    }
    else {
      if (attr(frame$alpha.dictionary, 'intercept')) {
        if (d == 1) {
          eta.alpha <- gamma
        }
        else {
          eta.alpha <-  gamma[1] + c(frame$input.matrix[, frame$alpha.dictionary, drop = FALSE] %*% gamma[-1])
        }
      }
      else {
        eta.alpha <- c(frame$input.matrix[, frame$alpha.dictionary, drop = FALSE] %*% gamma)
      }

      if (attr(frame$alpha.dictionary, 'offset'))
        eta.alpha <- frame$offset.matrix[, attr(frame$alpha.dictionary,
                                                'offset.index')] + eta.alpha
      logalpha <- linkinv(eta.alpha, log.p = TRUE)
    }

    #* Compute the lambda component
    if (r == 0) {
      if (attr(frame$lambda.dictionary, 'offset')) {
        eta.lambda <- frame$offset.matrix[, attr(frame$lambda.dictionary,
                                                 'offset.index')]
        loglambda <- linkinv(eta.lambda, log.p = TRUE)
      }
      else {
        eta.lambda <- Inf
        loglambda <- 0
      }
    }
    else {
      if (attr(frame$lambda.dictionary, 'intercept')) {
        if (r == 1) {
          eta.lambda <- delta
        }
        else {
          eta.lambda <- delta[1] + c(frame$input.matrix[, frame$lambda.dictionary, drop = FALSE] %*% delta[-1])
        }
      }
      else {
        eta.lambda <- c(frame$input.matrix[, frame$lambda.dictionary, drop = FALSE] %*% delta)
      }

      if (attr(frame$lambda.dictionary, 'offset'))
        eta.lambda <- frame$offset.matrix[, attr(frame$lambda.dictionary,
                                                 'offset.index')] + eta.lambda

      loglambda <- linkinv(eta.lambda, log.p = TRUE)
    }

    #* Compute the vector 'mu' of success probabilities
    logmus <- melinkinv (etas, sd = SDs, log.p = TRUE)
    rowSumlogmus <- rowSums(logmus)
    mu <- exp(loglambda + logalpha) + (1 - exp(logalpha)) * exp(loglambda + rowSumlogmus)
    validmu <- is.finite(mu) & (mu >= 0 & mu <= 1)

  })
}
