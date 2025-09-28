
toget_half_devfun <- function() {
  expression({
    # Slope indices
    slope.id <- frame$parindex$slopes
    nb.slops <- length(slope.id)

    #* (penalized-)deviance function to be minimized
    incriterion <- control$criterion
    if (identical(incriterion, 'MLJIC'))
      incriterion <- 'MLJ'
    switch(incriterion,
           ML = {
             half_devfun <- function (theta) {
               mutatedtheta <- mutate.params (theta, frame = frame,
                                              link = link,
                                              score = FALSE,
                                              information = FALSE,
                                              observed = FALSE,
                                              exp.lambda = exp.lambda)
               neg_devvalue <- mutatedtheta$loglike

               return(-neg_devvalue)
             }
           },
           MLJ = {
             if (frame$dims$npars > 1) {
               half_devfun <- function (theta) {
                 mutatedtheta <- mutate.params (theta, frame = frame,
                                                link = link,
                                                information = TRUE,
                                                observed = FALSE,
                                                exp.lambda = exp.lambda)

                 penalty <- try({
                   determinant(mutatedtheta$EFinfo, logarithm = TRUE)
                 }, silent = TRUE)

                 if (is(penalty, 'try-error')) {
                   return(exp(308) * frame$dims$nobs)
                 }

                 if (penalty$sign < 0) {
                   return(exp(308) * frame$dims$nobs)
                 }

                 penalty <- penalty$modulus[1]

                 neg_devvalue <- mutatedtheta$loglike + 0.5 * penalty

                 return(-neg_devvalue)
               }
             }
             else {
               half_devfun <- function (theta) {
                 mutatedtheta <- mutate.params (theta, frame = frame,
                                                link = link,
                                                information = TRUE,
                                                observed = FALSE,
                                                exp.lambda = exp.lambda)

                 if (mutatedtheta$EFinfo[1] < 0) {
                   return(exp(308) * frame$dims$nobs)
                 }

                 neg_devvalue <- mutatedtheta$loglike +
                   0.5 * logb(mutatedtheta$EFinfo[1])

                 return(-neg_devvalue)
               }
             }
           },
           MLPJ = {
             if (nb.slops > 1) {
               half_devfun <- function (theta) {
                 mutatedtheta <- mutate.params (theta, frame = frame,
                                                link = link,
                                                information = TRUE,
                                                observed = FALSE,
                                                exp.lambda = exp.lambda)

                 penalty <- try({
                   determinant(mutatedtheta$EFinfo[slope.id, slope.id],
                               logarithm = TRUE)
                 }, silent = TRUE)

                 if (is(penalty, 'try-error')) {
                   return(exp(308) * frame$dims$nobs)
                 }

                 if (penalty$sign < 0) {
                   return(exp(308) * frame$dims$nobs)
                 }

                 penalty <- penalty$modulus[1]

                 neg_devvalue <- mutatedtheta$loglike + 0.5 * penalty

                 return(-neg_devvalue)
               }
             }
             else if (nb.slops == 1) {
               half_devfun <- function (theta) {
                 mutatedtheta <- mutate.params (theta, frame = frame,
                                                link = link,
                                                information = TRUE,
                                                observed = FALSE,
                                                exp.lambda = exp.lambda)

                 if (mutatedtheta$EFinfo[slope.id, slope.id] < 0) {
                   return(exp(308) * frame$dims$nobs)
                 }

                 neg_devvalue <- mutatedtheta$loglike +
                   0.5 * logb(mutatedtheta$EFinfo[slope.id, slope.id])

                 return(-neg_devvalue)
               }
             }
             else {
               # Unlikely situation
               half_devfun <- function (theta) {
                 mutatedtheta <- mutate.params (theta, frame = frame,
                                                link = link,
                                                information = FALSE,
                                                observed = FALSE,
                                                exp.lambda = exp.lambda)

                 neg_devvalue <- mutatedtheta$loglike

                 return(-neg_devvalue)
               }
             }
           },
           {
             if (control$trace) {
               warning(paste0("currently recognized fitting criteria include:",
                              " 'ML' (Maximum Likelihood)",
                              ", 'MLJ' (ML with Jeffreys' priors) ",
                              ", and 'MLJIC' (MLJ with intercept correction)."))
             }
             stop(paste0("fitting criterion ", control$criterion, "not recognized."))
           }
    )
  })
}

toget_reduced.frame <- function() {
  expression({
    # Compute linear predictors (offsets for stage in the IC frame)
    mutated <- mutate.params (theta = theta,
                              frame = frame,
                              link = link,
                              score = FALSE,
                              information = FALSE,
                              observed = FALSE)

    # Build the IC frame
    IC_frame <- frame
    IC_frame$input.matrix <- NULL
    IC_frame$offset.matrix <- mutated$etas
    IC_frame$offset.matrix[, eta.id] <- t(t(mutated$etas[, eta.id, drop = FALSE]) - theta[int.id])
    IC_frame$me <- FALSE
    IC_frame$sd.input.matrix <- NULL

    ffun_IC <- function(j) {
      stg <- numeric(0L)
      atrbsj <- attributes(frame$stage.dictionary[[j]])
      atrbsj$names <- NULL
      attributes(stg) <- atrbsj
      attr(stg, "order") <- attr(stg, "assign") <- 0
      attr(stg, "offset") <- 1
      attr(stg, "offset.index") <- j
      attr(stg, "only.offset") <- !attr(stg, "intercept")
      attr(stg, "has.me") <- FALSE
      attr(stg, "contr.with.me") <- attr(stg, "me.index") <-
        attr(stg, "me.offset.index") <- numeric(0)
      attr(stg, "has.me.offset") <- attr(frame$stage.dictionary[[j]], "has.me") |
        attr(frame$stage.dictionary[[j]], "has.me.offset")

      return(stg)
    }
    IC_frame$stage.dictionary <- lapply(1:frame$dims$q,
                                        FUN = ffun_IC)

    stghasme <- sapply(IC_frame$stage.dictionary,
                       function(stg) {
                         attr(stg, "has.me.offset")
                       })
    IC_frame$me.offset <- frame$me | frame$me.offset
    if (IC_frame$me.offset)
      IC_frame$sd.offset.matrix <- mutated$SDs[, which(stghasme)]
    csumstghasme <- cumsum(stghasme)

    ffun_IC <- function(j) {
      stg <- IC_frame$stage.dictionary[[j]]
      if (attr(stg, "has.me.offset")) {
        attr(stg, "me.offset.index") <- csumstghasme[j]
      }
      return(stg)
    }
    IC_frame$stage.dictionary <- lapply(1:frame$dims$q,
                                        FUN = ffun_IC)
    if (length(IC_frame$alpha.dictionary)) {
      alattr <- attributes(IC_frame$alpha.dictionary)
      alattr$names <- NULL
      IC_frame$alpha.dictionary <- numeric(0)
      attributes(IC_frame$alpha.dictionary) <- alattr
    }
    attr(IC_frame$alpha.dictionary,"covs.names") <- character(0)
    attr(IC_frame$alpha.dictionary, "offset") <- TRUE
    attr(IC_frame$alpha.dictionary, "offset.index") <- q + 1

    if (length(IC_frame$lambda.dictionary)) {
      laattr <- attributes(IC_frame$lambda.dictionary)
      laattr$names <- NULL
      IC_frame$lambda.dictionary <- numeric(0)
      attributes(IC_frame$lambda.dictionary) <- laattr
    }
    attr(IC_frame$lambda.dictionary,"covs.names") <- character(0)
    attr(IC_frame$lambda.dictionary, "offset") <- TRUE
    attr(IC_frame$lambda.dictionary, "offset.index") <- q + 2
    IC_frame$stage.terms <- NULL
    IC_frame$auxy <- NULL
    IC_frame$dims$p <- length(IC_start)
    IC_frame$dims$pj <- rep(0, q)
    IC_frame$dims$d <- attr(IC_frame$alpha.dictionary, "intercept") + 0
    IC_frame$dims$r <- attr(IC_frame$lambda.dictionary, "intercept") + 0
    IC_frame$dims$pt <- IC_frame$dims$tj <-
      IC_frame$nb_covs <- IC_frame$nb_me.cols <- 0
    IC_frame$parnames <- names(IC_start)
    IC_frame$tables <- NULL
    IC_frame$frames <- NULL
    IC_frame$call0 <- IC_frame$call
    IC_frame$call <- NULL

    class(IC_frame) <- "msbm.frame"
  })
}

toget_ICfit <- function(){
  expression({
    # Deviance function to be minimized
    IC_half_devfun <- function (theta) {
      mutatedtheta <- mutate.params (theta, frame = IC_frame,
                                     link = link,
                                     score = FALSE,
                                     information = FALSE,
                                     observed = FALSE)

      neg_devvalue <- mutatedtheta$loglike

      return(-neg_devvalue)
    }

    # Call 'optim' for deviance minimization
    DevMinZation <- catch.conditions({
      stats::optim (par = IC_start,
                    fn = IC_half_devfun,
                    method = control$method,
                    lower = -Inf,
                    upper = Inf,
                    control = list(trace = FALSE,
                                   maxit = control$maxit,
                                   reltol = control$epsilon),
                    hessian = TRUE)
    })$value

    if (any(class(DevMinZation) %in% c("simpleError", "error", "condition", "try-error")) &
        !identical(control$method, "Nelder-Mead")) {
      DevMinZation <- catch.conditions({
        stats::optim (par = IC_start,
                      fn = IC_half_devfun,
                      method = "Nelder-Mead",
                      lower = -Inf,
                      upper = Inf,
                      control = list(trace = FALSE,
                                     maxit = control$maxit,
                                     reltol = control$epsilon),
                      hessian = TRUE)
      })$value
    }

    if (any(class(DevMinZation) %in% c("simpleError", "error", "condition", "try-error"))) {
      attr(IC, "IC.fail") <- TRUE
      attr(IC, "IC.fail.message") <- paste0("fitting for intercept correction failled: ",
                                            DevMinZation, ".")
      if (control$trace) {
        warning(attr(IC, "IC.fail.message"))
      }
    }
    else {
      IC <- TRUE
      attr(theta, 'before.IC') <- theta
      theta[int.id] <- DevMinZation$par

      attr(DevMinZation$counts, 'before.IC') <- iter
      iter <- DevMinZation$counts


      attr(converged, 'before.IC') <- converged
      converged[1] <- DevMinZation$convergence == 0

      attr(converged, 'code') <- conv.optim(DevMinZation$convergence)
      attr(converged, 'message') <- DevMinZation$message

      attr(criterion, "before.IC") <- criterion
      criterion[1] <- DevMinZation$value # half_devfun(theta)

      bef.cjacobian <- cjacobian
      # cjacobian <- c(numDerivjacobian(func = IC_half_devfun, x = theta[int.id], method = "Richardson")) #
      cjacobian[1:frame$dims$npars] <- c(numDerivjacobian(func = half_devfun, x = theta, method = "Richardson",
                                                          method.args = method.args))
      names(cjacobian) <- names(theta) # names(theta[int.id]) #
      attr(cjacobian, "before.IC") <- bef.cjacobian

      attr(DevMinZation$hessian, "before.IC") <- hessian
      #hessian <- DevMinZation$hessian #
      hessian <- numDerivhessian(func = half_devfun, x = theta, method = "Richardson") #
      attr(hessian, "before.IC") <- attr(DevMinZation$hessian, "before.IC")
      dimnames(hessian) <- list(names(theta), names(theta))
    }
  })
}

toget_dev.resids <- function() {
  expression({
    # Saturated model log-likelihood for each observation
    mumaxlik <- frame$y / weights
    maxlik.i <- numeric(nobs)
    upnext <- !is.na(mumaxlik) & (weights > 1)
    if (any(upnext)) {
      maxlik.i[upnext] <- stats::dbinom (x = frame$y[upnext],
                                         size = weights[upnext],
                                         prob = mumaxlik[upnext],
                                         log = TRUE)
    }

    # Model deviance
    dev.resids <- function (y, mu, size) {
      2 * (maxlik.i - stats::dbinom (x = y, size = size,
                                     prob = mu, log = TRUE))
    }
  })
}

toget_fitquantities <- function() {
  expression({
    # Compute logLik, score and information

    mutated <- mutate.params (theta = theta,
                              frame = frame,
                              link = link,
                              score = TRUE,
                              information = TRUE,
                              observed = TRUE)

    # Fitted model rank
    if (any(is.na(mutated$EFinfo))) {
      mrank <- dims$npars
      QR <- NA
      attr(QR, "solve") <- NA
      singular <- FALSE
    }
    else if (any(is.infinite(mutated$EFinfo))) {
      naid <- rowSums(mutated$EFinfo) + colSums(mutated$EFinfo)
      naid <- is.na(naid) | is.infinite(naid)
      mrank <- sum(!naid)
      QR <- NA
      attr(QR, "solve") <- NA
      singular <- mrank < frame$dims$npars
    }
    else {
      QR <- catch.conditions({
        qr(mutated$EFinfo) # chol would be faster, but requires a strictly positive definite matrix
      })$value

      if (any(class(QR) %in% c("simpleError", "error", "condition", "try-error"))) {
        QR <- NA
        if (any(is.na(mutated$EFinfo)) | any(is.infinite(mutated$EFinfo))) {
          naid <- rowSums(mutated$EFinfo) + colSums(mutated$EFinfo)
          naid <- is.na(naid) | is.infinite(naid)
          mrank <- sum(!naid)
        }
        else {
          mrank <- svd (mutated$EFinfo)$d
          mrank <- sum(abs(mrank) > control$epsilon)
        }
        singular <- TRUE
      }
      else {

        mrank <- QR$rank
        singular <- mrank < frame$dims$npars
        if (singular) {
          if (!singular.ok) {
            stop("singular fit encountered: consider droping any aliased coefficient in each model component")
          }

          if (control$trace)
            warning("singular Fisher information matrix")

          attr(QR, "solve") <- NA
        }
        else {
          attr(QR, "solve") <- catch.conditions({
            solve.qr(QR)
          })$value

          if (any(class(attr(QR, "solve")) %in% c("simpleError", "error", "condition", "try-error"))) {
            attr(QR, "solve") <- NA
          }
        }
      }
    }

    if (any(is.na(mutated$OFinfo))) {
      obsQR <- NA
      attr(obsQR, "solve") <- NA
    }
    else {
      obsQR <- catch.conditions({
        qr(mutated$OFinfo) # chol would be faster, but requires a strictly positive definite matrix
      })$value

      if (any(class(obsQR) %in% c("simpleError", "error", "condition", "try-error"))) {
        obsQR <- NA
        attr(obsQR, "solve") <- NA
      }
      else {
        if (is.na(mrank))
          mrank <- obsQR$rank

        if (obsQR$rank < dims$npars) {
          attr(obsQR, "solve") <- NA
        }
        else {
          attr(obsQR, "solve") <- catch.conditions({
            solve.qr(obsQR)
          })$value

          if (any(class(attr(obsQR, "solve")) %in% c("simpleError", "error", "condition", "try-error"))) {
            attr(obsQR, "solve") <- NA
          }
        }
      }
    }

    #validmu <- is.finite(mutated$mu) & (mutated$mu > 0 & mutated$mu < 1)
    resdf <- n.ok - mrank

    # Residual deviance
    dev.ind <- dev.resids (y = frame$y, mu = mutated$mu, size = weights) *
      frame$sample.weights
    dev <- sum(dev.ind)
    mbic <- - 2 * mutated$loglike
    mhqc <- 2 * frame$dims$npars
    maic <- mbic + mhqc
    mhqc <- mbic + mhqc * log(log(frame$dims$nobs))
    mbic <- mbic + frame$dims$npars * logb(frame$dims$nobs)
    yback <- as.numeric(frame$y)
    names(yback) <- names(frame$y)
    yresid <- yback - as.numeric(mutated$mu * weights)

    # Limit probabilities
    if ((frame$dims$d - attr(frame$alpha.dictionary, "intercept")) == 0 & !attr(frame$alpha.dictionary, "offset")) {

      if (mutated$etas[1,q+1] == -Inf & link$link %in% c('probit', 'logit')) {
        alpha.values <- 0
      }
      else {
        alpha.values <- link$linkinv (mutated$etas[1,q+1])
        if (link$link %in% c('probit', 'logit'))
          alpha.values[mutated$etas[1,q+1] == -Inf] <- 0
      }
    }
    else {
      if (all(mutated$etas[,q+1] == -Inf) & link$link %in% c('probit', 'logit')) {
        alpha.values <- numeric(length(mutated$etas[,q+1]))
      }
      else {
        alpha.values <- link$linkinv (mutated$etas[,q+1])
        if (link$link %in% c('probit', 'logit'))
          alpha.values[mutated$etas[1,q+1] == -Inf] <- 0
      }
    }
    if ((frame$dims$r - attr(frame$lambda.dictionary, "intercept")) == 0 & !attr(frame$lambda.dictionary, "offset")) {
      lambda.values <- link$linkinv (mutated$etas[1,q+2])
    }
    else {
      lambda.values <- link$linkinv (mutated$etas[,q+2])
    }
  })
}

toget_null.msbm.fit <- function(){
  expression({
    #* Build Half-deviance function (up to a constant)
    nulldev <- NA
    Do_a_nullfit <- FALSE
    multiple_fits <- FALSE
    nulltheta <- numeric(npars)
    names(nulltheta) <- names(start)
    cumpj <- cumsum(pj + intercepts)
    insert.int <- numeric(0)

    # alpha and lambda offsets
    if (attr(frame$alpha.dictionary, "offset")) {
      alpha.offset <- frame$offset.matrix[, attr(frame$alpha.dictionary, "offset.index")]

      if (attr(frame$alpha.dictionary, "intercept")) {
        nulltheta[frame$parindex$intercepts][sum(intercepts) + 1] <- 0
      }
    }
    else {
      alpha.offset <- - Inf
      if (attr(frame$alpha.dictionary, "intercept")) {
        nulltheta[frame$parindex$intercepts][sum(intercepts) + 1] <- -Inf
      }
    }
    log.alpha_o <- link$linkinv (alpha.offset, log.p = TRUE)

    if (attr(frame$lambda.dictionary, "offset")) {
      lambda.offset <- frame$offset.matrix[, attr(frame$lambda.dictionary, "offset.index")]
      if (attr(frame$lambda.dictionary, "intercept")) {
        nulltheta[frame$parindex$intercepts][sum(intercepts) + attr(frame$alpha.dictionary, "intercept") + 1] <- 0
      }
    }
    else {
      lambda.offset <- Inf
      if (attr(frame$lambda.dictionary, "intercept")) {
        nulltheta[frame$parindex$intercepts][sum(intercepts) + attr(frame$alpha.dictionary, "intercept") + 1] <- Inf
      }
    }
    log.lambda_o <- link$linkinv (lambda.offset, log.p = TRUE)
    min.p <- exp(log.alpha_o + log.lambda_o)
    slope.p <- (1 - exp(log.alpha_o))

    # Stage offsets, if any:  [stg1 stg2 ... stgq]
    all.offsets <- sapply(frame$stage.dictionary,
                          FUN = function(stg) { attr(stg, "offset") }) == 1
    fofun <- function(stg) {
      if (attr(stg, "offset")) {
        if (attr(stg, "has.me.offset")) {
          me.den <- frame$sd.offset.matrix[, attr(stg, "me.offset.index")]
          me.den <- sqrt(1 + (me.den/link$sigma)^2)
          frame$offset.matrix[, attr(stg, "offset.index")] / me.den
        }
        else
          frame$offset.matrix[, attr(stg, "offset.index")]
      }
      else if (attr(stg, "intercept")) {
        0
      }
      else {
        Inf
      }

    }
    stage.offsets <- lapply(frame$stage.dictionary, FUN = fofun)

    # Intercepts, if any (alpha and lambda component excluded)
    all.intercepts <- intercepts == 1
    stghasintercept <- which(all.intercepts)

    # No intercept in any mu_j?
    if (!any(all.intercepts)) {
      # No offset as well?
      if (!any(all.offsets)) {
        log.OffP <- link$linkinv (0, log.p = TRUE)

        # Slopes are undefined.
        # It is possible to find a least square solution (slopes) for
        # each stage, with zero as the response (the 0 used in log.OffP)
        # But the fitted eta_j will rarely be zero for all observations.
        # So we consider undefined slopes.
        nulltheta[frame$parindex$slopes[frame$parindex$slopes <= cumpj[q]]] <- NaN
      }
      else if (q == 1) {
        log.OffP <- link$linkinv (stage.offsets[[1]], log.p = TRUE)
      }
      else {
        log.OffP <- lapply(stage.offsets,
                           FUN = link$linkinv, log.p = TRUE)
        log.OffP <- do.call("cbind", log.OffP)
        log.OffP <- rowSums(log.OffP)

        for (joff in 1:q) {
          if (!all.offsets[joff]) {
            goNaN <- frame$parindex$slopes <= cumpj[joff]
            if (joff > 1)
              goNaN <- goNaN & frame$parindex$slopes > cumpj[joff - 1]

            nulltheta[frame$parindex$slopes[goNaN]] <- NaN
          }
        }
      }

      nullmu <- min.p + slope.p * exp(log.lambda_o + log.OffP)
      nulldev <- sum(dev.resids (frame$y, nullmu, weights) *
                       frame$sample.weights)
      nullrank <- 0
    }
    else {
      fitdev.resids <- function (y, mu, size) {
        - sum(stats::dbinom (x = y, size = size,
                             prob = mu, log = TRUE))
      }

      if (sum(!intercepts) > 1) {
        log.OffP <- lapply(stage.offsets[!all.intercepts],
                           FUN = link$linkinv, log.p = TRUE)
        log.OffP <- do.call("cbind", log.OffP)
        log.OffP <- rowSums(log.OffP)

        if (any(!all.intercepts & !all.offsets)) {
          for (joff in 1:q) {
            if (!all.intercepts[joff] & !all.offsets[joff]) {
              goNaN <- frame$parindex$slopes <= cumpj[joff]
              if (joff > 1)
                goNaN <- goNaN & frame$parindex$slopes > cumpj[joff - 1]

              nulltheta[frame$parindex$slopes[goNaN]] <- NaN
            }
          }
        }
      }
      else if (sum(!intercepts) == 1) {
        log.OffP <- link$linkinv (stage.offsets[!all.intercepts][[1]], log.p = TRUE)
        joff <- which(!intercepts)
        if (!all.offsets[joff]) {
          goNaN <- frame$parindex$slopes <= cumpj[joff]
          if (joff > 1)
            goNaN <- goNaN & frame$parindex$slopes > cumpj[joff - 1]

          nulltheta[frame$parindex$slopes[goNaN]] <- NaN
        }
      }
      else
        log.OffP <- 0
      slope.p <- (1 - exp(log.alpha_o)) * exp (log.lambda_o + log.OffP)

      all.int.off <- all.intercepts & all.offsets
      if (length(stghasintercept) == 1) {
        insert.int <- frame$parindex$intercepts[1]
        activestg.offset <- stage.offsets[[stghasintercept]]
        nullmu_fun <- function(beta) {
          log.activestg <- link$linkinv (activestg.offset + beta, log.p = TRUE)
          nullmu <- min.p + slope.p * exp(log.activestg)
          nullmu
        }

        null_devfun <- function(beta) {
          nullmu <- nullmu_fun (beta)
          nulldev <- sum(fitdev.resids (frame$y, nullmu, weights) *
                           frame$sample.weights)
          return(nulldev)
        }

        Do_a_nullfit <- TRUE
        startnull <- head(theta[frame$parindex$intercepts], 1)
        nullrank <- 1
      }
      else if (!any(all.int.off)) {
        # Stages with intercepts have no offset, so their product is
        # equivalent when an intercept is included in any of them
        insert.int <- frame$parindex$intercepts[1]
        attr(insert.int, "switch") <- list(current = insert.int,
                                           candidate = frame$parindex$intercepts[-1])
        attr(insert.int, "switch")$candidate <-
          attr(insert.int, "switch")$candidate[attr(insert.int, "switch")$candidate <= cumpj[q]]

        nulltheta[frame$parindex$intercepts[-1][frame$parindex$intercepts[-1] <= cumpj[q]]] <- Inf

        slope.p <- (1 - exp(log.alpha_o)) * exp (log.lambda_o + log.OffP)

        nullmu_fun <- function(beta) {
          log.activestg <- link$linkinv (beta, log.p = TRUE)
          nullmu <- min.p + slope.p * exp(log.activestg)
          nullmu
        }

        null_devfun <- function(beta) {
          nullmu <- nullmu_fun (beta)
          nulldev <- sum(fitdev.resids (frame$y, nullmu, weights) *
                           frame$sample.weights)
          return(nulldev)
        }

        Do_a_nullfit <- TRUE
        startnull <- 0
        nullrank <- 1
      }
      else {
        ffun <- function(k) {
          activestg.offsetk <- stage.offsets[[k]]

          notk <- stghasintercept[which(stghasintercept != k)]
          if (length(notk) == 1) {
            log.OffPk <- log.OffP + link$linkinv (stage.offsets[[notk]],
                                                  log.p = TRUE)
          }
          else if (length(notk) > 1) {
            log.OffPk <- lapply(stage.offsets[notk],
                                FUN = link$linkinv, log.p = TRUE)
            log.OffPk <- do.call("cbind", log.OffPk)
            log.OffPk <- rowSums(log.OffPk)
            log.OffPk <- log.OffP + log.OffPk
          }
          else
            log.OffPk <- log.OffP

          slope.p <- (1 - exp(log.alpha_o)) * exp (log.lambda_o + log.OffPk)

          nullmu_fun <- function(beta) {
            log.activestg <- link$linkinv (activestg.offsetk + beta, log.p = TRUE)
            nullmu <- min.p + slope.p * exp(log.activestg)
            nullmu
          }

          null_devfun <- function(beta) {
            nullmu <- nullmu_fun (beta)
            nulldev <- sum(fitdev.resids (frame$y, nullmu, weights) *
                             frame$sample.weights)
            return(nulldev)
          }

          startnull <- theta[frame$parindex$intercepts][sum(intercepts[1:k])]

          list(nullmu_fun = nullmu_fun,
               null_devfun = null_devfun,
               startnull = startnull)
        }
        null_fit.list <- lapply(stghasintercept, FUN = ffun)

        Do_a_nullfit <- TRUE
        multiple_fits <- TRUE
        nullrank <- 1
      }
    }

    if (Do_a_nullfit) {

      DoNullfit <- toget_DoNullFit()

      if(multiple_fits) {
        nstghasintercept <- length(stghasintercept)
        nullpars <- nulldevs <- numeric(length = nstghasintercept) + NA
        nullmus <- vector(mode = "list", length = nstghasintercept)

        for (k in 1:nstghasintercept) { # Loop over stages with intercepts
          nullmu_fun <- null_fit.list[[k]]$nullmu_fun
          null_devfun <- null_fit.list[[k]]$null_devfun
          startnull <- null_fit.list[[k]]$startnull
          eval(DoNullfit)
          nullpars[k] <- nullpar
          nulldevs[k] <- nulldev
          nullmus[[k]] <- nullmu
        }

        kmindev <- which.min(nulldevs)
        nulldev <- nulldevs[kmindev]
        nullmu <- nullmus[[kmindev]]
        nullpar <- nullpars[kmindev]

        # Position of 'nullpar' in nulltheta
        insert.int <- frame$parindex$intercepts[kmindev]
        if (nstghasintercept > 1) {
          notk <- stghasintercept[which(stghasintercept != kmindev)]
          hasnooffset <- which(!all.offsets)
          notkhasnooffset <- notk[notk %in% hasnooffset]
          if (length(notkhasnooffset)) {
            nulltheta[frame$parindex$intercepts[notkhasnooffset]] <- Inf
          }
        }
      }
      else {
        eval(DoNullfit)
      }

      nulltheta[insert.int] <- nullpar
      attr(nulltheta, "switch") <- attr(insert.int, "switch")
    }
    nulldf <- n.ok - nullrank
  })
}

toget_DoNullFit <- function () {
  expression({
    DevMinZation <- catch.conditions({
      stats::optim (par = startnull,
                    fn = null_devfun,
                    method = control$method,
                    lower = -Inf,
                    upper = Inf,
                    control = list(trace = FALSE,
                                   maxit = control$maxit,
                                   reltol = control$epsilon),
                    hessian = FALSE)
    })$value

    if (any(class(DevMinZation) %in% c("simpleError", "error", "condition", "try-error")) &
        !identical(control$method, "Nelder-Mead")) {
      DevMinZation <- catch.conditions({
        stats::optim (par = startnull,
                      fn = null_devfun,
                      method = "Nelder-Mead",
                      lower = -Inf,
                      upper = Inf,
                      control = list(trace = FALSE,
                                     maxit = control$maxit,
                                     reltol = control$epsilon),
                      hessian = FALSE)
      })$value
    }

    if (any(class(DevMinZation) %in% c("simpleError", "error", "condition", "try-error"))) {
      if (control$trace) {
        warning(paste0("fitting for null deviance calculation failled: ",
                       DevMinZation, "."))
      }
      nullpar <- nullmu <- nulldev <- NA
    }
    else {
      nullpar <- DevMinZation$par
      nullmu <- nullmu_fun (nullpar)
      nulldev <- sum(dev.resids (frame$y, nullmu, weights) *
                       frame$sample.weights)
    }
  })
}

toget_null.msbm.fit0 <- function(){
  expression({
    # Half-deviance function (up to a constant)
    nulldev <- NULL
    Do_a_nullfit <- FALSE
    if (q == 1 & !attr(frame$alpha.dictionary, "intercept") &
        !attr(frame$lambda.dictionary, "intercept")) {
      # If only one stage and both lambda and alpha have no intercept
      # Consider the unique stage as a standard binomial 'glm'
      # There are three cases then: intercept & offset; only intercept, and only offset
      lambda.offset <- if (attr(frame$lambda.dictionary, "offset"))
        frame$offset.matrix[, attr(frame$lambda.dictionary, "offset.index")]
      else
        Inf
      lambda.offset <- link$linkinv (lambda.offset, log.p = TRUE)

      alpha.offset <- if (attr(frame$alpha.dictionary, "offset"))
        frame$offset.matrix[, attr(frame$alpha.dictionary, "offset.index")]
      else
        - Inf
      alpha.offset <- link$linkinv (alpha.offset, log.p = TRUE)

      min.p <- exp(lambda.offset + alpha.offset)
      slope.p <- exp(lambda.offset) * (1 - exp(alpha.offset))

      if (attr(frame$stage.dictionary[[1]], "intercept")) {
        # Include ME of stage.offset if any
        if (attr(frame$stage.dictionary[[1]], "has.me.offset")) {
          me.den <- frame$sd.offset.matrix[, attr(frame$stage.dictionary[[1]], "me.offset.index")]
          me.den <- sqrt(1 + (me.den/link$sigma)^2)
        }
        else {
          me.den <- 1
        }

        if (attr(frame$stage.dictionary[[1]], "offset")) {
          stage.offset <- frame$offset.matrix[, attr(frame$stage.dictionary[[1]], "offset.index")]
        }
        else {
          stage.offset <- 0
        }

        fitdev.resids <- function (y, mu, size) {
          - sum(stats::dbinom (x = y, size = size,
                               prob = mu, log = TRUE))
        }

        nullmu_fun <- function(beta) {
          min.p + slope.p * link$linkinv ((stage.offset + beta) / me.den)
        }
        null_devfun <- function(beta) {
          nullmu <- nullmu_fun (beta)
          nulldev <- sum(fitdev.resids (frame$y, nullmu, weights) *
                           frame$sample.weights)
          return(nulldev)
        }

        Do_a_nullfit <- TRUE
        startnull <- theta[1]
        nullrank <- 1
      }
      else {
        # Only offsets (alpha, lambda, and the unique stage)
        if (attr(frame$stage.dictionary[[1]], "offset")) {
          stage.offset <- frame$offset.matrix[, attr(frame$stage.dictionary[[1]], "offset.index")]
        }
        else {
          if (attr(frame$lambda.dictionary, "offset") |
              attr(frame$alpha.dictionary, "offset")) {
            stage.offset <- Inf
          }
          else
            stage.offset <- 0
        }

        # Include ME of stage.offset if any
        if (attr(frame$stage.dictionary[[1]], "has.me.offset")) {
          me.den <- frame$sd.offset.matrix[, attr(frame$stage.dictionary[[1]], "me.offset.index")]
          me.den <- sqrt(1 + (me.den/link$sigma)^2)
        }
        else {
          me.den <- 1
        }

        nullmu <- min.p + slope.p * link$linkinv (stage.offset / me.den)
        nulldev <- sum(dev.resids (frame$y, nullmu, weights) *
                         frame$sample.weights)
        nullrank <- 0
      }
    }
    else if (attr(frame$lambda.dictionary, "intercept") &
             attr(frame$lambda.dictionary, "offset")) {
      # Null deviance fit: the null model ignore all components, except 'lambda'
      # If lambda is just an offset (e.g. lambda fixed to one) with no intercept,
      # the null model is just that offset.
      lambda.offset <- frame$offset.matrix[, attr(frame$lambda.dictionary, "offset.index")]
      fitdev.resids <- function (y, mu, size) {
        - sum(stats::dbinom (x = y, size = size,
                             prob = mu, log = TRUE))
      }
      nullmu_fun <- function(beta) {
        link$linkinv (lambda.offset + beta)
      }
      null_devfun <- function(beta) {
        nullmu <- nullmu_fun (beta)
        nulldev <- sum(fitdev.resids (frame$y, nullmu, weights) *
                         frame$sample.weights)
        return(nulldev)
      }

      Do_a_nullfit <- TRUE
      startnull <- theta[sum(pj + intercept) + d + 1]
      nullrank <- 1
    }
    else if (attr(frame$lambda.dictionary, "intercept")) {
      nullmu <- sum(frame$sample.weights * frame$y) / sum(frame$sample.weights * weights)
      nulldev <- sum(dev.resids (frame$y, nullmu, weights) *
                       frame$sample.weights)
      nullrank <- 1
    }
    else {
      if (attr(frame$lambda.dictionary, "offset")) {
        nullmu <- link$linkinv(frame$offset.matrix[, attr(frame$lambda.dictionary, "offset.index")])
      }
      else {
        nullmu <- link$linkinv(0)
      }
      nulldev <- sum(dev.resids (frame$y, nullmu, weights) *
                       frame$sample.weights)
      nullrank <- 0
    }

    if (Do_a_nullfit) {
      DevMinZation <- catch.conditions({
        stats::optim (par = startnull,
                      fn = null_devfun,
                      method = control$method,
                      lower = -Inf,
                      upper = Inf,
                      control = list(trace = FALSE,
                                     maxit = control$maxit,
                                     reltol = control$epsilon),
                      hessian = FALSE)
      })$value

      if (any(class(DevMinZation) %in% c("simpleError", "error", "condition", "try-error")) &
          !identical(control$method, "Nelder-Mead")) {
        DevMinZation <- catch.conditions({
          stats::optim (par = startnull,
                        fn = null_devfun,
                        method = "Nelder-Mead",
                        lower = -Inf,
                        upper = Inf,
                        control = list(trace = FALSE,
                                       maxit = control$maxit,
                                       reltol = control$epsilon),
                        hessian = FALSE)
        })$value
      }

      if (any(class(DevMinZation) %in% c("simpleError", "error", "condition", "try-error"))) {
        if (control$trace) {
          warning(paste0("fitting for null deviance calculation failled: ",
                         DevMinZation, "."))
        }
      }
      else {
        nullmu <- nullmu_fun (DevMinZation$par)
        nulldev <- sum(dev.resids (frame$y, nullmu, weights) *
                         frame$sample.weights)
      }
    }
    nulldf <- n.ok - nullrank
  })
}
