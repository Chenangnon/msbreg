
apply.IC.msbm <- function (object, ...) {
  if (object$IC) {
    return(object)
  }

  control <- object$control
  if (identical(control$criterion, "ML")) {
    object$IC <- TRUE
    return(object)
  }

  if (attr(object$IC, "IC.fail")) {
    warning(attr(IC, "IC.fail.message"))
    return(object)
  }

  #* Preliminary object required in the workspace
  frame <- model.frame(object)
  theta <- object$coefficients
  link <- object$link

  eval(get.msbm.dims())
  weights <- rep(frame$weights, length.out = nobs)
  n.ok <- nobs - sum(weights * frame$sample.weights == 0)

  iter <- object$iter
  converged <- object$converged
  criterion <- object$criterion
  hessian <- attr(criterion, "hessian")
  attr(criterion, "hessian") <- NULL

  #* Fitting for intercept correction
  if (control$trace) {
    cat("\n   ***  fitting for intercept correction ... \n")
  }

  # Indicators of regression intercepts in theta
  int.id <- frame$parindex$intercepts
  eta.id <- frame$parindex$eta.intercepts

  # Vector of intercepts only
  IC_start <- theta[int.id]

  # Reduce the model frame to estimate only intercepts
  eval(toget_reduced.frame())

  # Intercept correction
  eval(toget_ICfit())

  #* Residual deviance function
  eval (toget_dev.resids())

  #* Compute various quantities to return
  eval (toget_fitquantities())

  #* Update oitput
  object$coefficients <- theta
  object$fitted.values <- mutated$mu
  object$alpha.values  <- alpha.values
  object$lambda.values <- lambda.values
  object$rank <- mrank
  object$linear.predictors <- mutated$etas
  object$deviance <- dev
  object$aic <- maic
  object$bic <- mbic
  object$hic <- mhqc
  object$deviance.resid <- dev.ind
  object$iter <- iter
  object$y.resid <- yresid
  object$converged <- converged
  object$criterion <- structure(criterion,
                                hessian = hessian)
  object$fisher <- structure(mutated$EFinfo,
                             qr = QR,
                             observed = structure(mutated$OFinfo,
                                                  qr = qr(mutated$OFinfo)),
                             logLik = mutated$loglike,
                             score   = mutated$score,
                             score.i = mutated$score.i,
                             mu.coefficients = mutated$mu.theta,
                             c.alpha = mutated$c.alpha,
                             c.lambda = mutated$c.lambda)
  object$singular <- singular
  object$start <- structure(object$start, IC = IC_start)
  object$IC <- TRUE

  return(object)
}

cancel.IC.msbm <- function(object, ...) {
  if (!object$IC) {
    return(object)
  }

  frame <- model.frame(object)
  theta <- object$coefficients
  MLIC <- attr(object$start, "IC")
  if (length(MLIC) == length(frame$parindex$intercepts)) {

    eval(get.msbm.dims())
    weights <- rep(frame$weights, length.out = nobs)
    n.ok <- nobs - sum(weights * frame$sample.weights == 0)
    theta[frame$parindex$intercepts] <- MLIC

    #* Update fitting criterion terms
    object$IC <- FALSE
    attr(object$iter, "IC") <-
      attr(object$converged, 'IC') <-
      attr(object$criterion, "IC") <-
      attr(object$hessian, "IC") <- NULL
    link <- object$link

    #* Compute various quantities to return
    eval (toget_fitquantities())

    #* Update the fit
    object$coefficients <- theta
    object$fitted.values  <- mutated$mu
    object$alpha.values   <- alpha.values
    object$lambda.values  <- lambda.values
    object$rank  <- mrank
    object$linear.predictors  <- mutated$etas
    object$deviance  <- dev
    object$aic  <- maic
    object$bic  <- mbic
    object$hic  <- mhqc
    object$deviance.resid  <- dev.ind
    object$y.resid  <- yresid
    object$fisher  <- structure(mutated$EFinfo,
                                qr = QR,
                                observed = structure(mutated$OFinfo,
                                                     qr = qr(mutated$OFinfo)),
                                logLik = mutated$loglike,
                                score   = mutated$score,
                                score.i = mutated$score.i,
                                mu.coefficients = mutated$mu.theta,
                                c.alpha = mutated$c.alpha,
                                c.lambda = mutated$c.lambda)
    object$singular  <- singular
    attr(object$start, "IC")  <- NULL

    return(object)
  }
  else {
    stop("invalid 'msbm' object")
  }
}
