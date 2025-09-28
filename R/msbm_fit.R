# @rdname msbreg
#' @export msbm.fit
# @importFrom stats binomial
#'
#* Write "msbm.fit" which will use 'BB::BBsolve' to solve the score
# equation. Or use the alternative 'nleqslv::nleqslv' or pracma::fsolve
# The above two lines will become necessary when Firth's
# modified score will be implemented.
msbm.fit <- function (frame, y = frame$y,
                      start = NULL,
                      link = "logit",
                      control = list(),
                      singular.ok = TRUE) {
  stopifnot(length(y) == NROW(frame$input.matrix))
  control <- do.call("msbm.control", control)
  ynm <- names(y)
  frame$y <- as.numeric(y)
  names(frame$y) <- ynm

  if (length(frame$sample.weights) > 1) {
    frame$sample.weights <- frame$dims$nobs *
    frame$sample.weights / sum(frame$sample.weights)
  }

  #* Get Binomial link function
  eval(get.linkfun())

  #* Initialize the model parameter
  if (is.null(start)) {
    if (identical(control$start, "none") & is.null(control$pars)) {
      start <- initialize.msbm.fit(frame, link = link,
                                   control = control)
      start <- as.numeric(start)
      names (start) <- frame$parnames
    }
  }
  else {
    stopifnot(is.numeric(start))
    stopifnot(length(start) %in% c(frame$dims$npars - 1, frame$dims$npars, frame$dims$npars + 1))

    # Deal with simple updates of 'lambda.formula' with a non-updated 'start' argument
    if (length(start) == frame$dims$npars - 1) {
      if (frame$dims$r == 1) {
        start <- c(start, `(Intercept).lambda` = link$linkfun(.99))
      }
      else {
        stopifnot(length(start) == frame$dims$npars)
      }
    }
    else if (length(start) == frame$dims$npars + 1) {
      if (frame$dims$r == 0) {
        start <- start[-length(start)]
      }
      else {
        stopifnot(length(start) == frame$dims$npars)
      }
    }

    nmstart <- names(start)
    if (length(nmstart) != frame$dims$npars) {
      if (length(nmstart) > 0)
        names(start) <- c(nmstart, frame$parnames[-c(1:length(nmstart))])
      else
        names(start) <- frame$parnames
    }
  }

  #* (penalized-)deviance function to be minimized
  exp.lambda <- FALSE
  eval(toget_half_devfun())

  #* Add box constraint if necessary
  # Default limits
  lowerb <- -Inf
  upperb <- Inf

  # Indicators of regression slopes (in theta)
  #* Set dimensions from mbm.frame: q, pj (for j = 1:q), d, r
  eval(get.msbm.dims())
  weights <- rep(frame$weights, length.out = nobs)
  n.ok <- nobs - sum(weights * frame$sample.weights == 0)

  # Determine slopes with restricted ranges
  if (nb.slops) {
    # Replicate 'slope.signs' to the length of 'slope.id'
    control$slope.signs <- rep(control$slope.signs,
                               length.out = nb.slops)

    # Slopes restricted to the negative orphan
    negative.orphan <- control$slope.signs < 0
    if (any(c(negative.orphan))) {
      upperb <- rep(Inf, frame$dims$npars)
      upperb[slope.id[negative.orphan]] <- 0
      if (!is.null(start)) {
        start[slope.id[negative.orphan]] <- pmin(-control$epsilon,
                                                 start[slope.id[negative.orphan]])
      }
    }

    # Slopes restricted to the positive orphan
    positive.orphan <- control$slope.signs > 0
    if (any(c(positive.orphan))) {
      lowerb <- rep(-Inf, frame$dims$npars)
      lowerb[slope.id[positive.orphan]] <- 0
      if (!is.null(start)) {
        start[slope.id[positive.orphan]] <- pmax(control$epsilon,
                                                 start[slope.id[negative.orphan]])
      }
    }
  }

  #* Call 'optim' for deviance minimization
  if (control$trace) {
    if (identical(control$criterion, "ML"))
      cat("\n   ***  deviance minimization ... \n")
    else
      cat("\n   ***  penalized-deviance minimization ... \n")
  }

  # Search with batch starting values, if any
  starts <- NULL
  finalopt <- 1
  if (!identical(control$start, "none") | !is.null(control$pars)) {
    if (!is.null(control$pars)) {
      stopifnot (NCOL(control$pars) == frame$dims$npars)
    }

    if (control$trace) {
      cat("\n   ***  initial search over ", NROW(control$pars), " initial values  ... \n")
    }

    DevMinZation0 <- catch.conditions({
      optim.selfstart (frame$dims$npars,
                       design = control$start,
                       star.points = control$star.points,
                       inf = control$inf,
                       pars = control$pars,
                       fn = half_devfun,
                       method = control$method,
                       lower = lowerb,
                       upper = upperb,
                       control = list(trace = max(0, control$trace-1),
                                      maxit = control$maxit,
                                      reltol = control$epsilon),
                       hessian = FALSE)
    })$value

    # Warn if no solution
    if (any(class(DevMinZation0) %in% c("simpleError", "error",
                                        "condition", "try-error"))) {
      if (control$trace) {
        warning(paste0("batch fitting failled: ",
                       DevMinZation0))

      }

      # Initialize theta
      if (is.null(start)) {
        start <- initialize.msbm.fit(frame, link = link,
                                     control = control)
        start <- as.numeric(start)
        names (start) <- frame$parnames
        if (nb.slops) {
          if (any(c(negative.orphan))) {
            start[slope.id[negative.orphan]] <- pmin(-control$epsilon,
                                                     start[slope.id[negative.orphan]])
          }

          if (any(c(positive.orphan))) {
            start[slope.id[positive.orphan]] <- pmax(control$epsilon,
                                                     start[slope.id[negative.orphan]])
          }
        }
      }
    }
    else {
      starts <- DevMinZation0[c("table", "starts", "k")]

      start <- DevMinZation0$starts[DevMinZation0$k,]
      names (start) <- frame$parnames
    }
  }

  # Do a prior simplex search, if no batch of starting values was used/succeeded
  if (is.null(starts)) {

    if (control$trace &
        identical(control$start, "none") &
        is.null(control$pars)) {
      cat("\n     +  initial search ... \n\n")
    }

    DevMinZation0 <- catch.conditions({
      stats::optim (par = start,
                    fn = half_devfun,
                    method = "Nelder-Mead",
                    lower = lowerb,
                    upper = upperb,
                    control = list(trace = max(0, control$trace-1),
                                   maxit = max(25, control$maxit/2),
                                   reltol = control$epsilon),
                    hessian = FALSE)
    })$value


    # Try a "SANN" search, if "Nelder-Mead" failed
    # Do a prior warming/cooling search
    # BUT not sure if SANN is useful at all, the method depends heavily on the setting
    if (any(class(DevMinZation0) %in% c("simpleError", "error",
                                        "condition", "try-error"))) {
      if (!identical(control$method, "SANN")) {
        DevMinZation0 <- catch.conditions({
          stats::optim (par = start,
                        fn = half_devfun,
                        method = "SANN",
                        lower = lowerb,
                        upper = upperb,
                        control = list(trace = max(0, control$trace-1),
                                       maxit = min(max(1000, control$maxit/2), 10000),
                                       temp = 20, fnscale = 1e-2,
                                       reltol = control$epsilon),
                        hessian = FALSE)
        })$value
      }

      # Stop if no solution
      if (any(class(DevMinZation0) %in% c("simpleError", "error",
                                          "condition", "try-error"))) {

        finalopt <- 0
        DevMinZation0 <- list(par = start,
                              value = half_devfun (start),
                              counts = c(0, 0),
                              convergence = 1000,
                              message = "fitting failled: Reduce model size?",
                              hessian = NULL)

        if (control$trace) {
          warning(paste0("fitting failled: ",
                         DevMinZation0,
                         ". \nReduce model size, or try another optimizer (see '?msbm.control')"))
        }

      }
    }
  }
  start2 <- DevMinZation0$par
  if (!is.null(start)) {
    names(start2) <- names(start)
  }
  else {
    names(start2) <- frame$parnames
  }

  # Final optimization
  if (finalopt) {
    if (control$trace) {
      cat("\n     +  final optimization ... \n\n")
    }

    DevMinZation <- catch.conditions({
      stats::optim (par = start2,
                    fn = half_devfun,
                    method = control$method,
                    lower = lowerb,
                    upper = upperb,
                    control = list(trace = max(0, control$trace-1),
                                   maxit = control$maxit,
                                   fnscale = 1e-10,
                                   reltol = control$epsilon),
                    hessian = TRUE)
    })$value

    if (any(class(DevMinZation) %in% c("simpleError", "error", "condition",
                                       "try-error"))) {
      finalopt <- 0

      if (control$trace) {
        warning(paste0("fitting failled: ",
                       DevMinZation,
                       ". Try another optimizer (see '?msbm.control')"))
      }
    }
  }

  if (!finalopt) {
    DevMinZation <- DevMinZation0
    DevMinZation0$counts <- c(0, 0)

    DevMinZation$hessian <- catch.conditions({
      numDerivhessian (func = half_devfun, x = DevMinZation0$par)
    })$value

    if (any(class(DevMinZation$hessian) %in% c("simpleError", "error", "condition",
                                       "try-error"))) {
      DevMinZation$hessian <- matrix(NA, nrow = npars, ncol = npars)
    }
  }

  # Is UNIT asymptote search allowed or requested?
  unit.request <- FALSE
  if (is.null(control$unit.lambda)) {
    control$unit.lambda <- attr(frame$lambda.dictionary, "intercept") &
      !attr(frame$lambda.dictionary, "offset")
  }
  else {
    control$unit.lambda <- as.logical(control$unit.lambda)[1]
    if (!attr(frame$lambda.dictionary, "intercept") & control$unit.lambda) {
      if (control$trace) {
        warning(" ignoring 'control$unit.lambda' because there is no intercept in the linear predictor of 'lambda'. \n")
      }

      control$unit.lambda <- FALSE
    }
    unit.request <- control$unit.lambda
  }

  # Perform unit asymptote search
  unitresults <- NULL
  method.args <- list(eps=control$epsilon, d=0.0001,
                      zero.tol=sqrt(.Machine$double.eps/7e-7),
                      r=4, v=2, show.details=FALSE)
  if (control$unit.lambda) {
    if (control$trace) {
      cat("\n          ... exploring unit-asymptote parameter region ... \n\n")
    }

    # Initial parameter values
    unitstart0 <- start[1:(frame$dims$npars - frame$dims$r)]
    unitstart <- DevMinZation$par[1:(frame$dims$npars - frame$dims$r)]

    # Parameter bounds
    if (length(lowerb) > 1) {
      unitlowerb <- lowerb[1:(frame$dims$npars - frame$dims$r)]
    }
    else {
      unitlowerb <- lowerb
    }

    if (length(upperb) > 1) {
      unitupperb <- upperb[1:(frame$dims$npars - frame$dims$r)]
    }
    else {
      unitupperb <- upperb
    }

    # Fixed values of lambda link scale parameters under unit-asymptote assumption
    appendpars <- DevMinZation$par[-c(1:(frame$dims$npars - frame$dims$r))]
    appendpars[1] <- Inf

    #* (penalized-)deviance function to be minimized
    switch(incriterion,
           ML = {
             unithalf_devfun <- function (theta) {
               mutatedtheta <- mutate.params (c(theta, appendpars), frame = frame,
                                              link = link,
                                              score = FALSE,
                                              information = FALSE,
                                              observed = FALSE)
               neg_devvalue <- mutatedtheta$loglike

               return(-neg_devvalue)
             }
           },
           MLJ = {
             if (frame$dims$npars > 1) {
               unithalf_devfun <- function (theta) {
                 mutatedtheta <- mutate.params (c(theta, appendpars), frame = frame,
                                                link = link,
                                                direct.lambda = TRUE,
                                                information = TRUE,
                                                observed = FALSE)

                 penalty <- try({
                   determinant(mutatedtheta$EFinfo, logarithm = TRUE)
                 }, silent = TRUE)

                 if (is(penalty, 'try-error')) {
                   return(exp(308) * nobs)
                 }

                 if (penalty$sign < 0) {
                   return(exp(308) * nobs)
                 }

                 penalty <- penalty$modulus[1]

                 neg_devvalue <- mutatedtheta$loglike + 0.5 * penalty

                 return(-neg_devvalue)
               }
             }
             else {
               # Unlikely situation
               unithalf_devfun <- function (theta) {
                 mutatedtheta <- mutate.params (c(theta, appendpars), frame = frame,
                                                link = link,
                                                direct.lambda = TRUE,
                                                information = TRUE,
                                                observed = FALSE)

                 if (mutatedtheta$EFinfo[1] < 0) {
                   return(exp(308) * nobs)
                 }

                 neg_devvalue <- mutatedtheta$loglike +
                   0.5 * logb(mutatedtheta$EFinfo[1])

                 return(-neg_devvalue)
               }
             }
           },
           MLPJ = {
             slope.id.unit <- slope.id[slope.id <= (frame$dims$npars - frame$dims$r)]
             if (length(slope.id.unit) > 1) {
               unithalf_devfun <- function (theta) {
                 mutatedtheta <- mutate.params (c(theta, appendpars), frame = frame,
                                                link = link,
                                                direct.lambda = TRUE,
                                                information = TRUE,
                                                observed = FALSE)

                 penalty <- try({
                   determinant(mutatedtheta$EFinfo[slope.id.unit, slope.id.unit],
                               logarithm = TRUE)
                 }, silent = TRUE)

                 if (is(penalty, 'try-error')) {
                   return(exp(308) * nobs)
                 }

                 if (penalty$sign < 0) {
                   return(exp(308) * nobs)
                 }

                 penalty <- penalty$modulus[1]

                 neg_devvalue <- mutatedtheta$loglike + 0.5 * penalty

                 return(-neg_devvalue)
               }
             }
             else if (length(slope.id.unit) == 1) {
               unithalf_devfun <- function (theta) {
                 mutatedtheta <- mutate.params (c(theta, appendpars), frame = frame,
                                                link = link,
                                                direct.lambda = TRUE,
                                                information = TRUE,
                                                observed = FALSE)

                 if (mutatedtheta$EFinfo[slope.id.unit, slope.id.unit] < 0) {
                   return(exp(308) * nobs)
                 }

                 neg_devvalue <- mutatedtheta$loglike +
                   0.5 * logb(mutatedtheta$EFinfo[slope.id.unit, slope.id.unit])

                 return(-neg_devvalue)
               }
             }
             else {
               unithalf_devfun <- function (theta) {
                 mutatedtheta <- mutate.params (c(theta, appendpars), frame = frame,
                                                link = link,
                                                information = FALSE,
                                                observed = FALSE)

                 neg_devvalue <- mutatedtheta$loglike

                 return(-neg_devvalue)
               }
             }
           }
    )

    # Optimization with both starting values
    # From 'start' ...
    unitDevMinZation0 <- catch.conditions({
      stats::optim (par = unitstart0,
                    fn = unithalf_devfun,
                    method = "Nelder-Mead",
                    lower = unitlowerb,
                    upper = unitupperb,
                    control = list(trace = max(0, control$trace-1),
                                   maxit = max(25, control$maxit/2),
                                   reltol = control$epsilon),
                    hessian = TRUE)
    })$value

    if (any(class(unitDevMinZation0) %in% c("simpleError", "error", "condition",
                                           "try-error"))) {
      if (!identical(control$method, "Nelder-Mead")) {
        unitDevMinZation0 <- catch.conditions({
          stats::optim (par = unitstart0,
                        fn = unithalf_devfun,
                        method = control$method,
                        lower = unitlowerb,
                        upper = unitupperb,
                        control = list(trace = max(0, control$trace-1),
                                       maxit = control$maxit,
                                       fnscale = 1e-10,
                                       reltol = control$epsilon),
                        hessian = TRUE)
        })$value
      }
    }
    else {
      Dcounts0 <- unitDevMinZation0$counts
      unitDevMinZation0 <- catch.conditions({
        stats::optim (par = unitDevMinZation0$par,
                      fn = unithalf_devfun,
                      method = control$method,
                      lower = unitlowerb,
                      upper = unitupperb,
                      control = list(trace = max(0, control$trace-1),
                                     maxit = control$maxit,
                                     fnscale = 1e-10,
                                     reltol = control$epsilon),
                      hessian = TRUE)
      })$value

      if (!any(class(unitDevMinZation0) %in% c("simpleError", "error", "condition",
                                              "try-error"))) {
        unitDevMinZation0$counts <- unitDevMinZation0$counts + Dcounts0
      }
    }

    # From the best value so far
    unitDevMinZation <- catch.conditions({
      stats::optim (par = unitstart,
                    fn = unithalf_devfun,
                    method = "Nelder-Mead",
                    lower = unitlowerb,
                    upper = unitupperb,
                    control = list(trace = max(0, control$trace-1),
                                   maxit = control$maxit,
                                   reltol = control$epsilon),
                    hessian = TRUE)
    })$value

    if (any(class(unitDevMinZation) %in% c("simpleError", "error", "condition",
                                            "try-error"))) {
      if (!identical(control$method, "Nelder-Mead")) {
        unitDevMinZation <- catch.conditions({
          stats::optim (par = unitstart,
                        fn = unithalf_devfun,
                        method = control$method,
                        lower = unitlowerb,
                        upper = unitupperb,
                        control = list(trace = max(0, control$trace-1),
                                       maxit = control$maxit,
                                       reltol = control$epsilon),
                        hessian = TRUE)
        })$value
      }
    }
    else {
      Dcounts0 <- unitDevMinZation$counts
      unitDevMinZation <- catch.conditions({
        stats::optim (par = unitDevMinZation$par,
                      fn = unithalf_devfun,
                      method = control$method,
                      lower = unitlowerb,
                      upper = unitupperb,
                      control = list(trace = max(0, control$trace-1),
                                     maxit = control$maxit,
                                     fnscale = 1e-10,
                                     reltol = control$epsilon),
                      hessian = TRUE)
      })$value

      if (!any(class(unitDevMinZation) %in% c("simpleError", "error", "condition",
                                               "try-error"))) {
        unitDevMinZation$counts <- unitDevMinZation$counts + Dcounts0
      }
    }

    # Any success?
    if (any(class(unitDevMinZation0) %in% c("simpleError", "error", "condition",
                                           "try-error")) &
        any(class(unitDevMinZation) %in% c("simpleError", "error", "condition",
                                           "try-error"))) {
      if (control$trace) {
        warning(paste0("search of a unit asymptote fit failled: ",
                       unitDevMinZation,
                       if (unit.request) ". Try a direct fit with 'lambda.formula = ~ 0'"))
      }
    }
    else {
      # Pick the best solution
      if (any(class(unitDevMinZation) %in% c("simpleError", "error", "condition",
                                             "try-error"))) {
        unitDevMinZation <- unitDevMinZation0
      }
      else if (!any(class(unitDevMinZation0) %in% c("simpleError", "error", "condition",
                                                  "try-error"))) {
        if (unitDevMinZation0$value < unitDevMinZation$value) {
          unitDevMinZation <- unitDevMinZation0
        }
      }

      # Pick again the best with regard to the full model in 'DevMinZation'
      if (unitDevMinZation$value < DevMinZation$value) {
        DevMinZation$par[c(1:(frame$dims$npars - frame$dims$r))] <- unitDevMinZation$par
        DevMinZation$par[-c(1:(frame$dims$npars - frame$dims$r))] <- appendpars

        DevMinZation$counts <- DevMinZation$counts + unitDevMinZation$counts
        DevMinZation$convergence <- unitDevMinZation$convergence

        DevMinZation$value <- unitDevMinZation$value
        hpar <- DevMinZation$par
        hpar[frame$dims$npars - frame$dims$r + 1] <- exp(709) # Finite substitute to 'Inf'
        DevMinZation$hessian <- catch.conditions({
          numDerivhessian (func = half_devfun, x = hpar)
        })$value

      }
      unitresults <- unitDevMinZation[c('par', 'value', 'hessian', 'counts', 'convergence')]

      if (control$gradient) {
        unitresults$gradient <- c(numDerivjacobian(func = unithalf_devfun, x = unitresults$par,
                                                   method = "Richardson", method.args = method.args))
        names(unitresults$gradient) <- names(unitresults$par)
      }
    }
  }

  # Save convergence status
  iter <- DevMinZation$counts + DevMinZation0$counts
  converged <- DevMinZation$convergence == 0
  if (length(DevMinZation$convergence) == 0) {
    DevMinZation$convergence <- NA
    converged <- FALSE
    attr(converged, 'code') <- NA
  }
  else
    attr(converged, 'code') <- conv.optim(DevMinZation$convergence)
  attr(converged, 'message') <- DevMinZation$message
  boundary <- (DevMinZation$par == lowerb) | (DevMinZation$par == upperb)

  # Extracting basic results
  theta <- as.numeric(DevMinZation$par)
  names(theta) <- names(start)
  if (is.null(names(theta))) {
    names(theta) <- frame$parnames
  }
  criterion <- DevMinZation$value
  hessian <- DevMinZation$hessian
  if (control$gradient) {
    cgradient <- c(numDerivjacobian(func = half_devfun, x = theta,
                                    method = "Richardson", method.args = method.args))
    names(cgradient) <- names(theta)
  }
  else
    cgradient <- NULL

  # Apply intercept correction if required
  IC <- FALSE
  attr(IC, "IC.fail") <- FALSE
  IC_start <- NULL
  # Indicators of regression intercepts in theta
  int.id <- frame$parindex$intercepts
  if (identical(control$criterion, "MLJIC")) {
    if (control$trace) {
      cat("\n   ***  fitting for intercept correction ... \n")
    }

    # Indicators of linear predictors
    eta.id <- frame$parindex$eta.intercepts

    # Vector of intercepts only
    IC_start <- theta[int.id]
    if (any(is.infinite(IC_start)))
      IC_start[is.infinite(IC_start)] <- exp(709)

    # Reduce the model frame to estimate only intercepts
    eval(toget_reduced.frame())

    # Fitting for intercept correction
    eval(toget_ICfit())
  }

  #* Residual deviance function
  eval (toget_dev.resids())

  #* Compute various quantities to return
  eval (toget_fitquantities())

  #* Null deviance calculation
  if (control$trace) {
    cat("\n   ***  null deviance calculation ... \n\n")
  }
  eval (toget_null.msbm.fit())

  # Check if nulldev > dev
  if (!any(is.na(c(nulldev, dev)))) {
    if (nulldev < dev & control$trace) {
      warning(paste0("residual deviance of fitted model larger than null deviance: ",
                     "check model stage structure, and parameter bounds"))

      all(is.finite(nulltheta))
    }
  }

  # converged[1] <- converged[1] & (mrank == dims$npars)
  # Format output
  out <- list(coefficients = theta,
              fitted.values = mutated$mu,
              alpha.values  = alpha.values,
              lambda.values = lambda.values,
              rank = mrank,
              null.rank = nullrank,
              link = link,
              linear.predictors = mutated$etas,
              deviance = dev,
              aic = maic,
              bic = mbic,
              hqc = mhqc,
              deviance.resid = dev.ind,
              null.deviance = nulldev,
              iter = iter,
              nobs = frame$dims$nobs,
              weights = weights,
              sample.weights = frame$sample.weights,
              df.residual = resdf,
              df.null = nulldf,
              y = yback,
              y.resid = yresid,
              converged = converged,
              boundary = boundary,
              criterion = structure(criterion,
                                    gradient = cgradient,
                                    hessian = hessian), # objective = half_devfun
              fisher = structure(mutated$EFinfo,
                                 qr = QR,
                                 observed = structure(mutated$OFinfo,
                                                      qr = obsQR),
                                 logLik = mutated$loglike,
                                 score   = mutated$score,
                                 score.i = mutated$score.i,
                                 mu.coefficients = mutated$mu.theta,
                                 c.alpha = mutated$c.alpha,
                                 c.lambda = mutated$c.lambda),
              singular = singular,
              start = structure(start, IC = IC_start),
              starts = starts,
              null.theta = nulltheta,
              IC = IC,
              me = frame$me,
              me.offset = frame$me.offset,
              me.sd = if (frame$me | frame$me.offset) mutated$SDs,
              unit.asymptote = if(control$unit.lambda) unitresults)

  return(out)
}


conv.optim <- function(code) {
  paste(code, ": ",
        switch(as.character(code),
               `1` = "the iteration limit 'maxit' had been reached",
               `10` = "degeneracy of the Nelder-Mead simplex",
               `51` = "warning from the 'L-BFGS-B' method",
               `52` = "error from the 'L-BFGS-B' method",
               `1000` = "optimization failled: returning initial parameter values. Reduce model size?",
               "successful completion"))
}

#    unitframe <- frame
#    unitframe$dims$r <- 0
#    unitframe$lambda.dictionary <- numeric(0)
#    attr(unitframe$lambda.dictionary, 'intercept') <- 0
#    attr(unitframe$lambda.dictionary, 'covs.names') <- character(0)
#    attr(unitframe$lambda.dictionary, 'offset') <- FALSE
#
#    unitframe$parnames <- frame$parnames[1:(frame$dims$npars - frame$dims$r)]
#    unitframe$parindex$intercepts <- frame$parindex$intercepts[-1]
#    unitframe$parindex$eta.intercepts <- frame$parindex$eta.intercepts[-1]
