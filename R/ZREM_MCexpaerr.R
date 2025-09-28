# Requires that all elements of 'betas' have the same length
#
#' @importFrom utils write.table

MCexpaerr <- function (ns = c(100, 200, 300, 500, 750, 1000, 1500, 2000),
                       lambdas = c(0.25, 0.375, 0.5, 0.625, 0.75,
                                   0.8,  0.85,  0.9, 0.95,  1),
                       betas = list(c(3, -1), c(3, -2)),
                       xdistrs = list('unif', 'norm', 'lnorm'),
                       edistrs = list('bridge', 'norm'),
                       esigmas = list(0.1, 0.5, 1),
                       link = "logit",
                       msbm.control = list(),
                       criterions = "MLJ",
                       unit.lambda = TRUE, # Include fit with lambda = 1?
                       unit.all = TRUE, # FALSE ==> Only fit assuming lambda = 1 for no.error and unit.error
                       no.error = TRUE,    # Include fit with input.me = list()?
                       half.error = TRUE,  # Include fit with input.me = list(errx / 2)?
                       double.error = TRUE,# Include fit with input.me = list(errx * 2)?
                       B = 2000, seed = NULL,
                       cl = NULL, chunk.size = NULL,
                       savepath = NULL,
                       basic.stats = "all",
                       sup.stats = function(x) if (x$dims$r == 0) extractscores(score.test(x)),
                       sup.args = NULL, fisher.matrix = TRUE,
                       adjust = FALSE, sandwich = FALSE,
                       R2.method = c("KL", "Nagelkerke", "COR"),
                       cor.method = "pearson", method = "sim.local", ...) {
  # Take unique elements of simulation parameters
  ns <- unique(ns)
  lambdas <- unique(lambdas)
  betas <- unique(betas)
  xdistrs <- unique(xdistrs)
  esigmas <- unique(esigmas)
  stopifnot(all(unlist(esigmas) >= 0))

  # Ensure all betas have the same length
  lbeta <- sapply(betas, length)
  if (length(lbeta) > 1)
    stopifnot(all(lbeta[-1] == lbeta[1]))

  #* Save/Set the random generator seed
  eval(setsave.RNGstate())

  #* Control argument
  msbm.control <- do.call( "msbm.control", msbm.control)
  msbm.control <- msbm.control[names(msbm.control) != 'criterion']
  if (length(criterions) == 0) {
    criterions <- "MLJ"
  }
  else {
    stopifnot(all(criterions %in% c("ML", "MLJ", "MLJIC")))
  }
  criterions <- unique(criterions)

  pgrid <- expand.grid(n = ns,
                       lambda = lambdas,
                       jbeta = 1:length(betas),
                       kxdistr = 1:length(xdistrs),
                       ledistr = 1:length(edistrs),
                       mesigma = 1:length(esigmas))

  Nsettings <- NROW(pgrid)
  cat(paste0('\nTotal settings: ', Nsettings, '.\n'))
  onesetting <- function (jgr) {

    cat(paste0(jgr, '.'))

    nlambdajbetakxdistr <- pgrid[jgr,]
    n <- as.numeric(nlambdajbetakxdistr[1])
    lambda <- as.numeric(nlambdajbetakxdistr[2])
    j <- as.numeric(nlambdajbetakxdistr[3])
    k <- as.numeric(nlambdajbetakxdistr[4])
    l <- as.numeric(nlambdajbetakxdistr[5])
    m <- as.numeric(nlambdajbetakxdistr[6])

    beta <- betas[[j]]
    nbbeta <- length(beta)
    intcpt <- seq(1, nbbeta, by = 2)
    slopes <- seq(2, nbbeta, by = 2)
    if (is.null(names(beta))) {
      names(beta) <- c(rbind(paste0("beta0", 1:length(intcpt)),
                             paste0("beta", 1:length(slopes))))
    }

    xdistr <- xdistrs[[k]]
    wslopes <- beta[slopes]
    wslopes[wslopes == 0] <- -1 # Assume beta[2] == -1
    wmin <- (5 - beta[intcpt]) / wslopes
    wmax <- (-5 - beta[intcpt]) / wslopes

    edistr <- edistrs[[l]]
    esigma <- esigmas[[m]]

    # simulate data
    mcdata <- as.data.frame (
      sim.predictors (n = n, xdistr = xdistr,
                      xmin = wmin, xmax = wmax)
    )
    wnms <- colnames(mcdata)
    nvars <- NCOL(mcdata)

    mcformulax <- paste0(wnms, collapse = " | ")
    mcformula <- as.formula(paste0("y ~ ", mcformulax))
    mcformulax <- as.formula(paste0("~ ", mcformulax))
    eta.lambda <- switch(link,
                         logit = qlogis(lambda),
                         qnorm(lambda))

    mcframe <- msbm.frame(mcformulax,
                          lambda.formula = ~ 1,
                          data = mcdata)
    theta <- c(beta, delta = eta.lambda)

    mcdataerr <- sim.xerrors (x = mcdata, edistr = edistr, esigma = esigma)
    xnms <- paste0(wnms, "err")
    colnames(mcdataerr) <- xnms

    mcformula_err <- paste0(colnames(mcdataerr), collapse = " | ")
    mcformula_err <- as.formula(paste0("y ~ ", mcformula_err))

    sdmat <- matrix(esigma, nrow = n, ncol = nvars)
    sdnm <- paste0("err", wnms)
    colnames(sdmat) <- sdnm
    input.me <- vector(mode = "list", length = nvars)
    for (j in 1:nvars) {
      input.me[[j]] <- as.formula(paste0("~", sdnm[j]))
    }
    names(input.me) <- xnms

    mcdata <- cbind(mcdata, mcdataerr, sdmat)
    mcdata$y <- simulate(mcframe, nsim = 1, link = link, theta = theta)
    mcdata$y <- as.numeric(mcdata$y)

    if (half.error) {
      sdmat <- matrix(esigma/2, nrow = n, ncol = nvars)
      sdnm <- paste0("herr", wnms)
      colnames(sdmat) <- sdnm
      hinput.me <- vector(mode = "list", length = nvars)
      for (j in 1:nvars) {
        hinput.me[[j]] <- as.formula(paste0("~", sdnm[j]))
      }
      names(hinput.me) <- xnms
      mcdata <- cbind(mcdata, sdmat)
    }

    if (double.error) {
      sdmat <- matrix(esigma * 2, nrow = n, ncol = nvars)
      sdnm <- paste0("derr", wnms)
      colnames(sdmat) <- sdnm
      dinput.me <- vector(mode = "list", length = nvars)
      for (j in 1:nvars) {
        dinput.me[[j]] <- as.formula(paste0("~", sdnm[j]))
      }
      names(dinput.me) <- xnms
      mcdata <- cbind(mcdata, sdmat)
    }

    # Define the true model
    thetastart <- as.numeric(theta)
    if (is.infinite(eta.lambda))
      thetastart <- c(as.numeric(beta), delta = 36)
    mcmsbmfit <- msbreg (mcformula,
                         lambda.formula = ~ 1,
                         data = mcdata, frame = TRUE,
                         start = thetastart,
                         control = msbm.control,
                         criterion = 'MLJ')

    # Define fitting models
    fitmodels <- list()
    km <- 0
    for (kmj in 1:length(criterions)) {
      if (unit.lambda) {
        if (no.error) {
          km <- km + 1
          fitmodels[[km]] <- msbreg (mcformula_err,
                                     lambda.formula = ~ 0,
                                     data = mcdata, frame = TRUE,
                                     start = as.numeric(beta),
                                     control = msbm.control,
                                     criterion = criterions[kmj])
          fitmodels[[km]]$call$criterion <- criterions[kmj]
        }

        if (half.error & unit.all) {
          km <- km + 1
          fitmodels[[km]] <- msbreg (mcformula_err,
                                     input.me = hinput.me,
                                     lambda.formula = ~ 0,
                                     data = mcdata, frame = TRUE,
                                     start = as.numeric(beta),
                                     control = msbm.control,
                                     criterion = criterions[kmj])
          fitmodels[[km]]$call$criterion <- criterions[kmj]
        }

        km <- km + 1
        fitmodels[[km]] <- msbreg (mcformula_err,
                                   input.me = input.me,
                                   lambda.formula = ~ 0,
                                   data = mcdata, frame = TRUE,
                                   start = as.numeric(beta),
                                   control = msbm.control,
                                   criterion = criterions[kmj])
        fitmodels[[km]]$call$criterion <- criterions[kmj]

        if (double.error & unit.all) {
          km <- km + 1
          fitmodels[[km]] <- msbreg (mcformula_err,
                                     input.me = dinput.me,
                                     lambda.formula = ~ 0,
                                     data = mcdata, frame = TRUE,
                                     start = as.numeric(beta),
                                     control = msbm.control,
                                     criterion = criterions[kmj])
          fitmodels[[km]]$call$criterion <- criterions[kmj]
        }
      }

      if (no.error) {
        km <- km + 1
        fitmodels[[km]] <- msbreg (mcformula_err,
                                   lambda.formula = ~ 1,
                                   data = mcdata, frame = TRUE,
                                   start = thetastart,
                                   control = msbm.control,
                                   criterion = criterions[kmj])
        fitmodels[[km]]$call$criterion <- criterions[kmj]
      }

      if (half.error) {
        km <- km + 1
        fitmodels[[km]] <- msbreg (mcformula_err,
                                   input.me = hinput.me,
                                   lambda.formula = ~ 1,
                                   data = mcdata, frame = TRUE,
                                   start = thetastart,
                                   control = msbm.control,
                                   criterion = criterions[kmj])
        fitmodels[[km]]$call$criterion <- criterions[kmj]
      }

      km <- km + 1
      fitmodels[[km]] <- msbreg (mcformula_err,
                                 input.me = input.me,
                                 lambda.formula = ~ 1,
                                 data = mcdata, frame = TRUE,
                                 start = thetastart,
                                 control = msbm.control,
                                 criterion = criterions[kmj])
      fitmodels[[km]]$call$criterion <- criterions[kmj]

      if (double.error) {
        km <- km + 1
        fitmodels[[km]] <- msbreg (mcformula_err,
                                   input.me = dinput.me,
                                   lambda.formula = ~ 1,
                                   data = mcdata, frame = TRUE,
                                   start = thetastart,
                                   control = msbm.control,
                                   criterion = criterions[kmj])
        fitmodels[[km]]$call$criterion <- criterions[kmj]
      }
    }

    # Run B MC simulations
    Btime <- system.time({
      Braw <- catch.conditions({
        MCmsbm (msbm = mcmsbmfit, theta = theta,
                models = fitmodels, nsim = B,
                basic.stats = basic.stats, sup.stats = sup.stats,
                sup.args = sup.args, fisher.matrix = fisher.matrix,
                adjust = adjust, sandwich = sandwich,
                R2.method = R2.method, cor.method = cor.method,
                method = method)
      })$value
    })

    if (any(class(Braw) %in% c("simpleError", "error",
                               "condition", "try-error"))) {
      return(list(sumry = NA, Braw = NA))
    }

    # Save results
    Bresults <- list(n = rep(n, NROW(Braw$estimates)),
                     lambda = lambda, xdistr = xdistr,
                     edistr = edistr, esigma = esigma)
    Bresults <- c(Bresults, as.list(beta))
    Bresults <- list(Braw = do.call("cbind", Bresults))

    Bresults$Braw <- as.data.frame(Bresults$Braw)
    Bresults$Braw$n <- as.numeric(Bresults$Braw$n)
    Bresults$Braw$lambda <- as.numeric(Bresults$Braw$lambda)
    Bresults$Braw$esigma <- as.numeric(Bresults$Braw$esigma)

    for (clj in 1:nbbeta) {
      Bresults$Braw[, 5 + clj] <- as.numeric(Bresults$Braw[, 5 + clj])
    }
    Bresults$Braw <- cbind(Bresults$Braw, Braw$estimates)

    cat("\n Braw$estimates \n ")
    print(str(Braw$estimates))
    print(colnames(Braw$estimates))

    Bresults$sumry <- c(colMeans(Braw$estimates, na.rm = TRUE),
                        colVars(Braw$estimates, na.rm = TRUE),
                        Elapsed = Btime[3])
    Bresults$sumry <- as.data.frame(rbind(Bresults$sumry))
    colnames (Bresults$sumry) <- c(paste0("Mean.", colnames(Braw$estimates)),
                                   paste0("Var.", colnames(Braw$estimates)),
                                   "Elapsed")
    Bresults$sumry <- cbind(Bresults$Braw[1, 1:(5 + nbbeta)],
                            Bresults$sumry)

    return(Bresults)
  }

  Bouts <- matteLapply (X = 1:Nsettings,
                        FUN = onesetting,
                        cl = cl, chunk.size = chunk.size)
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }

  sumry <- lapply(Bouts, FUN = function(x) {
    if (is.list(x))
      x$sumry
    else
      NA
  })

  nsimcol <- sapply(sumry, FUN = length)
  if (!all(nsimcol == max(nsimcol))) {
    maxsim <- which.max(nsimcol)[1]
    sumry <- c(sumry[[maxsim]], sumry[-maxsim])
  }

  sumrynms <- names(sumry[[1]])
  sumry <- lapply(sumry, FUN = function(x) {
    names(x) <- sumrynms[1:length(x)]
    return(x)
  })

  sumry <- as.data.frame(do.call('rbind', sumry))

  Braw <- lapply(Bouts, FUN = function(x) {
    if (is.list(x))
      x$Braw
    else
      NA
  })

  nsimcol <- sapply(Braw, FUN = length)
  if (!all(nsimcol == max(nsimcol))) {
    maxsim <- which.max(nsimcol)[1]
    Braw <- c(Braw[[maxsim]], Braw[-maxsim])
  }

  Brawnms <- colnames(Braw[[1]])
  Braw <- lapply(Braw, FUN = function(x) {
    colnames(x) <- Brawnms[1:NCOL(x)]
    return(x)
  })

  Braw <- as.data.frame(do.call('rbind', Braw))

  if (!is.null(savepath)) {
    savepath <- paste0(savepath, 'batchrun_')
    write.table(Braw, file  = paste0(savepath, 'Braw.txt'))
    write.table(sumry, file  = paste0(savepath, 'sumry.txt'))
  }

  return(list(sumry = sumry, Braw = Braw))

}

sim.xerrors <- function (x, edistr = "bdridge", esigma = 0.1) {
  n <- NROW(x)
  p <- NCOL(x)

  if (length(edistr) == 1 & length(esigma) == 1) {
    switch(edistr,
           bdridge = {
             ephi <- 1 / sqrt(3 * ((esigma^2)/pi^2) + 1)
             errmat <- matrix (rbridge(n * p, location = 0, scale = ephi),
                               nrow = n, ncol = p)
           },
           {
             errmat <- matrix (rnorm(n * p, mean = 0, sd = esigma),
                               nrow = n, ncol = p)
           })

  }
  else {
    edistr <- rep(edistr, length.out = p)
    esigma <- rep(esigma, length.out = p)

    errmat <- sapply(1:p, FUN = function (j) {
      switch(edistr[j],
             bdridge = {
               ephij <- 1 / sqrt(3 * ((esigma[j]^2)/pi^2) + 1)
               errmatj <- rbridge(n, location = 0, scale = ephij)
             },
             {
               errmatj <- rnorm(n, mean = 0, sd = esigma[j])
             })
      errmatj
    })

  }

  out <- x + errmat

  cnm <- colnames(x)
  if (length(cnm) == p)
    colnames(out) <- paste0(cnm, "err")

  out
}
