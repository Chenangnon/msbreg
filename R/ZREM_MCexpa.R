
#
# Requires that all elements of 'betas' have the same length
#
#' @importFrom utils write.table
#' @importFrom utils str
# ... = additional arguments passed to MCmsbm
MCexpa <- function (ns = c(100, 200, 300, 500, 750, 1000, 2000),
                    lambdas = c(0.25, 0.5, 0.75, 0.9, 1),
                    betas = list(c(3, -1)),
                    xdistrs = list('unif', 'norm', 'lnorm'),
                    link = "logit",
                    msbm.control = list(),
                    criterions = "MLJ",
                    unit.lambda = TRUE,
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
                       kxdistr = 1:length(xdistrs))

  Nsettings <- NROW(pgrid)
  cat(paste0('\nTotal settings: ', Nsettings, '.\n'))
  onesetting <- function (jgr) {

    cat(paste0(jgr, '.'))

    nlambdajbetakxdistr <- pgrid[jgr,]
    n <- as.numeric(nlambdajbetakxdistr[1])
    lambda <- as.numeric(nlambdajbetakxdistr[2])
    j <- as.numeric(nlambdajbetakxdistr[3])
    k <- as.numeric(nlambdajbetakxdistr[4])

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

    # simulate data
    mcdata <- as.data.frame (
      sim.predictors (n = n, xdistr = xdistr,
                      xmin = wmin, xmax = wmax)
    )

    mcformulax <- paste0(colnames(mcdata), collapse = " | ")
    mcformula <- as.formula(paste0("y ~ ", mcformulax))
    mcformulax <- as.formula(paste0("~ ", mcformulax))
    eta.lambda <- switch(link,
                         logit = qlogis(lambda),
                         qnorm(lambda))

    mcframe <- msbm.frame(mcformulax,
                          lambda.formula = ~ 1,
                          data = mcdata)
    theta <- c(beta, delta = eta.lambda)


    mcdata$y <- simulate(mcframe, nsim = 1, link = link, theta = theta)
    mcdata$y <- as.numeric(mcdata$y)

    # Define the true model
    thetastart <- theta
    if (is.infinite(eta.lambda))
      thetastart <- c(beta, delta = 36)
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
        km <- km + 1
        fitmodels[[km]] <- msbreg (mcformula,
                                   lambda.formula = ~ 0,
                                   data = mcdata, frame = TRUE,
                                   start = beta,
                                   control = msbm.control,
                                   criterion = criterions[kmj])
        fitmodels[[km]]$call$criterion <- criterions[kmj]
      }

      km <- km + 1
      if (identical(criterions[kmj], "MLJ")) {
        fitmodels[[km]] <- mcmsbmfit
      }
      else {
        fitmodels[[km]] <- msbreg (mcformula,
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
                method = method, ...)
      })$value
    })

    if (any(class(Braw) %in% c("simpleError", "error",
                                          "condition", "try-error"))) {
      return(list(sumry = NA, Braw = NA))
    }

    # Save results
    Bresults <- list(n = rep(n, NROW(Braw$estimates)),
                     lambda = lambda, xdistr = xdistr)
    Bresults <- c(Bresults, as.list(beta))
    Bresults <- list(Braw = do.call("cbind", Bresults))

    Bresults$Braw <- as.data.frame(Bresults$Braw)
    Bresults$Braw$n <- as.numeric(Bresults$Braw$n)
    Bresults$Braw$lambda <- as.numeric(Bresults$Braw$lambda)
    for (clj in 1:nbbeta) {
      Bresults$Braw[, 3 + clj] <- as.numeric(Bresults$Braw[, 3 + clj])
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
    Bresults$sumry <- cbind(Bresults$Braw[1, 1:(3 + nbbeta)],
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

  Braw <- as.data.frame(do.call('rbind', Braw))

  if (!is.null(savepath)) {
    savepath <- paste0(savepath, 'batchrun_')
    write.table(Braw, file  = paste0(savepath, 'Braw.txt'))
    write.table(sumry, file  = paste0(savepath, 'sumry.txt'))
  }

  return(list(sumry = sumry, Braw = Braw))

}

sim.local <- function (statistics, nsim, stat.names = NULL, ...) {

  if (!is.function(statistics))
    statistics <- get(statistics, mode = "function",
                      envir = parent.frame())

  outmat <- lapply(1:nsim, FUN = statistics)
  nsimcol <- sapply(outmat, FUN = length)
  if (!all(nsimcol == max(nsimcol))) {
    maxsim <- which.max(nsimcol)[1]
    outmat <- c(outmat[[maxsim]], outmat[-maxsim])
  }
  outmat <- do.call("rbind", outmat)

  if (NCOL(outmat) == length(stat.names))
    colnames(outmat) <- stat.names

  list(estimates = as.matrix(outmat), raw = NULL)
}
