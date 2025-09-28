
#
# Requires that all elements of 'betas' have the same length
#
#' @importFrom utils write.table
#' @importFrom utils str
# ... = additional arguments passed to MCmsbm
MCexpaprime <- function (ns = c(100, 200, 300, 500, 750, 1000, 2000),
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
                         sandwich = FALSE,
                         R2.method = c("KL", "Nagelkerke", "COR"),
                         resample.x = FALSE,
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
    stopifnot(all(criterions %in% c("ML", "MLPJ", "MLJ", "MLJIC")))
  }
  criterions <- unique(criterions)
  stopifnot(link %in% c('logit', 'probit'))

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

    # simulate predictors
    mcdata0 <- as.data.frame (
      sim.predictors (n = n, xdistr = xdistr,
                      xmin = wmin, xmax = wmax)
    )
    colnames(mcdata0) <- 'x'

    eta.lambda <- switch(link,
                         logit = qlogis(lambda),
                         qnorm(lambda))
    mcframe <- msbm.frame(~ x,
                          lambda.formula = ~ 1,
                          data = mcdata0)
    ## Compute true success probabilities
    theta <- c(beta, delta = eta.lambda)
    mu.values <- mutate.params(theta = theta,
                               frame = mcframe, link = link,
                               score = FALSE, information = FALSE,
                               observed = FALSE)$mu

    thetastart <- theta
    if (is.infinite(eta.lambda))
      thetastart <- c(beta, delta = 36)

    ## Statistics from each fit
    statnames <- list()
    statnames$GLM <- c('beta01', 'beta1', 'se.beta01', 'se.beta1', 'mad.mu', 'rmse.mu', 'R2.mu',
                       'null.deviance', 'deviance', 'logLik', 'converged')

    # ML fits
    if ('ML' %in% criterions) {
      if (unit.lambda) {
        statnames$MLunit <- c(statnames$GLM[-11], 'criterion', 'converged')
      }

      statnames$ML <- c('beta01', 'beta1', 'delta', 'se.beta01', 'se.beta1', 'se.delta', 'mad.mu',
                        'rmse.mu', 'R2.mu', 'null.deviance', 'deviance', 'logLik', 'criterion', 'converged')
    }

    # MLJ fits
    if ('MLJ' %in% criterions) {
      if (unit.lambda) {
        statnames$MLJunit <- c(statnames$GLM[-11], 'criterion', 'converged')
      }

      statnames$MLJ <- c('beta01', 'beta1', 'delta', 'se.beta01', 'se.beta1', 'se.delta', 'mad.mu',
                         'rmse.mu', 'R2.mu', 'null.deviance', 'deviance', 'logLik', 'criterion', 'converged')
    }

    # MLPJ fits
    if ('MLPJ' %in% criterions) {
      if (unit.lambda) {
        statnames$MLPJunit <- c(statnames$GLM[-11], 'criterion', 'converged')
      }

      statnames$MLPJ <- c('beta01', 'beta1', 'delta', 'se.beta01', 'se.beta1', 'se.delta', 'mad.mu',
                         'rmse.mu', 'R2.mu', 'null.deviance', 'deviance', 'logLik', 'criterion', 'converged')
    }

    # MLJIC fits
    if ('MLJIC' %in% criterions) {
      if (unit.lambda) {
        statnames$MLJICunit <- c(statnames$GLM[-11], 'criterion', 'converged')
      }

      statnames$MLJIC <- c('beta01', 'beta1', 'delta', 'se.beta01', 'se.beta1', 'se.delta', 'mad.mu',
                           'rmse.mu', 'R2.mu', 'null.deviance', 'deviance', 'logLik', 'criterion', 'converged')
    }
    nstats <- sapply(statnames, FUN = length)

    ## Function to extract model summaries
    get.stats_k <- function (k, x) {
      x <- x[[k]]

      if (any(class(x) %in% c("simpleError", "error", "condition", "try-error"))) {
        return(rep(NA, nstats[k]))
      }

      if (k == 1) {
        sumx <- catch.conditions({
          summary(x)
        })$value
      }
      else {
        sumx <- catch.conditions({
          summary(x, fisher.matrix = fisher.matrix,
                  sandwich = sandwich,
                  rsquared.method = R2.method,
                  cor.method = cor.method)
        })$value
      }


      if (any(class(sumx) %in% c("simpleError", "error", "condition", "try-error"))) {
        return(rep(NA, nstats[k]))
      }

      if (k == 1) {
        out.stats <- c(sumx$coefficients[,1],
                       sumx$coefficients[,2],
                       mad.mu = mean(abs(x$fitted.values - mu.values), na.rm = TRUE),
                       rmse.mu = sqrt(mean((x$fitted.values - mu.values)^2, na.rm = TRUE)),
                       R2.mu = cor(x$fitted.values, mu.values, method = cor.method)^2,
                       null.deviance = x$null.deviance %||% NA, deviance = x$deviance %||% NA,
                       logLik = logLik (x)[1], converged = x$converged %||% NA)
      }
      else {
        out.stats <- c(sumx$coefficients[,1],
                       sumx$coefficients[,2],
                       mad.mu = mean(abs(sumx$mu - mu.values), na.rm = TRUE),
                       rmse.mu = sqrt(mean((sumx$mu - mu.values)^2, na.rm = TRUE)),
                       R2.mu = cor(sumx$mu, mu.values, method = cor.method)^2,
                       null.deviance = x$null.deviance %||% NA, deviance = x$deviance %||% NA,
                       logLik = logLik (x)[1], criterion = x$criterion[1],
                       converged = x$converged %||% NA)
      }

      if (length(out.stats) != nstats[k]) {
        return(rep(NA, nstats[k]))
      }

      return(out.stats)
    }

    ## Function for one MC run
    one_setting_run <- function(i) {

      if (resample.x) {
        # simulate predictors
        mcdata0 <- as.data.frame (
          sim.predictors (n = n, xdistr = xdistr,
                          xmin = wmin, xmax = wmax)
        )
        colnames(mcdata0) <- 'x'

        eta.lambda <- switch(link,
                             logit = qlogis(lambda),
                             qnorm(lambda))
        mcframe <- msbm.frame(~ x,
                              lambda.formula = ~ 1,
                              data = mcdata0)

        ## Compute true success probabilities
        mu.values <- mutate.params(theta = theta,
                                   frame = mcframe, link = link,
                                   score = FALSE, information = FALSE,
                                   observed = FALSE)$mu
      }


      # simulate response
      mcdata <- mcdata0
      mcdata$y <- simulate(mcframe, nsim = 1, link = link, theta = theta)
      mcdata$y <- as.numeric(mcdata$y)

      ### Get fitted models
      fitmodels <- list()

      # Fit a traditional GLM
      nfits <- 1
      fitmodels$GLM <- try({
        glm (y ~ x, family = binomial(link),
             data = mcdata)
      }, silent = TRUE)

      # Fit MSB models
      # ML fits
      if ('ML' %in% criterions) {
        if (unit.lambda) {
          fitmodels$MLunit <- try({
            msbreg (y ~ x,
                    lambda.formula = ~ 0, link = link,
                    data = mcdata, frame = TRUE,
                    start = beta,
                    control = msbm.control,
                    criterion = 'ML')
          }, silent = TRUE)

          startML <- thetastart
          if (lambda > 0.99 & !identical(class(fitmodels$MLunit), 'try-error')) {
            startML[1:2] <- fitmodels$MLunit$coefficients[1:2]
            startML[3] <- 37
          }

          fitmodels$ML <- try({
            msbreg (y ~ x,
                    lambda.formula = ~ 1, link = link,
                    data = mcdata, frame = TRUE,
                    start = startML,
                    control = msbm.control,
                    criterion = 'ML')
          }, silent = TRUE)

        }
        else {
          fitmodels$ML <- try({
            msbreg (y ~ x,
                    lambda.formula = ~ 1, link = link,
                    data = mcdata, frame = TRUE,
                    start = thetastart,
                    control = msbm.control,
                    criterion = 'ML')
          }, silent = TRUE)
        }
      }

      # MLJ fits
      if ('MLJ' %in% criterions) {
        if (unit.lambda) {
          fitmodels$MLJunit <- try({
            msbreg (y ~ x,
                    lambda.formula = ~ 0, link = link,
                    data = mcdata, frame = TRUE,
                    start = beta,
                    control = msbm.control,
                    criterion = 'MLJ')
          }, silent = TRUE)

          startML <- thetastart
          if (lambda > 0.99 & !identical(class(fitmodels$MLJunit), 'try-error')) {
            startML[1:2] <- fitmodels$MLJunit$coefficients[1:2]
            startML[3] <- 37
          }

          fitmodels$MLJ <- try({
            msbreg (y ~ x,
                    lambda.formula = ~ 1, link = link,
                    data = mcdata, frame = TRUE,
                    start = startML,
                    control = msbm.control,
                    criterion = 'MLJ')
          }, silent = TRUE)

        }
        else {
          fitmodels$MLJ <- try({
            msbreg (y ~ x,
                    lambda.formula = ~ 1, link = link,
                    data = mcdata, frame = TRUE,
                    start = thetastart,
                    control = msbm.control,
                    criterion = 'MLJ')
          }, silent = TRUE)
        }

      }

      # MLPJ fits
      if ('MLPJ' %in% criterions) {
        if (unit.lambda) {
          fitmodels$MLPJunit <- try({
            msbreg (y ~ x,
                    lambda.formula = ~ 0, link = link,
                    data = mcdata, frame = TRUE,
                    start = beta,
                    control = msbm.control,
                    criterion = 'MLPJ')
          }, silent = TRUE)

          startML <- thetastart
          if (lambda > 0.99 & !identical(class(fitmodels$MLPJunit), 'try-error')) {
            startML[1:2] <- fitmodels$MLPJunit$coefficients[1:2]
            startML[3] <- 37
          }

          fitmodels$MLPJ <- try({
            msbreg (y ~ x,
                    lambda.formula = ~ 1, link = link,
                    data = mcdata, frame = TRUE,
                    start = startML,
                    control = msbm.control,
                    criterion = 'MLPJ')
          }, silent = TRUE)

        }
        else {
          fitmodels$MLPJ <- try({
            msbreg (y ~ x,
                    lambda.formula = ~ 1, link = link,
                    data = mcdata, frame = TRUE,
                    start = thetastart,
                    control = msbm.control,
                    criterion = 'MLPJ')
          }, silent = TRUE)
        }

      }

      # MLJIC fits
      if ('MLJIC' %in% criterions) {
        if (unit.lambda) {
          fitmodels$MLJICunit <- try({
            msbreg (y ~ x,
                    lambda.formula = ~ 0, link = link,
                    data = mcdata, frame = TRUE,
                    start = beta,
                    control = msbm.control,
                    criterion = 'MLJIC')
          }, silent = TRUE)
        }

        fitmodels$MLJIC <- try({
          msbreg (y ~ x,
                  lambda.formula = ~ 1, link = link,
                  data = mcdata, frame = TRUE,
                  start = thetastart,
                  control = msbm.control,
                  criterion = 'MLJIC')
        }, silent = TRUE)
      }

      ## Get model fit summaries
      get.stats <- function (x) {
        outk <- lapply(1:length(x),
                       FUN = get.stats_k, x = x)
        outk <- do.call("c", outk)

        c(outk, mean.y = mean(mcdata$y, na.rm = TRUE))
      }

      out <- get.stats (fitmodels)

#      if (any(is.na(out))) {
#        browser()
#      }

      out
    }

    # Run B MC simulations
    outmat <- catch.conditions({
      lapply(1:B, FUN = one_setting_run)
    })$value

    if (any(class(outmat) %in% c("simpleError", "error",
                                 "condition", "try-error"))) {
      return(NA)
    }

    nsimcol <- sapply(outmat, FUN = length)
    if (!all(nsimcol == max(nsimcol))) {
      maxsim <- which.max(nsimcol)[1]
      outmat <- c(outmat[[maxsim]], outmat[-maxsim])
    }
    outmat <- do.call("rbind", outmat)

    # Vector of statistic names
    stat.names <- lapply(1:length(statnames), FUN = function(k) {
      paste0('M', k, '.', statnames[[k]])
    })
    stat.names <- c(do.call("c", stat.names), 'mean.y')

    if (NCOL(outmat) == length(stat.names))
      colnames(outmat) <- stat.names

    # Save results
    Bresults <- list(xdistr = rep(xdistr, NROW(outmat)),
                     n = n, lambda = lambda)
    Bresults <- c(Bresults, as.list(beta))
    Bresults <- do.call("cbind", Bresults)
    Bresults <- as.data.frame(Bresults)

    Bresults <- cbind(Bresults, outmat)
    Bresults$n <- as.numeric(Bresults$n)
    Bresults$lambda <- as.numeric(Bresults$lambda)
    Bresults$beta01 <- as.numeric(Bresults$beta01)
    Bresults$beta1 <- as.numeric(Bresults$beta1)

    cat("\n estimates \n ")
    print(str(outmat))
    print(colnames(outmat))

    return(Bresults)
  }

  Braw <- matteLapply (X = 1:Nsettings,
                       FUN = onesetting,
                       cl = cl, chunk.size = chunk.size)
  if (!is.null(cl)) {
    parallel::stopCluster(cl = cl)
  }

  nsimcol <- sapply(Braw, FUN = length)
  if (!all(nsimcol == max(nsimcol))) {
    maxsim <- which.max(nsimcol)[1]
    Braw <- c(Braw[[maxsim]], Braw[-maxsim])
  }

  Braw <- as.data.frame(do.call('rbind', Braw))

  if (!is.null(savepath)) {
    savepath <- paste0(savepath, 'primerun_')
    write.table(Braw, file  = paste0(savepath, 'Braw.txt'))
  }

  return(Braw)

}


