#'
#' Monte Carlo Simulations
#'
#' Simulate data, fit Multistage Binomial (MSB) models and
#' extract some estimates/statistics: repeat this \code{nsim}
#' times and return the replicates of estimates/statistics.
#' Or return a function to replicate the estimates/statistics.
#'
#' @param msbm a description of the true MSB model (from which
#' data will be generated). This can be any of:
#'
#' \itemize{
#' \item a *MSB model frame* of class \code{"msbm.frame"};
#' \item a *MSB model fit* of class \code{"msbm"};
#' \item a *MSB model summary* of class \code{"summary.msbm"}; or,
#' \item a *character* string giving a \link{call} to \link[msbreg]{msbreg}
#' to obtain a model fit (of class \code{"msbm"}).
#' }
#'
#' @param theta numeric vector, model parameters to generate data.
#' The argument \code{theta} **should always** be supplied.
#' For \code{MCmsbm}, when the argument \code{msbm}
#' is of class \code{"msbm"} (or a call which evaluates to that
#' class of object), the default is the \code{$coefficients}
#' component of \code{msbm}.
#'
#' @param models either a list of \code{"msbm"} (or \code{"msbm.frame"})
#' class objects, or a character vector (each element describing a \link{call}
#' to \link[msbreg]{msbreg}). This specifies model(s) to be fitted to each
#' simulated dataset. If \code{NULL} (the default), the data model \code{msbm}
#' is used.
#'
#' @param link a \code{"link"} object, or a function that produces
#' a \code{"link"} class object (or a character giving the name of
#' such a function). See further description in \link[msbreg]{msbreg}.
#'
#' @param basic.stats character vector indicating basic (pre-defined)
#' quantities to include in the set of returned statistics, in addition
#' to the vector component \code{$coefficients} and related standard errors
#' from \code{"msbm"} fits.
#' Currently, \code{basic.stats} can include any or many of:
#' \describe{
#' \item{\code{"mean.y"}}{ sample average of the binary response
#' (\code{y/weights} where \code{weights} denotes the vector of
#' number of trials);}
#' \item{\code{"mad.mu"}}{ mean absolute deviation of the fitted
#' success probabilities from the true success probabilities;}
#' \item{\code{"rmse.mu"}}{ root mean square error of the fitted
#' success probabilities from the true success probabilities;}
#' \item{\code{"R2.mu"}}{ squared correlation between the fitted
#' success probabilities from the true success probabilities;}
#' \item{\code{"null.deviance"}, \code{"deviance"}, \code{"aic"},
#' \code{"bic"}, \code{"hqc"}}{ the
#' corresponding components from a \code{"msbm"} fit;}
#' \item{\code{"logLik"}}{ the log-likelihood value returned by
#' \link[stats]{logLik};}
#' \item{\code{"criterion"}}{ the optimal value of the fitting criterion
#' (the \code{"criterion"} component of a \code{"msbm"} fit);}
#' \item{\code{"r.squared"}}{ the pseudo-\eqn{R^2} value(s) returned by
#' \link[msbreg]{rsquared};}
#' \item{\code{"converged"}}{ logical indicating convergence;}
#' \item{\code{"emp.lambda"}}{ an empirical estimate of the model component
#' \eqn{\lambda} (maximum success probability of the binary response) computed
#' using \link[msbreg]{infer.lambda};}
#' \item{\code{adjust.lambda}}{ an adjusted empirical estimate of \eqn{\lambda},
#' the raw empirical estimate is adjusted by setting regression coefficients
#' to their estimates;}
#' \item{\code{"all"}}{ all of the above choices.}
#' }
#'
#' Defaults to \code{"all"}.
#'
#' @param sup.stats a function, or a character naming a function to compute
#' some interest statistics from an \code{"msbm"} class object.
#' If \code{sup.args} is \code{NULL}, the function \code{sup.stats} must
#' take one argument and return a numeric vector. Otherwise, \code{sup.stats}
#' must take two arguments and return a numeric vector.
#'
#' @param sup.args (list of) supplemental argument(s) for the function
#' \code{sup.stats}. If supplied (and not \code{NULL}), \code{sup.args}
#' is passed as a second argument to \code{sup.stats} (an object of class
#' \code{"msbm"} being the first argument).
#'
#' @param fisher.matrix logical, should the Fisher (expected) information
#' matrix be used to obtain standard errors? Defaults to \code{TRUE}.
#' The alternative is the observed information matrix.
#'
#' @param adjust,sandwich additional arguments to adjust the default
#' standard error estimators, passed to \link[msbreg]{summary.msbm}.
#'
#' @param R2.method character, indicates the pseudo-\eqn{R^2} statistic
#' to include into the basic statistics. Only used when
#' \code{basic.stats} includes \code{"r.squared"} (or \code{"all"}).
#' See argument \code{method} of \link[msbreg]{rsquared} for possible
#' values of \code{R2.method}.
#'
#' @param cor.method character indicating which correlation coefficient
#' to consider (only used when \code{basic.stats} includes \code{"r.squared"}
#' or \code{"all"}, and \code{R2.method = 'COR'}).
#'
#' @param method character, the method to be used to run the simulations.
#' The default \code{method = "sim.boot"} uses the function \link[boot]{boot}
#' from the package \code{boot} to run the Monte Carlo simulations as
#' a parametric bootstrap. See **Details**.
#'
#' @param nsim numeric, number of simulations (replicates) to run.
#' Defaults to \code{nsim = 999}.
#'
#' @param seed a single value (random seed), interpreted as an integer,
#' or \code{NULL}. It is passed to \link{set.seed} when not \code{NULL}.
#'
#' @param parallel,ncpus,cl named arguments to allow parallel computations.
#' See **Details**.
#'
#' @param ... further named arguments passed to or from other methods.
#' Currently, none is used.
#'
#' @param statistics a function (or a character string naming a function)
#' which takes a dataset with the same shape as the one in argument
#' \code{msbm} and returns a vector of statistic(s) of interest.
#'
#' @param data the dataset used to create \code{msbm} (the \code{data}
#' argument in \code{msbm} when the latter is a character).
#'
#' @param stat.names character, optional vector of names of all statistics
#' returned by the function \code{statistics}. Only used when the length
#' of \code{stat.names} matches the length of the output from \code{statistics}.
#'
# @usage
# MCmsbm (msbm, theta = msbm$coefficients,
#         models = NULL, link = msbm$link, basic.stats = "all",
#         sup.stats = NULL, sup.args = NULL, fisher.matrix = TRUE,
#         R2.method = c("KL", "Nagelkerke", "COR"), cor.method = "pearson",
#         method = "sim.boot", nsim = 999, seed = NULL,
#         parallel = c("no", "multicore", "snow"),
#         ncpus = getOption("boot.ncpus", 1L), cl = NULL, ...)
#
# sim.boot (statistics, theta, data, nsim,
#           parallel, ncpus, cl, stat.names = NULL, ...)
#'
#' @details
#' The function \code{MCmsbm} provides a frame to run simulations
#' on MSB models with minimum coding requirements.
#' Depending on the argument \code{method}, the function performs
#' one of the two actions:
#' \describe{
#' \item{\code{method = "sim.frame"}:}{ builds a function named \code{statistics}
#' that generates a dataset (using \code{msbm} and \code{theta}),
#' fits all fitting models (given in \code{models}), and returns a vector of
#' statistics including the \code{$coefficients} component from each fitting
#' model, the related estimates of standard errors, any additional statistics
#' specified through the \code{character} argument \code{basic.stats}, and any
#' further statistics specified through the \code{function} argument
#' \code{sup.stats};}
#' \item{\code{method = "sim.boot"}}{ (the current default) performs the action
#' corresponding to \code{method = "sim.frame"} and then calls the function
#' \link[boot]{boot} from the \code{boot} package to run \code{nsim}
#' Monte Carlo simulations (that is, a *parametric bootstrap*) based on
#' the built function \code{statistics}.}
#' }
#'
#' The argument \code{method} can however be replaced by (the name of) a
#' user-supplied function that takes the same arguments as \code{sim.boot}
#' and returns a list including at least a component \code{$estimates} (a
#' matrix with \code{nsim} rows and \eqn{p} columns with \eqn{p} the number
#' of element returned by the function \code{statistics}).
#'
#' The option \code{method = "sim.frame"} can be used to prepare a
#' generating function for Monte Carlo simulations using a simple call
#' to the \link{apply}-like family of functions. Its use also extends to
#' running Monte Carlo simulations using the function
#' \link[MonteCarlo]{MonteCarlo} from package \code{MonteCarlo}
#' (see **Examples**).
#'
#' The argument \code{seed} is used to set the \code{.Random.seed} of
#' \code{R}, i.e. the random number generator (RNG) state for random number
#' generation. Although the simulation uses random numbers, a call to the
#' function will not alter the RGN state if \code{.Random.seed}
#' exists in the \code{.GlobalEnv}. If \code{.GlobalEnv} did not exist,
#' it is first created by calling \code{runif(1)}. In either case, the
#' RGN state (i.e. \code{.Random.seed}) right before bootstrapping is
#' performed is saved and restored after bootstrapping.
#'
#' On multicore plateformes, when \code{method = "sim.boot"}, consider
#' using arguments \code{parallel}, \code{ncpus}, and \code{cl} to allow
#' parallel computations which can be faster. This is currently based on
#' the \code{R} package \code{parallel} or \code{snow} (see \link[boot]{boot}).
#' Arguments \code{parallel}, \code{ncpus}, and \code{cl} are also
#' passed to any user-supplied \code{method} which may allow parallel
#' computations.
#'
#' @return When \code{method = "sim.frame"} an object inheriting from class
#' \code{"msbm.sim.frame"} which is a list with the following elements:
#'
# \describe{
#' \item{\code{statistics}}{ a function that generates data, fits all models,
#' and returns estimates/statistics from all model fits;}
#' \item{\code{stat.names}}{ character vector giving the names of the
#' estimates/statistics that the function \code{statistics} returns;}
#' \item{\code{nstats}}{ numeric vector giving the number of estimates/statistics
#' from each model fit;}
#' \item{\code{npars}}{ numeric vector giving the number of estimated
#' parameters for each model fit;}
#' \item{\code{frame}}{ model frame (\code{"msbm.frame"} class object) of
#' of the true data model;}
#' \item{\code{link}}{ the used \code{link} object;}
#' \item{\code{theta}}{ numeric vector of true model parameters;}
#' \item{\code{models}}{ list of \code{call}s of model fits.}
# }
#'
#' When \code{method = "sim.boot"} or supplied by the user, the returned
#' object inherits from class \code{"msbm.sim"} which has the elements
#' returned for \code{method = "sim.frame"} and the additional elements:
#'
#' \item{\code{estimates}}{ numeric matrix of estimates/statistics with
#' \code{nsim} rows and \code{sum(nstats)} columns;}
#' \item{\code{RNGstate}}{ the state of the random number generator
#' right before the simulation;}
#' \item{\code{raw}}{ list of additional information.
#' For \code{method = "sim.boot"}, an object of class \code{"boot"} without
#' its \code{$t} component (the latter is the matrix \code{estimates}).}
#'
#' An object of class \code{"msbm.sim"} has a \code{print} method for a quick
#' display of summary statistics for the columns of the \code{$estimates}
#' component.
#'
#' @export MCmsbm
#' @aliases sim.boot
#' @export sim.boot
#'
#' @seealso \link[msbreg]{bootfit.msbm} for bootstrapping a multistage
#' binomial model fit.
#'
#' @examples
#' ## Generate some data for the simulations
#' set.seed(167)
#' mcdf <- data.frame (x1 = rexp(1000, 0.5),
#'                     x2 = rnorm(1000))
#' mcdf$y <- sim.msb (formula = ~ x1 | x2,
#'                    data = mcdf,
#'                    link = "probit",
#'                    theta = c(3, -1, 3, -1, 1.6449),
#'                    seed = NULL)$y
#'
#' head(mcdf)
#'
#' ## Fitting the true model
#' msbfit <- msbreg(formula = y ~ x1 | x2,
#'                  link = 'probit', data = mcdf,
#'                  criterion = 'MLJ', frame = TRUE)
#'
#' summary (msbfit)
#' \dontrun{
#' ## Run nsim = 10 MC simulations: this takes about 6 seconds
#' MCMSBfit = MCmsbm (msbm = msbfit,
#'                    theta = c(3, -1, 3, -1, 1.6449),
#'                    nsim = 10, seed = 1)
#'
#' ## Summary of returned statistics
#' MCMSBfit
#'
#' Table = data.frame(Truth = c(3, -1, 3, -1, 1.6449),
#'                    Mean = colMeans(MCMSBfit$estimates[,1:5]),
#'                    SD = colSds(MCMSBfit$estimates[,1:5]),
#'                    Mean.se = colMeans(MCMSBfit$estimates[,6:10]))
#' Table$`bias (%)` <- 100 * (Table$Mean - Table$Truth) / abs(Table$Truth)
#' Table$`rrmse (%)` <- 100 * (Table$SD^2 + (Table$Mean - Table$Truth)^2) /
#'                       abs(Table$Truth)
#' Table$`bias.se (%)` <- 100 * (Table$Mean.se - Table$SD) / abs(Table$SD)
#'
#' print(Table, digits = 4)
#'
#' ### Note: we need a large number of replicates (say nsim = 1000),
#' # to approach the true behavior of the estimates. With nsim = 1000,
#' # the simulations take about 20 minutes on MacOS (16GB) M1 Pro and we get:
#' #
#' #                     Truth    Mean      SD Mean.se bias (%) rrmse (%) bias.se (%)
#' # (Intercept).1       3.000  3.0569 0.31077  0.3223    1.896     3.327       3.713
#' # x1                 -1.000 -1.0171 0.09901  0.1022   -1.712     1.010       3.207
#' # (Intercept).2       3.000  2.9282 0.89664  0.8390   -2.393    26.970      -6.427
#' # x2                 -1.000 -0.9887 0.42621  0.4041    1.125    18.178      -5.190
#' # (Intercept).lambda  1.645  1.7248 0.17298  0.2078    4.858     2.207      20.155
#' ###
#' ### See `?tsbmmcrun` for all results
#'
#' ## Run nsim = 10 MC simulations with the true model
#' #  but also fit a model with lambda = 1
#' # Fit the model with lambda = 1
#' msbfit0 <- msbreg(formula = y ~ x1 | x2,
#'                   lambda.formula = ~ 0,
#'                   link = 'probit', data = mcdf,
#'                   criterion = 'MLJ', frame = TRUE)
#'
#' # Run simulations
#' MCMSBfit2 = MCmsbm (msbm = msbfit,
#'                     theta = c(3, -1, 3, -1, 1.6449),
#'                     models = list(msbfit0, msbfit),
#'                     nsim = 10, seed = 1)
#'
#' # Summary of returned statistics
#' MCMSBfit2
#' }
#'
MCmsbm <- function (msbm,
                    theta = msbm$coefficients,
                    models = NULL, # Fitting model(s) to be fitted to simulated data
                    link = msbm$link,
                    basic.stats = "all", # Statistics to return
                    sup.stats = NULL, sup.args = NULL, # Statistics to return
                    fisher.matrix = TRUE,
                    adjust = FALSE, sandwich = FALSE,
                    R2.method = c("KL", "Nagelkerke", "COR"),
                    cor.method = "pearson",
                    method = "sim.boot",
                    nsim = 999, seed = NULL,
                    parallel = c("no", "multicore", "snow"),
                    ncpus = getOption("boot.ncpus", 1L), cl = NULL, ...) {
  #* Matched call
  mcall <- match.call()

  #* Save/Set the random generator seed
  eval(setsave.RNGstate())

  #* True model (model frame, TRUE parameter vector) and data
  envCall <- parent.frame()
  eval(toget_TRUEModelParData())

  #* Link function
  if (!missing(link)) {
    eval(get.linkfun())
  }
  else if (!inherits(msbm, "msbm")) {
    link <- msbreg::logit()
  }

  #* Basic statistics to be extracted from each fit
  eval(check.basic.sim.stats())

  #* Fitting models
  eval(toget_FittingModels())
  nobjects <- length(objects)
  for (modk in seq_len(nobjects)) {
    if (any(class(objects[[modk]]) %in% c("simpleError", "error",
                                       "condition", "try-error"))) {
      stop(paste0("initial model evaluation failled: check element number ", modk,
                  " of argument 'models'"))
    }
  }
  npars <- sapply(objects, FUN = function(objectj) objectj$dims$npars)

  #* Matched arguments
  R2.method <- match.arg(R2.method)
  cor.method <- match.arg(cor.method)

  #* Sup. stats
  eval(toget_SUP.sim.stats())

  #* Compute true success probabilities
  mu.values <- mutate.params(theta = theta,
                             frame = mf, link = link,
                             score = FALSE, information = FALSE,
                             observed = FALSE)$mu

  #* Function to extract statistics from a fit
  eval(toget_get_get.stats_0())
  statsALL <- lapply (1:nobjects,
                      FUN = get.stats_0,
                      x = objects)
  NAstats <- sapply(statsALL, FUN = function(x) {
    length(x) == 1 & all(is.na(x))
  })
  if (any(NAstats)) {
    stop(paste0("initial model evaluation failled: check element number ", which(NAstats),
                " of argument 'models'"))
  }
  nstats <- sapply(statsALL, length)
  eval(toget_get.sim.stats())

  #* Function to generate data replicates based on the model frame and theta;
  #* and fit the "fitting" models
  fitkfun <- function(callk, data) {
    callk$data <- data
    callk <- as.call(callk)
    outk <- catch.conditions({
      eval(callk, envir = envCall)
    })$value
  }
  calculate.sim.statistics <- function (...) {
    datak <- data
    datak[[mf$auxy$ycolname]] <- as.numeric(simulate(mf, theta = theta, link = link))

    newfits <- lapply(model.calls, FUN = fitkfun, data = datak)

    gotstats <- get.stats (newfits)

    gotstats
  }

  # Names of estimates/statistics
  eval(toget_ParsNames.stats())

  #* Call the workhorse function
  if (identical(method, "sim.frame")) {
    out <- structure(list(statistics = calculate.sim.statistics,
                          stat.names = outnames,
                          nstats = nstats,
                          frame = mf,
                          link = link,
                          theta = theta,
                          models = model.calls,
                          call = mcall),
                     class = "msbm.sim.frame")
  }
  else {

    out <- eval(call(if (is.function(method)) "method" else method,
                     statistics = calculate.sim.statistics,
                     theta = theta, frame = mf, link = link,
                     data = data, nsim = nsim,
                     parallel = parallel, ncpus = ncpus, cl = cl,
                     stat.names = outnames))

    out <- structure(c(out,
                       list(statistics = calculate.sim.statistics,
                            stat.names = outnames,
                            nstats = nstats,
                            npars = npars,
                            frame = mf,
                            link = link,
                            theta = theta,
                            models = model.calls,
                            call = mcall,
                            RNGstate = RNGstate)),
                     class = "msbm.sim")
  }

  out
}

toget_TRUEModelParData <- function() {
  expression({
    #* True model (model frame, TRUE parameter vector) and data
    if (is.character(msbm)) {
      # Turn a character argument 'msbm' into a 'msbm' class object (fit)
      msbm <- str2lang (msbm)
      if (is.null(msbm$link))
        msbm$link <- if(!missing(link)) link else "logit"
      msbm[[1L]] <- quote(msbreg::msbreg)
      msbme <- catch.conditions({
        eval(msbm, envir = envCall)
      })$value
      if (any(class(msbme) %in% c("simpleError", "error",
                               "condition", "try-error"))) {
        for (jpframe in 2:100) {
          msbme <- catch.conditions({
            eval(msbm, parent.frame(jpframe))
          })$value
          if (!any(class(msbme) %in% c("simpleError", "error",
                                      "condition", "try-error"))) {
            envCall <- parent.frame(jpframe)
            break
          }
        }

        if (any(class(msbme) %in% c("simpleError", "error",
                                     "condition", "try-error"))) {
          cat("\n jpframe: ", jpframe, "\n")
          stop(msbme)
        }
      }
      msbm <- msbme; rm(msbme)
    }

    if (inherits(msbm, "msbm") | inherits(msbm, "summary.msbm")) {
      if (inherits(msbm, "summary.msbm")) {
        msbm <- eval(msbm$call, envir = envCall)
      }
      mf <- model.frame(msbm)
      if (!missing(theta)) {
        stopifnot(is.numeric(theta))
        stopifnot(length(theta) == length(msbm$coefficients))
      }
      data <- msbm$data
    }
    else if(inherits(msbm, "msbm.frame")) {
      mf <- msbm
      if(missing(theta))
        stop("argument 'theta' must be supplied if 'msbm' is of class `msbm.frame`")
      stopifnot(is.numeric(theta))
      stopifnot(length(theta) == length(msbm$parnames))

      if (!is.null(msbm$call$data)) {
        data <- eval(msbm$call$data, envir = envCall)
      }
      else {
        data <- environment(msbm$call$formula)
        stopifnot(!is.null(data))
      }
    }
  })
}

check.basic.sim.stats <- function() {
  expression({
    #* Basic statistics to be extracted from each fit
    available.stats <- c("all", "emp.lambda", "adjust.lambda", "mean.y", "mad.mu", "rmse.mu",
                         "R2.mu", "null.deviance", "deviance", "logLik", "criterion",
                         "aic", "bic", "hqc", "r.squared", "converged")
    if (!is.null(basic.stats)) {
      stopifnot(is.character(basic.stats))
      stopifnot(all(basic.stats %in% available.stats))
      if (any(basic.stats %in% "all"))
        basic.stats <- available.stats[-1]
    }
    else
      basic.stats <- character(0)
    available.stats <- available.stats[-1]
    add.stats <- available.stats %in% basic.stats
  })
}

toget_FittingModels <- function() {
  expression({
    #* Fitting models
    if (is.null(models)) {
      nbModels <- 1
      model.calls <- msbm$call
      if(inherits(msbm, "msbm.frame")) {
        model.calls[[1L]] <- quote(msbreg::msbreg)
      }
      model.calls <- list(model.calls)

      # Evaluate models in 'models': try model fitting using each call in 'model.calls'
      if (inherits(msbm, "msbm")) {
        objects <- list(msbm)
      }
      else {
        objects <- catch.conditions({
          eval(model.calls[[1]], envir = envCall)
        })$value

        if (any(class(objects) %in% c("simpleError", "error",
                                    "condition", "try-error"))) {
          stop("initial model evaluation failled: check argument 'models'")
        }

        objects <- list(objects)

      }
    }
    else {
      objects <- NULL
      nbModels <- length(models)
      if (is.character(models)) {
        # Get model calls from characters
        model.calls <- vector(mode = "list", length = nbModels)
        for (k in 1:nbModels) {
          model.calls[[k]] <- str2lang (models[k])
          if (is.null(model.calls[[k]]$link))
            model.calls[[k]]$link <- if(!missing(link)) link else "logit"
          model.calls[[k]][[1L]] <- quote(msbreg::msbreg)
        }

        # Get model fits from model calls
        models <- vector(mode = "list", length = nbModels)
        for (k in 1:nbModels) {
          models[[k]] <- catch.conditions({
            eval(model.calls[[k]], envir = envCall)
          })$value

          if (any(class(models[[k]]) %in% c("simpleError", "error",
                                        "condition", "try-error"))) {
            stop(paste0("initial model evaluation failled: check element number ", k,
                        " of argument 'models'"))
          }
        }
        objects <- models
      }

      if(!is.list(models)) {
        stop("'models' must be a list or a character vector")
      }

      # Extract model calls from fits
      model.calls <- lapply(models, FUN = function(x) {
        x <- x$call

        if(inherits(x, "msbm.frame")) {
          x[[1L]] <- quote(msbreg::msbreg)
        }

        return(x)
      })

      # Evaluate models in 'models': try model fitting using each call in 'model.calls'
      if (is.null(objects)) {
        objects <- vector(mode = "list", length = nbModels)
        for (k in 1:nbModels) {
          if (inherits(models[[k]], "msbm")) {
            objects[[k]] <- models[[k]]
          }
          else {
            objects[[k]] <- catch.conditions({
              eval(model.calls[[k]], envir = envCall)
            })$value

            if (any(class(objects[[k]]) %in% c("simpleError", "error",
                                              "condition", "try-error"))) {
              stop(paste0("initial model evaluation failled: check element number ", k,
                          " of argument 'models'"))
            }
          }
        }
      }
    }
  })
}

toget_SUP.sim.stats <- function() {
  expression({
    if (is.null(sup.stats)) {
      NULLsup.stats <- TRUE
      sup.stats <- function(x, ...) NULL
    }
    else {
      NULLsup.stats <- FALSE
      if (is.character(sup.stats))
        sup.stats <- get(sup.stats, mode = "function", envir = parent.frame())
      if (!is.function(sup.stats))
        stop("the supplied 'sup.stats' is not a (name of a) function")

      if (!is.null(sup.args)) {
        sup.stats0 <- sup.stats
        sup.stats <- function(x, ...) {
          sup.stats0 (x, sup.args)
        }
      }
    }
  })
}

toget_get_get.stats_0 <- function () {
  expression({
    get.stats_0 <- function (k, x) {
      x <- x[[k]]

      if (any(class(x) %in% c("simpleError", "error", "condition", "try-error"))) {

        return(NA)
      }

      sumx <- catch.conditions({
        summary(x, fisher.matrix = fisher.matrix,
                adjust = adjust, sandwich = sandwich,
                rsquared.method = R2.method,
                cor.method = cor.method)
      })$value

      if (any(class(sumx) %in% c("simpleError", "error", "condition", "try-error"))) {
        return(NA)
      }

      if ("logLik" %in% available.stats[add.stats])
        sumx$logLik <- logLik (x)[1]

      if ("criterion" %in% available.stats[add.stats])
        sumx$criterion <- x$criterion[1]

      if ("mad.mu" %in% available.stats[add.stats])
        sumx$mad.mu <- mean(abs(sumx$mu - mu.values), na.rm = TRUE)

      if ("rmse.mu" %in% available.stats[add.stats])
        sumx$rmse.mu <- sqrt(mean((sumx$mu - mu.values)^2, na.rm = TRUE))

      if ("R2.mu" %in% available.stats[add.stats])
        sumx$R2.mu <- cor(sumx$mu, mu.values, method = cor.method)^2

      if ("mean.y" %in% available.stats[add.stats])
        sumx$mean.y <- mean(sumx$mu + x$y.resid/x$weights, na.rm = TRUE)

      if (any(c("emp.lambda", "adjust.lambda") %in% available.stats[add.stats]))
        emp.L <- empirical.lambda (x, ...)

      if ("emp.lambda" %in% available.stats[add.stats])
        sumx$emp.lambda <- attr(emp.L, "mu_ab")

      if ("adjust.lambda" %in% available.stats[add.stats])
        sumx$adjust.lambda <- emp.L[1]

      if (!NULLsup.stats) {
        sup.out <- if (is.null(sup.args)) sup.stats (x) else  sup.stats (x, sup.args)
      }
      else
        sup.out <- NULL

      out.stats <- c(sumx$coefficients[,1],
                     sumx$coefficients[,2],
                     if (any(add.stats)) unlist(lapply(available.stats[add.stats],
                                                       FUN = function (stat) {
                                                         sumx[stat]
                                                       })),
                     sup.out)

      return(out.stats)
    }
  })
}

toget_get.sim.stats <- function () {
  expression({
    get.stats_k <- function (k, x) {
      x <- x[[k]]

      if (any(class(x) %in% c("simpleError", "error", "condition", "try-error"))) {

        return(rep(NA, nstats[k]))
      }

      sumx <- catch.conditions({
        summary(x, fisher.matrix = fisher.matrix,
                adjust = adjust, sandwich = sandwich,
                rsquared.method = R2.method,
                cor.method = cor.method)
      })$value

      if (any(class(sumx) %in% c("simpleError", "error", "condition", "try-error"))) {
        return(rep(NA, nstats[k]))
      }

      if ("logLik" %in% available.stats[add.stats])
        sumx$logLik <- logLik (x)[1]

      if ("criterion" %in% available.stats[add.stats])
        sumx$criterion <- x$criterion[1]

      if ("mad.mu" %in% available.stats[add.stats])
        sumx$mad.mu <- mean(abs(sumx$mu - mu.values), na.rm = TRUE)

      if ("rmse.mu" %in% available.stats[add.stats])
        sumx$rmse.mu <- sqrt(mean((sumx$mu - mu.values)^2, na.rm = TRUE))

      if ("R2.mu" %in% available.stats[add.stats])
        sumx$R2.mu <- cor(sumx$mu, mu.values, method = cor.method)^2

      if ("mean.y" %in% available.stats[add.stats])
        sumx$mean.y <- mean(sumx$mu + x$y.resid/x$weights, na.rm = TRUE)

      if (any(c("emp.lambda", "adjust.lambda") %in% available.stats[add.stats]))
        emp.L <- empirical.lambda (x, ...)

      if ("emp.lambda" %in% available.stats[add.stats])
        sumx$emp.lambda <- attr(emp.L, "mu_ab")

      if ("adjust.lambda" %in% available.stats[add.stats])
        sumx$adjust.lambda <- emp.L[1]

      if (!NULLsup.stats) {
        sup.out <- if (is.null(sup.args)) sup.stats (x) else  sup.stats (x, sup.args)
      }
      else
        sup.out <- NULL

      out.stats <- c(sumx$coefficients[,1],
                     sumx$coefficients[,2],
                     if (any(add.stats)) unlist(lapply(available.stats[add.stats],
                                                       FUN = function (stat) {
                                                         sumx[stat]
                                                       })),
                     sup.out)

      if (length(out.stats) != nstats[k]) {
        return(rep(NA, nstats[k]))
      }

      return(out.stats)
    }

    get.stats <- function (x) {
      outk <- lapply(1:length(x),
                     FUN = get.stats_k, x = x)
      do.call("c", outk)
    }
  })
}

toget_ParsNames.stats <- function() {
  expression({
    fpnamesfun <- function (k) {
      outnamesk <- names(statsALL[[k]])
      thetak0 <- objects[[k]]$coefficients
      nparsk <- length(objects[[k]]$coefficients)

      if (length(outnamesk) == length((statsALL[[k]]))) {

        outnamesk[(nparsk + 1):(2 * nparsk)] <-
          paste0("se.", outnamesk[(nparsk + 1):(2 * nparsk)])
      }
      else {

      outnamesk <- names(objects[[k]]$coefficients)
      if (length(outnamesk) != nparsk) {
        outnamesk <- model.frame(objects[[k]])$parnames

        if (length(outnamesk) != nparsk) {
          outnamesk <- paste0("coef", if (nobjects > 1) k,
                              ".", 1:nparsk)
        }
      }
      outnamesk <- c(outnamesk, paste0("se.", outnamesk))
      if (any(add.stats)) {
        outnamesk <- c(outnamesk, available.stats[add.stats])
        if ("r.squared" %in% available.stats) {
          outnamesk <- c(outnamesk[1:(length(outnamesk) - 1)],
                         paste0(R2.method, ".R2"))
        }
      }

      if (!NULLsup.stats) {
        valsup.statk <- if (is.null(sup.args)) sup.stats (objects[[k]])
        else sup.stats (objects[[k]], sup.args)
        nsup.statk <- length(valsup.statk)
        if (nsup.statk) {
          if (length(names(valsup.statk)) == nsup.statk)
            outnamesk < c(outnamesk, names(valsup.statk))
          else
            outnamesk < c(outnamesk,
                          paste0("sup", if (nobjects > 1) k,
                                 ".", 1:nsup.statk))
        }
      }
      }

      if (length(statsALL) > 1)
        outnamesk <- paste0("M", k, ".", outnamesk)

      return(make.unique(outnamesk))
    }

    outnames <- lapply(1:length(statsALL), FUN = fpnamesfun)
    outnames <- do.call("c", outnames)
    outnames <- make.unique(outnames)
  })
}
