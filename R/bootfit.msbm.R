#' Bootstrap a Multistage Binomial Model Fit
#'
#' This is a method of \link[msbreg]{bootfit} for class \code{"msbm"} objects.
#' It bootstraps a Multistage Binomial (MSB) model fit and returns an object
#' inheriting from class "\link[boot]{boot}".
#'
#' @param object an object of class \code{"msbm"}, typically returned
#' by \link[msbreg]{msbreg}.
#'
#' @param coefficients numeric vector, model parameters to generate data
#' for parametric bootstrap (when \code{sim = 'parametric'}).
#' The default is the \code{$coefficients} component of \code{object}.
#'
#' @param R,sim,stype,simple,parallel,ncpus,cl named arguments passed to \link[boot]{boot}.
#'
#' @param seed a single value (random seed), interpreted as an integer,
#' or \code{NULL}. It is passed to \link{set.seed} when not \code{NULL}.
#' Same use as for \link[msbreg]{bootfit} (see **Details** therein).
#'
#' @param basic.stats character vector indicating basic (pre-defined)
#' statistics to include in the set of bootstrapped statistics, in addition
#' to the vector component \code{$coefficients} of \code{"msbm"} fits.
#' The default (\code{NULL}) means
#' none. Currently, \code{basic.stats} can include any or many of:
#'
#' \describe{
#' \item{"\code{r.squared}"}{ to include some pseudo-\eqn{R^2} measures;}
#' \item{"\code{alpha}"}{ to include the \code{$alpha} component of
#' \code{summary(object)};}
#' \item{"\code{lambda}"}{ to include the \code{$lambda} component of
#' \code{summary(object)};}
#' \item{"\code{alpha.lambda}"}{ to include the \code{$alpha.lambda} component
#' of \code{summary(object)};}
#' \item{"\code{mu}"}{ to include the \code{$mu} component
#' of \code{summary(object)};}
#' \item{"\code{se.fit}"}{ to include estimated standard errors for estimated
#' model parameters and quantities such as \code{alpha}, \code{lambda},
#' \code{alpha.lambda}, and \code{mu} (standard errors are useful for
#' instance for building a Studentized bootstrap confidence interval,
#' see \link[boot]{boot.ci}).}
#' }
#'
#' @param sup.stats a function, or a character naming a function to compute
#' some interest statistics from the \code{"msbm"} class input.
#' If \code{sup.args} is \code{NULL}, the function \code{sup.stats} must
#' take one argument and return a numeric vector. Otherwise, \code{sup.stats}
#' must take two arguments and return a numeric vector.
#'
#' @param sup.args (list of) supplemental argument(s) for the function
#' \code{sup.stats}. If supplied (and not \code{NULL}), \code{sup.args}
#' is passed as a second argument to \code{sup.stats} (an object of the
#' same shape as \code{object} being the first argument).
#'
#' @param fisher.matrix logical, should the Fisher (expected) information
#' matrix be used to obtain standard errors? Defaults to  \code{TRUE}.
#' The alternative, \code{fisher.matrix = FALSE}, is the observed
#' information matrix.
#'
#' @param sandwich logical, should the sandwich variance-covariance
#' estimator be considered for standard error calculation?
#' Defaults to \code{FALSE}.
#'
#' @param R2.method character, indicates the pseudo-\eqn{R^2} statistic
#' to include into the bootstrapped statistics. Only used when
#' \code{basic.stats} includes \code{"r.squared"}. See argument \code{method}
#' of \link[msbreg]{rsquared} for possible values of \code{R2.method}.
#'
#' @param cor.method character indicating which correlation coefficient
#' to consider (only used when \code{basic.stats} includes \code{"r.squared"},
#' and \code{R2.method = "COR"}).
#'
#' @param ... further named arguments passed to or from other methods.
#' Currently, the method takes and passes named arguments \code{strata},
#' \code{L}, \code{m}, and \code{weights} to \link[boot]{boot}
#' (see documentation therein).
#'
# @usage
# ## S3 method for class 'msbm'
# bootfit (object, basic.stats = NULL, sup.stats = NULL,
#          sup.args = NULL, fisher.matrix = TRUE,
#          R2.method = c("KL", "Nagelkerke", "COR"),
#          cor.method = "pearson", R = 999, seed = NULL,
#          sim = c("parametric", "ordinary", "balanced",
#                  "permutation", "antithetic"), stype = 'i',
#          simple = FALSE, parallel = c("no", "multicore", "snow"),
#          ncpus = getOption("boot.ncpus", 1L), cl = NULL, ...)
#
#' @details
#' This \link[msbreg]{bootfit} method wraps \link[boot]{boot} to specifically
#' handle \code{"msbm"} class objects. When \code{sim = 'parametric'}
#' (*parametric bootstrap*), the response variable in \code{object} is
#' resampled from a binomial distribution with success probability
#' determined by the estimated model parameters in \code{object}, while
#' re-using covariates in the original data used to fit the model.
#' For *non-parametric bootstrap*, both the response and the covariates are
#' resampled from the original data.
#'
#' The argument \code{sup.stats} is intended to allow more flexibility
#' in the specification of statistics: any quantity of interest can be
#' extracted from the fit \code{object} to obtain its bootstrap
#' distribution.
#'
#' On multicore plateformes, consider using arguments \code{parallel},
#' \code{ncpus}, and \code{cl} to allow parallel computations which can
#' be faster. This is currently based on the \code{R} package \code{parallel}
#' or \code{snow} (see \link[boot]{boot}).
#'
# \code{R = NULL} corresponds to jackknife approximation to the bootstrap,
#  that is, the usual *delete-one* estimates: one observation is deleted
#  from the original dataset at a time and the model parameters are
#  estimated using the remaining \code{nobs - 1} observations.
#'
#' @return An object inheriting from class \code{"boot"} (as defined
#' from package \code{boot}). See \link[msbreg]{bootfit}.
#'
#' @seealso \link[msbreg]{bootfit} for the default method.
#' Also see \link[msbreg]{MCmsbm} for Monte Carlo simulations
#' with possibly many multistage binomial models.
#'
#' @examples
#' ##* Load packages
#' require(boot)
#' require(msbreg)
#'
#' ##* Infertility data
#' data ("infert", package = "datasets")
#'
#' ## Logistic regression fit to the infert data
#' GLMres <- msbreg (case ~ spontaneous,
#'                   lambda.formula = ~ 0, # maximum success probability set to 1
#'                   data = infert, criterion = "ML")
#'
#' summary (GLMres)
#'
#' ## MSB model fit
#' MSBres <- msbreg (case ~ spontaneous,
#'                   lambda.formula = ~ 1, # unknown constant maximum success probability
#'                   start = c(GLMres$coefficients, 1),
#'                   data = infert, criterion = "ML")
#'
#' sumMSBres <- summary (MSBres)
#' sumMSBres
#'
#' AIC (GLMres, MSBres)
#' rsquared(GLMres, MSBres, adjust.size = TRUE)
#'
#' #* By the larger AIC (as compared with the simple GLM)
#' # and the not larger adjusted r-squared, the additional
#' # parameter lambda does not seem justified
#' #* Is adding 'lambda' overfitting?
#'
#' # Wald confidence interval (normality assumed at logit scale)
#' plogis(sumMSBres$coefficients[3,1] +
#'        c(-1, 1) * sumMSBres$coefficients[3,2] * qnorm(0.975))
#'
#' # This interval is useless (too wide)
#' # We could boostrap the fit
#' \dontrun{
#' set.seed(167)
#' bMSBres = bootfit(MSBres, R = 999) # (takes about 3 minutes)
#'
#' # Build bootstrap percentile confidence interval for lambda
#' boot.ci (bMSBres, index = 3, conf = 0.95, type = "perc",
#'          h = MSBres$link$linkinv)
#'
#' # Intervals :
#' # Level     Percentile
#' #   95%   [0.5321,  1.0000]
#'
#' # The data can only be used to infer that lambda is
#' # most likely above 0.5.
#'
#' # Estimating lambda is indeed overfitting,
#' # but reason advocates lambda = 1 against any
#' # other value above 0.5.
#'
#' # Based on a profiled likelihood graphic,
#' # we can use the ML estimate while keeping the variability in mind
#'  profile.logLike <- function(lambda) {
#'    f0 <- function(lambda0) {
#'      data <- infert
#'      data$eta.lambda <- qlogis(lambda0)
#'      res <- msbreg (case ~ spontaneous,
#'                     lambda.formula = ~ offset(eta.lambda) + 0,
#'                     start = GLMres$coefficients,
#'                     data = data, criterion = "ML",
#'                     frame = TRUE)
#'      return(logLik(res)[1])
#'    }
#'  sapply(lambda, f0)
#' }
#'
#' # (takes about 35 seconds)
#' curve(profile.logLike, from = 0.75, to = 0.95,
#'       n = 500, xlab = expression(lambda), ylab = "logLik")
#'
#' # Of course, we can look for additional
#' # predictor(s) in the 'infert' dataset
#' }
#'
#' # Simulated data example
#' data("test1data")
#'
#' # Fit a MSB model
#' msbfit = msbreg (cbind(y, Total - y) ~ x1 + offset(off1) | x2,
#'                  control = list(criterion = "MLJ"),
#'                  data = test1data)
#' summary(msbfit)
#'
#' \dontrun{
#' # Bootstrap the model fit (takes about 20 minutes)
#' bootmsbfit <- bootfit (msbfit, R = 999)
#' bootmsbfit
#'
#' # PARAMETRIC BOOTSTRAP
#' #
#' # Call:
#' # bootfit(object = msbfit, R = 999)
#' #
#' # Bootstrap Statistics :
#' #         original        bias     std. error
#' # t1*     2.739186  -0.003356267   0.11203355
#' # t2*    -1.831737   0.002217852   0.06162532
#' # t3*     3.273386   0.008053197   0.11695760
#' # t4*    -2.100589  -0.002623420   0.04072792
#' # t5*     1.857221   0.009574422   0.15732079
#' # t6*  -703.325885 -19.358176807  17.37753997
#' # t7*   428.224752  28.138879811  28.26254848
#' # t8* 20477.886114  38.158915179 275.78735080
#'
#' # Build bootstrap percentile confidence interval for the slope of 'x1'
#' boot.ci (bootmsbfit, index = 2, conf = 0.95, type = "perc")
#'
#' #
#' # BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#' # Based on 999 bootstrap replicates
#' #
#' # CALL :
#' # boot.ci(boot.out = bootmsbfit, conf = 0.95, type = "perc", index = 2)
#' #
#' # Intervals :
#' # Level     Percentile
#' # 95%   (-1.953, -1.719 )
#' # Calculations and Intervals on Original Scale
#'
#' # Build bootstrap percentile confidence intervals for 'lambda'
#' # (using the inverse link function to map logit(lambda) to lambda)
#' boot.ci (bootmsbfit, index = 5, conf = 0.95, type = "perc",
#'          h = msbfit$link$linkinv)
#'
#' # BOOTSTRAP CONFIDENCE INTERVAL CALCULATIONS
#' # Based on 999 bootstrap replicates
#' #
#' # CALL :
#' # boot.ci(boot.out = bootmsbfit, conf = 0.95, type = "perc", index = 5,
#' #         h = msbfit$link$linkinv)
#' #
#' # Intervals :
#' # Level     Percentile
#' # 95%   ( 0.8277,  0.9017 )
#' # Calculations and Intervals on  Transformed Scale
#' #
#' }
#'
# detach(package:msbreg)
# detach(package:boot)
#
# @exportS3Method bootfit msbm
#' @rdname bootfit.msbm
#' @exportS3Method msbreg::bootfit
# @export bootfit.msbm
bootfit.msbm <-
  function (object,
            coefficients = coefficients,
            basic.stats = NULL,
            sup.stats = NULL, sup.args = NULL,
            fisher.matrix = TRUE, sandwich = FALSE,
            R2.method = c("KL", "Nagelkerke", "COR"),
            cor.method = "pearson",
            R = 999, seed = NULL,
            sim = "parametric", stype = 'i',
            simple = FALSE, parallel = c("no", "multicore", "snow"),
            ncpus = getOption("boot.ncpus", 1L), cl = NULL,
             ...) {
    #* Matched call and computations
    mcall <- match.call()
    out <- bootfit.msbmcore (object = object, coefficients = coefficients,
                             basic.stats = basic.stats,
                             sup.stats = sup.stats, sup.args = sup.args,
                             fisher.matrix = fisher.matrix, sandwich = sandwich,
                             R2.method = R2.method, cor.method = cor.method,
                             R = R, seed = seed, sim = sim, stype = stype,
                             simple = simple, parallel = parallel, ncpus = ncpus,
                             cl = cl,...)

    # Update the 'call' component
    out$boot.call <- out$call
    out$call <- mcall

    return(out)
  }

### setMethod commented out to prevent the warning:
### Undocumented S4 methods:
###      generic 'bootfit' and siglist 'msbm'
#setMethod ("bootfit",
#           signature = "msbm",
#           definition = bootfit.msbm)

bootfit.msbmcore <-
  function (object,
            coefficients = object$coefficients,
            basic.stats = NULL,
            sup.stats = NULL, sup.args = NULL,
            fisher.matrix = TRUE, sandwich = FALSE,
            R2.method = c("KL", "Nagelkerke", "COR"),
            cor.method = "pearson",
            R = 999, seed = NULL,
            sim = c("parametric", "ordinary", "balanced",
                    "permutation", "antithetic"), stype = 'i',
            simple = FALSE, parallel = c("no", "multicore", "snow"),
            ncpus = getOption("boot.ncpus", 1L), cl = NULL, ...) {
    #* Save/Set the random generator seed
    eval(setsave.RNGstate())

    #* Basic statistics to be extracted from each fit
    eval(check.basic.stats())

    #* Matched arguments
    sim <- match.arg(sim)
    sim <- sim[1]
    stype <- match.arg(stype)
    R2.method <- match.arg(R2.method)
    cor.method <- match.arg(cor.method)

    #* Build the 'data' argument to 'boot::boot'
    original.data <- as.data.frame(object$data)
    nobs <- object$dims$nobs
    stopifnot(NROW(original.data) == nobs)
    if(is.null(object$frame))
      object$frame <- model.frame(object)
    original.data$y <- object$frame$y
    if (missing(coefficients)) {
      theta0 <- object$coefficients
    }
    else {
      if (identical(sim[1], "parametric")) {
        if (length(coefficients) == length(object$coefficients)) {
          if (!all(coefficients == object$coefficients)) {
            mmatrix <- model.matrix(object)
            if(NCOL(mmatrix) == length(coefficients)) {
              stop(paste0("argument 'coefficients'of wrong length: expected = ",
                          NCOL(mmatrix), ", observed = ", length(coefficients)))
            }
            object$fitted.values <- mmatrix %*% coefficients
          }
        }
        else {
          stop(paste0("argument 'coefficients'of wrong length: expected = ",
                      length(object$coefficients), ", observed = ", length(coefficients)))
        }
      }
      theta0 <- coefficients
    }
    link <- object$link
    control <- object$control
    method <- object$method
    envCall <- parent.frame()

    if (is.null(sup.stats)) {
      sup.stats <- function(x, ...) NULL
      nstat.values <- 0
    }
    else {
      if (is.character(sup.stats))
        sup.stats <- get(sup.stats, mode = "function", envir = parent.frame())
      if (!is.function(sup.stats))
        stop("the supplied 'sup.stats' is not a (name of a) function")

      if (is.null(sup.args))
        stat.values0 <- sup.stats (object)
      else
        stat.values0 <- sup.stats (object, sup.args)
      if (!is.numeric(stat.values0)) {
        stop("function 'sup.stats' must return a numeric vector")
      }
      nstat.values <- length(stat.values0)
      names.stat.values <- names(stat.values0)
      if(length(names.stat.values) != nstat.values) {
        names.stat.values <- paste0("stat.value", 1:nstat.values)
      }
    }

    #* A function to format different types of indices
    eval(make.get.indices ())

    #* A function to extract statistics from a fit
    eval(toget_get.stats())

    # Parameter names
    outnames <- names(theta0)
    if (length(outnames) != length(theta0)) {
      outnames <- paste0("coef.", 1:length(theta0))
    }
    secoef <- if (se.fit) paste0("se.", outnames)
    if (any(add.stats)) {
      sumobject <- summary(object)
      sout <- lapply(available.add.stats[add.stats],
                     FUN = function (stat) {
                       n1 <- length(sumobject[stat])
                       return(if (n1) paste0(stat, 1:n1))
                     })
      outnames <- c(outnames, unlist(sout))

      sesout <- if (se.fit) unlist(
        lapply(available.add.stats[add.stats],
               FUN = function (stat) {
                 n2 <- length(sumobject[paste0("se.", stat)])
                 return(if (n2) paste0("se.", stat, 1:n2))
               })
      )
    }
    outnames <- c(outnames, "logLik", "deviance", "null.deviance",
                  if (r.squared) R2.method,
                  if (nstat.values) names.stat.values)

    if (se.fit) {
      outnames <- c(outnames, secoef, if (any(add.stats)) sesout)
    }

    outnames <- make.unique(outnames)

    #* Build the 'statistic' and 'ran.gen' arguments to 'boot::boot'
    fitcall <- object$call


    #* Function 'ran.gen'
    if (identical(sim, "parametric")) {
      boot.statistic <- function (data, ...) {
        newcall <- fitcall
        newcall$data <- data
        newcall <- as.call(newcall)
        newfit <- eval(newcall, envir = envCall)
        get.stats (newfit)
      }

      ran.gen <- function(data, mle) {
        #mfcall <- object$frame$call
        #mfcall$data <- data
        #mf <- eval(mfcall, envir = envCall)
        data[[object$frame$auxy$ycolname]] <- as.numeric(simulate(object$frame,
                                                                  theta = mle,
                                                                  link = link))
        return(data)
      }
    }
    else {
      boot.statistic <- function (data, index, ...) {
        newcall <- fitcall
        newcall$data <- data[get.indices (index), ]
        newcall <- as.call(newcall)
        newfit <- eval(newcall, envir = envCall)
        get.stats (newfit)
      }

      ran.gen <- function(d, p) d # Default of the "boot::boot" function
    }

    #* Call the boot::boot function
    #simple <- simple & identical(sim, "ordinary")
    out <- boot::boot (data = original.data, statistic = boot.statistic,
                       R = R, sim = sim, stype = stype, ran.gen = ran.gen,
                       mle = theta0, simple = simple, ...,
                       parallel = parallel, ncpus = ncpus, cl = cl)

    if (length(outnames) == NCOL(out$t))
      colnames(out$t) <- outnames
    out$stats.names <- outnames

    class(out) <- c(class(out), "bootfit.msbm")

    out
  }

check.basic.stats <- function() {
  expression({
    #* Basic statistics to be extracted from each fit
    available.add.stats <- c("mu", "alpha", "lambda", "alpha.lambda")
    if (!is.null(basic.stats)) {
      stopifnot(is.character(basic.stats))
      stopifnot(all(basic.stats %in% c("r.squared",
                                       available.add.stats,
                                       "se.fit")))
    }
    add.stats <- available.add.stats %in% basic.stats
    r.squared <- "r.squared" %in% basic.stats
    se.fit <- "se.fit" %in% basic.stats
  })
}

toget_get.stats <- function () {
  expression({
    get.stats <- function (x) {
      sumx <- summary(x, fisher.matrix = fisher.matrix,
                      sandwich = sandwich, se.values = se.fit)
      out.stats <- c(sumx$coefficients[,1],
                     if (any(add.stats)) unlist(lapply(available.add.stats[add.stats],
                                                       FUN = function (stat) {
                                                         sumx[stat]
                                                       })),
                     logLik(x),
                     x$deviance,
                     x$null.deviance,
                     if (r.squared) unlist(rsquared(x, method = R2.method,
                                                    cor.method = cor.method)[1:length(R2.method)]),
                     if (nstat.values & is.null(sup.args)) sup.stats (x),
                     if (nstat.values & !is.null(sup.args)) sup.stats (x, sup.args)
      )

      if (se.fit) {
        out.stats <- c(out.stats,
                       sumx$coefficients[,2],
                       if (any(add.stats)) unlist(lapply(available.add.stats[add.stats],
                                                         FUN = function (stat) {
                                                           sumx[paste0("se.", stat)]
                                                         })))
      }

      return(out.stats)
    }
  })
}

make.get.indices <- function() {
  expression({
    switch(stype,
           i = {
             get.indices <- identity
           },
           f = {
             get.indices <- function(findex) {
               posindex <- which(findex > 0)
               index <- lapply(posindex,
                               FUN = function(j) {
                                 rep(j, length.out = findex[j])
                               })
               index <- unlist(index)

               return(index)
             }
           },
           w = {
             get.indices <- function(windex) {

               findex <- floor(nobs * windex / sum(windex))
               Delta.nobs <- nobs - sum(findex)
               if (Delta.nobs > 0) {
                 forder <- order(windex - findex/nobs, na.last = NA,
                                 decreasing = TRUE)
                 findex[forder[1:Delta.nobs]] <-
                   findex[forder[1:Delta.nobs]] + 1
               }

               posindex <- which(findex > 0)
               index <- lapply(posindex,
                               FUN = function(j) {
                                 rep(j, length.out = findex[j])
                               })
               index <- unlist(index)

               return(index)
             }
           })
  })
}
