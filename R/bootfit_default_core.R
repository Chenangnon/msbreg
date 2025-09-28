
bootfit_default_core <-
  function (object, basic.stats = NULL,
            sup.stats = NULL, sup.args = NULL, R = 999, seed = NULL,
            sim = c("ordinary", "parametric", "balanced", "permutation", "antithetic"),
            stype = c("i", "f", "w"), yname, simple = FALSE,
            parallel = c("no", "multicore", "snow"),
            ncpus = getOption("boot.ncpus", 1L),
            cl = NULL, ...,
            .force.coef = TRUE,
            envCall = parent.frame(n=2)) {
    #* Save/Set the random generator seed
    eval(setsave.RNGstate())

    #* Basic statistics to be extracted from each fit
    eval(check.def.basic.stats())

    # The call to bootstrap
    fitcall <- getCall(object)
    if (is.null(fitcall)) {
      stop("no 'call' component found in 'object'")
    }
    original.data <- object$data
    if (is.null(original.data))
      original.data <- fitcall$data
    if (is.null(original.data)) {
      stop("no 'data' component found in 'object' and 'object$call'")
    }
    original.data <- as.data.frame(original.data)
    nobs <- NROW(original.data)
    simple <- as.logical(simple)
    stopifnot(is.logical(simple))

    #* Possible data simulation (resampling) methods
    sim.ok <- c("parametric", "ordinary", "balanced", "permutation", "antithetic")
    stopifnot(all(sim %in% sim.ok))
    sim <- sim[1]
    stype <- match.arg(stype)

    #* Check if there is a simulate method for object
    if (identical(sim, "parametric")) {
      has.simulate.method <- check.has.simulate.method (object)
      if (!has.simulate.method) {
        warning("sim = 'parametric' only possible when 'object' is of a class with a 'simulate' method ")
        stop(paste0("cannot find a 'simulate' method for objects of class ", class(object), "."))
      }

      if (missing(yname)) {
        yname <- catch.conditions({
          as.character(object$call$formula)[[2]]
        })$value

        if (any(class(yname) %in% c("simpleError", "error", "condition", "try-error"))) {

          yvec <- catch.conditions({
            original.data[["y"]]
          })$value

          if (any(class(yvec) %in% c("simpleError", "error", "condition", "try-error"))) {
            stop("supply the argument 'yname'")
          }

          if (length(yvec) != nobs) {
            stop("supply the argument 'yname'")
          }
        }
        else {
          yname <- get.ycolname (yname)

          yvec <- catch.conditions({
            original.data[[yname]]
          })$value

          if (any(class(yvec) %in% c("simpleError", "error", "condition", "try-error"))) {
            print(yvec)
            stop("supply the argument 'yname'")
          }

          if (length(yvec) != nobs) {
            stop("supply the argument 'yname'")
          }
        }
      }
      else {
        yvec <- catch.conditions({
          original.data[[yname]]
        })$value

        if (any(class(yvec) %in% c("simpleError", "error", "condition", "try-error"))) {
          print(yvec)
          stop("incorrect 'yname' argument")
        }

        if (length(yvec) != nobs) {
          stop("incorrect 'yname' argument")
        }
      }
    }

    #* A function to format different types of indices
    eval(make.get.indices ())

    #* Basic statistics to be bootstrapped
    eval(toget_get.def.basic.stats())

    #* Additional statistics to be bootstrapped
    eval(toget_get.def.sup.stats())

    #* 'get.stats' function
    if (is.null(get.basic.stats) & is.null(get.sup.stats)) {
      stop("no statistic to bootstrap")
    }
    else if (is.null(get.sup.stats)) {
      get.stats <- get.basic.stats
      outnames <- basic.names
    }
    else if (is.null(get.basic.stats)) {
      get.stats <- get.sup.stats
      outnames <- sup.names
    }
    else {
      get.stats <- function(x, ...) {
        c(get.basic.stats(x, ...), get.sup.stats(x, ...))
      }
      outnames <- c(basic.names, sup.names)
    }

    #* Statistic function for boot::boot
    if (identical(sim, "parametric")) {
      statistic <- function (data, index) {
        newcall <- fitcall
        newcall$data <- data[index, ]
        newcall <- as.call(newcall)
        newfit <- eval(newcall, envir = envCall)
        get.stats (newfit)
      }

      ran.gen <- function(data, mle) {
        data[[yname]] <- unlist(simulate(object, nsim = 1))
        return(data)
      }
    }
    else {
      statistic <- function (data, index, ...) {
        newcall <- fitcall
        newcall$data <- data[get.indices (index), ]
        newcall <- as.call(newcall)
        newfit <- eval(newcall, envir = envCall)
        get.stats (newfit)
      }

      ran.gen <- function(d, p) d # Default of the "boot::boot" function
    }

    #* Call the boot::boot function
    simple <- simple & identical(sim, "ordinary")
    out <- boot::boot (data = original.data,
                       statistic = statistic,
                       R = R,
                       sim = sim, stype = stype,
                       ran.gen = ran.gen,
                       mle = 0,
                       simple = simple, ...,
                       parallel = parallel,
                       ncpus = ncpus, cl = cl)

    if (length(outnames) == NCOL(out$t))
      colnames(out$t) <- outnames
    out$stats.names <- outnames

    out
  }

# stype = 'f' means the usual bootstrap but with frequency representing
# the number of times each unit is drawn.
# stype = 'w' means the relative frequencies (weights) of each element of the
# initial data in the bootstrap sample
# See: https://stats.stackexchange.com/questions/606385/how-to-use-the-argument-stype-in-boot-package-in-r

check.def.basic.stats <- function() {
  expression({
    #* Basic statistics to be extracted from each fit
    available.basic.stats <- c("all", "coef", "logLik", "deviance", "null.deviance")
    if (!is.null(basic.stats)) {
      stopifnot(is.character(basic.stats))
      basic.stats[basic.stats %in% "coefficients"] <- "coef"
      stopifnot(all(basic.stats %in% available.basic.stats))
      if ("all" %in% basic.stats) {
        basic.stats <- available.basic.stats[-1]
      }
    }
    else if (is.null(sup.stats)) {
      basic.stats <- c("coef", "logLik", "deviance", "null.deviance")
    }
    selected.basic.stats <- available.basic.stats[which(available.basic.stats %in% basic.stats)]
  })
}

toget_get.def.basic.stats <-
  function() {
    expression({
      if (!is.null(basic.stats) | is.null(sup.stats)) {
        theta0 <- object$coefficients
        has.coefs <- FALSE
        if ("coef" %in% selected.basic.stats | .force.coef) {
          if(is.null(theta0) | !is.numeric(theta0)) {
            theta0 <- catch.conditions({
              stats::coef(object)
            })$value
            if(is.null(theta0)) {
              if (is.null(sup.stats))
                stop("'object' must have a 'coefficients' component or a 'coef' method when 'sup.stats = NULL'")

              coefvalues <- function(x) NULL
            }
            else if (!is.numeric(theta0)) {
              if (is.null(sup.stats)) {
                print(theta0)
                stop ("'coef(object)' must return a numeric vector of model coefficients when 'sup.stats = NULL'")
              }

              coefvalues <- function(x) NULL
            }
            else {

              coefvalues <- function(x) {
                stats::coef(x)
              }
              has.coefs <- TRUE
            }
          }
          else {
            coefvalues <- function(x) {
              x$coefficients
            }
            has.coefs <- TRUE
          }
        }
        else {
          coefvalues <- function(x) NULL
        }

        if ("logLik" %in% selected.basic.stats) {
          Llik0 <- object$logLik[1]
          if(is.null(Llik0) | !is.numeric(Llik0)) {
            Llik0 <- catch.conditions({
              stats::logLik(object)[1]
            })$value
            if(is.null(Llik0) | !is.numeric(Llik0)) {
              Likvalues <- function(x) {
                NULL
              }
            }
            else {
              Likvalues <- function(x) {
                stats::logLik(x)[1]
              }
            }
          }
          else {
            Likvalues <- function(x) {
              x$logLik[1]
            }
          }
        }
        else {
          Likvalues <- function(x) NULL
        }

        if ("deviance" %in% selected.basic.stats) {
          Dev0 <- object$deviance[1]
          if(is.null(Dev0) | !is.numeric(Dev0)) {
            Dev0 <- catch.conditions({
              stats::deviance(object)[1]
            })$value
            if(is.null(Dev0) | !is.numeric(Dev0)) {
              Devvalues <- function(x) NULL
            }
            else {
              Devvalues <- function(x) {
                stats::deviance(x)[1]
              }
            }
          }
          else {
            Devvalues <- function(x) {
              x$deviance[1]
            }
          }
        }
        else {
          Devvalues <- function(x) NULL
        }

        if ("null.deviance" %in% selected.basic.stats) {
          nullDev0 <- object$null.deviance
          if (is.null(nullDev0) | !is.numeric(nullDev0)) {
            nullDevvalues <- function(x) NULL
          }
          else {
            nullDevvalues <- function(x) {
              x$null.deviance[1]
            }
          }
        }
        else {
          nullDevvalues <- function(x) NULL
        }

        get.basic.stats <- function(x, ...) {
          c(coefvalues(x),
            logLik = Likvalues(x),
            deviance = Devvalues(x),
            null.deviance = nullDevvalues (x))
        }

        # Parameter names
        if (has.coefs) {
          basic.names <- names(theta0)
          if (length(basic.names) != length(theta0)) {
            basic.names <- paste0("coef.", 1:length(theta0))
          }
        }
        else {
          basic.names <- NULL
        }

        if ("logLik" %in% selected.basic.stats) {
          if (!is.null(Llik0) & is.numeric(Llik0)) {
            basic.names <- c(basic.names, "logLik")
          }
        }

        if ("deviance" %in% selected.basic.stats) {
          if (!is.null(Dev0) & is.numeric(Dev0)) {
            basic.names <- c(basic.names, "deviance")
          }
        }

        if ("null.deviance" %in% selected.basic.stats) {
          if (!is.null(nullDev0) & is.numeric(nullDev0)) {
            basic.names <- c(basic.names, "null.deviance")
          }
        }
      }
      else {
        get.basic.stats <- NULL
        basic.names <- NULL
      }
    })
  }

toget_get.def.sup.stats <-
  function() {
    expression({
      if (!is.null(sup.stats)) {
        if (is.character(sup.stats))
          sup.stats <- get(sup.stats, mode = "function", envir = parent.frame())
        if (!is.function(sup.stats)) {
          print(sup.stats)
          stop("'sup.stats' function not found")
        }

        if (is.null(sup.args)) {
          stats0 <- sup.stats (object)
          if (is.list(stats0)) {
            get.sup.stats <- function(x, ...) {
              unlist(sup.stats(x))
            }
          }
          else if (is.numeric(stats0)) {
            get.sup.stats <- function(x, ...) {
              sup.stats (x)
            }
          }
          else {
            print(stats0)
            stop("'sup.stats()' output not recognized ('list' or 'numeric' expected)")
          }
        }
        else {
          stats0 <- sup.stats (object, sup.args)

          if (is.list(stats0)) {
            get.sup.stats <- function(x, ...) {
              unlist(sup.stats(x, sup.args))
            }
          }
          else if (is.numeric(stats0)) {
            get.sup.stats <- function(x, ...) {
              sup.stats(x, sup.args)
            }
          }
          else {
            print(stats0)
            stop("'sup.stats(, sup.args)' output not recognized ('list' or 'numeric' expected)")
          }
        }

        # Parameter names
        sup.names <- names(stats0)
        if (length(sup.names) != length(stats0)) {
          sup.names <- paste0("stat.", 1:length(stats0))
        }

      }
      else {
        get.sup.stats <- NULL
        sup.names <- NULL
      }
    })
  }

check.has.simulate.method <-
  function(object) {
    has.simulate.method <- FALSE

    sim0 <- catch.conditions({
      stats::simulate(object, nsim = 1)
    })$value

    if (!any(class(sim0) %in% c("simpleError", "error",
                                "condition", "try-error"))) {
      ftd <- fitted(object)
      n <- length(ftd)
      if (is.data.frame(sim0)) {
        if (all(dim(sim0) == c(n, 1))) {
          has.simulate.method <- TRUE
        }
      }
      else if (is.list(sim0)) {
        if (length(sim0) == 1) {
          if (length(sim0[[1]]) == n) {
            has.simulate.method <- TRUE
          }
        }
      }
    }

    return(has.simulate.method)
  }
