#
#  See: https://www.r-bloggers.com/2016/11/why-optim-is-out-of-date/
#'
#' Self-starting Multivariable Optimization
#'
#' Generate points from a specified design to start
#' a multivariable optimization with \link[stats]{nlminb} or \link[stats]{optim}.
#' Starting values can be requested from a \code{Central-Composite Design}
#' (CCD) or a \code{Box-Behnken Design} (BBD) over the sample space.
#'
#' @param npars numeric (scalar), number of parameters taken by the function
#' to be optimized (\code{fn}).
# \code{npars}\eqn{\geq 2} is required.
#'
#' @param inf numeric (scalar), limits of the sample space to substitute to
#' any \code{Inf} in \code{lower} or \code{upper}. This substitution is only
#' used to generate the starting values, the supplied \code{lower} or
#' \code{upper} arguments are passed unchanged to \link[stats]{nlminb}
#' or \link[stats]{optim}.
#'
#' @param design character, one of \code{"bbd"} (use all points from a BBD
#' as starting parameter values), \code{"ccd"} (use all points from a CCD
#' as starting values), and \code{"one"} (use one starting point: the
#' center of the sample space).
#'
#' @param pars optional matrix of initial values (rows of \code{pars})
#' for the parameters to be optimized over. If \code{pars} is supplied,
#' all previous arguments are ignored.
#'
#' @param fn,gr,...,method,lower,upper,control,hessian arguments passed
#' to \link[stats]{nlminb} or \link[stats]{optim}.
#'
#' @param chunk.size,cl arguments passed to \link[parallel]{parLapply}.
#' @param LB logical, should \link[parallel]{parLapplyLB} be used instead
#' of \link[parallel]{parLapply}?
#'
#' @param star.points character indicating how axial points
#' (start points) are introduced in the design. One of \code{"rotatable"},
#' \code{"orthogonal"}, \code{"spherical"}.
#'
#' @details
#' Note that arguments after ... must be matched exactly.
#'
#' The function aims to self-start the General-purpose Optimization
#' routine \link[stats]{nlminb} or \link[stats]{optim}. This is mostly
#' useful for non regular/smooth objective function \code{fn} with the
#' simplex method \code{method = "Nelder-Mead"}. Trading-off time for
#' precision, \code{optim.selfstart} increases the chance of finding
#' the global optimum of the objective by starting from a wide range of
#' points in the sample space, depending on the argument \code{design}.
#'
#' The argument \code{pars} is the matrix counterpart of the
#' argument \code{par} of \link[stats]{optim} to allow more flexibility:
#' a user-determined set of starting points can be supplied.
#'
#' In all cases, \link[stats]{nlminb} or \link[stats]{optim} is
#' repeatedly called with each starting value, and the best optimization
#' result is returned.
#'
#' On multicore plateformes, consider using argument \code{chunk.size},
#' \code{cl} and \code{LB} o allow parallel computations which can be faster.
#' This is currently based on the \code{R} package \code{parallel}.
#'
#' @return a list as returned by \link[stats]{optim}, with the following
#' additional elements:
# \describe{
#'
#' \item{\code{starts}}{ matrix of used starting parameter values;}
#'
#' \item{\code{table}}{ matrix of optimization results for each starting
#' point: the objective function value (first column), the optimal vector
#' of parameters, and the numbers of function and gradient evaluations;}
#'
#' \item{\code{k}}{ the number of the row (in \code{starts}) of an
#' optimal starting parameter, that is, starting parameter value which
#' led to the optimal objective function value.}
# }
#'
#' @seealso \link[stats]{optim} for the workhorse routine.
#'
#' @export optim.selfstart
#'
#' @importFrom stats optim
#' @importFrom stats optimHess
#' @importFrom parallel getDefaultCluster
#'
#' @examples
#' require(graphics)
#'
#'
#' #* Example taken from optim help page:
#' #  Modified Rosenbrock Banana function
#' fr <- function(x) {
#'  x1 <- x[1]
#'  x2 <- x[2]
#'  100 * (x2 - floor(x1)*x1)^2 + (1 - floor(x1))^2
#' }
#'
#' # The minimum value of fr is zero, attained at x = c(a, a)
#' # for any a in the semi-open [1, 2).
#'
#' # As expected derivative based method fails to find the minimum value 0
#' # starting from x = c(-1.2, 1) as on optim help page
#' optim(par = c(-1.2,1), fn = fr, method = "BFGS")[c('value', 'par')]
#'
#'  # But the simplex method also fails
#' optim(par = c(-1.2,1), fn = fr, method = "Nelder-Mead")[c('value', 'par')]
#'
#' # Try multiple starting values with the simplex method
#' optim.selfstart(2, fn = fr, method = "Nelder-Mead",
#'                 design = 'ccd')[c('value', 'par')]
#'
#'
optim.selfstart <- function (npars = NCOL(pars),
                             design = c("bbd", "ccd", "one"),
                             star.points = c("rotatable",
                                             "orthogonal",
                                             "spherical"),
                             inf = 7, # qlogis(.999)
                             pars = NULL,
                             fn, gr = NULL, ...,
                             method = c("Nelder-Mead", "CG",
                                        "BFGS", "L-BFGS-B",
                                        "SANN", "Brent",
                                        "port", "nlminb"),
                             lower = -Inf, upper = Inf,
                             control = list(),
                             hessian = FALSE,
                             chunk.size = NULL,
                             cl = parallel::getDefaultCluster(),
                             LB = FALSE) {
  mcall <- match.call()
  if(missing(pars) & missing(npars)) {
    stop("at least one of arguments 'pars' and 'npars' must be supplied")
  }

  if (is.null(pars)) {
    stopifnot(is.numeric(npars), is.finite(npars))
    stopifnot(floor(npars) == npars)
    stopifnot(is.numeric(inf), is.finite(inf))
    inf <- abs(inf)
    eps <- control$epsilon
    if (is.null(eps))
      eps <- sqrt(.Machine$double.eps)

    lows <- lower
    lows[lows == -Inf] <- min(-inf, lows[is.finite(lows)])
    stopifnot(is.finite(lows))
    lows <- lows + eps
    lows <- rep(lows, length.out = npars)
    upps <- upper
    upps[upps == Inf] <- max(inf, upps[is.finite(upps)])
    stopifnot(is.finite(upps))
    upps <- upps - eps
    upps <- rep(upps, length.out = npars)

    design <- match.arg(design)
    switch(design,
           one = {
             starts <- rbind((lows + upps) / 2)
           },
           ccd = {
             star.points <- match.arg(star.points)
             starts <- ccd.start (npars = npars,
                                  lows = lows,
                                  upps = upps,
                                  star.points = star.points)
           },
           bbd = {
             starts <- bbd.start (npars = npars,
                                  lows = lows,
                                  upps = upps)
           })
  }
  else {
    starts <- as.matrix(pars)
    stopifnot(NCOL(starts) == npars)
    stopifnot(is.numeric(starts))
  }

  method <- match.arg(method)
  if (identical(tolower(method), "port"))
    method <- "nlminb"
  switch(method,
         nlminb = {
           nlmcontrol <- list(iter.max = if (!is.null(control$maxit)) control$maxit,
                           trace = if (!is.null(control$trace)) control$trace,
                           abs.tol = if (!is.null(control$abstol)) control$abstol,
                           rel.tol = if (!is.null(control$reltol)) control$reltol)

           if (!is.null(control$fnscale)) {
             if (control$fnscale < 0) {
               fn0 <- fn
               fn <- function(x, ...) {
                 - fn0 (x, ...)
               }
             }
           }

           infunc <- function (k) {
             outk <- catch.conditions({
               stats::nlminb(par = starts[k,],
                             objective = fn, gradient = gr, ...,
                             lower = lower,
                             upper = upper,
                             control = nlmcontrol)
             })$value

             if (!any(class(outk) %in% c("simpleError", "error", "condition", "try-error"))) {
               outk$value <- outk$objective
               if (!is.null(control$fnscale)) {
                 if (control$fnscale < 0) {
                   outk$value <- - outk$value
                 }
               }
               outk$objective <- NULL
               outk$counts <- outk$evaluations
               outk$evaluations <- NULL
             }

             outk
           }
         },
         {
           infunc <- function (k) {
             catch.conditions({
               stats::optim(par = starts[k,],
                     fn = fn, gr = gr, ..., method = method,
                     lower = lower,
                     upper = upper,
                     control = control,
                     hessian = FALSE)
             })$value
           }
         })

  list.optim <- matteLapply(1:NROW(starts),
                            FUN = infunc,
                            chunk.size = chunk.size,
                            cl = cl,
                            LB = LB)
  optvalues <- lapply(list.optim, FUN = function(x) {
    if (any(class(x) %in% c("simpleError", "error", "condition", "try-error"))) {
      NA
    }
    else {
      c(value = x$value, x$par, x$counts)
    }
  })
  optvalues <- do.call("rbind", optvalues)

  if (all(is.na(optvalues[,1]))) {
    stop(paste0("all '", design, "' starting points failled"))
  }

  if (NROW(starts) > 1) {
  if (!is.null(control$fnscale)) {
    if (control$fnscale < 0)
      k <- which.max (optvalues[,1])[1]
    else
      k <- which.min (optvalues[,1])[1]
  }
  else {
    k <- which.min(optvalues[,1])[1]
  }
  }
  else
    k <- 1

  switch(method,
         nlminb = {
           main <- stats::nlminb(par = optvalues[k, 2:(npars + 1)],
                                 objective = fn, gradient = gr, ...,
                                 lower = lower,
                                 upper = upper,
                                 control = nlmcontrol)
           main$value <- main$objective
           main$objective <- NULL
           main$counts <- main$evaluations
           main$evaluations <- NULL
           if (hessian) {
             main$hessian <- stats::optimHess (par = main$par,
                                               fn = fn,
                                               gr = gr, ...,
                                               control = control)
           }
           if (!is.null(control$fnscale)) {
             if (control$fnscale < 0) {
               main$value <- - main$value
               if (hessian) {
                 main$hessian <- - main$hessian
               }
             }
           }
         },
         {
           main <- stats::optim(par = optvalues[k, 2:(npars + 1)],
                                fn = fn, gr = gr, method = method,
                                lower = lower,
                                upper = upper,
                                control = control,
                                hessian = hessian)
         })

  main$counts <- main$counts + list.optim[[k]]$counts

  main$table <- optvalues
  main$starts <- starts
  main$k <- k

  structure(main, call = mcall)

}
# @param eps positive real \eqn{\epsilon}, tolerance for reaching zero.

ccd.start <- function (npars,
                       lows = -10, upps = 10,
                       factorial = TRUE,
                       center = TRUE,
                       axial = TRUE,
                       star.points = c("rotatable",
                                       "orthogonal",
                                       "spherical")) {
  stopifnot(is.numeric(npars))
  stopifnot(is.finite(npars))
  stopifnot(floor(npars) == npars)

  if (is.null(center))
    center <- !missing(lows) | !missing(upps)

  if (is.null(axial)) {
    axial <- TRUE
  }

  if (is.null(factorial)) {
    factorial <- (npars <= 4) | !any(c(center, axial))
  }

  if (!any(c(factorial, center, axial))) {
    stop("at leat one of 'factorial', 'center', and 'axial' must be TRUE")
  }

  # Factorial points in the design
  if (factorial) {
    bfactorial <- lapply(1:npars, FUN = function(k) c(-1, 1))
    bfactorial <- as.matrix(do.call("expand.grid", bfactorial))
    rownames(bfactorial) <- paste0("Factorial.", 1:(2^npars))
    colnames(bfactorial) <- paste0("par", 1:npars)
  }
  else
    bfactorial <- NULL

  # Central points in the design
  if (center) {
    bcenter <- rbind(numeric(npars))
    rownames(bcenter) <- "Center"
    colnames(bcenter) <- paste0("par", 1:npars)
  }
  else
    bcenter <- NULL

  # Axial points in the design
  if (axial) {
  star.points <- match.arg(star.points)
  switch(star.points,
         orthogonal = {
           f <- 2^npars
           t <- 2 * npars + 1
           Q <- (sqrt(f+t) - sqrt(f))^2
           alpha <- (Q * f / 4)^(0.25)
         },
         rotatable = {
           alpha <- 2^(npars/4)
         },
         spherical = {
           alpha <- sqrt(npars)
         })

  if (npars > 1) {
    baxial <- alpha * diag(npars)
    baxial <- rbind(baxial, - baxial)[c(seq(1, 2 * npars - 1, 2),
                                        seq(2, 2 * npars, 2)), , drop = FALSE]
  }
  else {
    baxial <- rbind(alpha,
                    -alpha)
  }

  rownames(baxial) <- paste0("Axial.", 1:(2*npars))
  colnames(baxial) <- paste0("par", 1:npars)
  }
  else
    baxial <- NULL

  # The full design
  ccdesign <- rbind(bcenter,
                  bfactorial,
                  baxial)
  ccdesign <- unique(ccdesign) / alpha

  # Bounds
  Lows <- rep(lows, length.out = npars)
  upps <- rep(upps, length.out = npars)
  lows <- pmin(Lows, upps)
  upps <- pmax(Lows, upps)
  units <- (upps - lows)
  ccdesign <- t(lows + ((t(ccdesign) + 1) / 2) * units)

  structure(ccdesign,
            npars = npars,
            star.points = star.points,
            alpha = alpha,
            class = "ccd.start")
}

bbd.start <- function (npars,
                       lows = -10, upps = 10) {
  stopifnot(is.numeric(npars))
  stopifnot(is.finite(npars))
  stopifnot(floor(npars) == npars)

  bfactorial2 <- matrix(c(-1, -1,
                         1, -1,
                         -1, 1,
                         1, 1), nrow = 4, ncol = 2,
                        byrow = TRUE)

  if (npars > 1) {
  index <- do.call("rbind",
                   lapply(1:(npars - 1),
                          FUN = function(k) {
                            cbind(k, (k + 1):npars)
                          }))


  nblocks <- choose(npars, 2)
  bbdesign <- do.call("rbind",
                     lapply(1:nblocks,
                            FUN = function(k) {
                              bk <- matrix(0, nrow = 4, ncol = npars)
                              bk[, index[k,]] <- bfactorial2

                              bk
                            }))
  }
  else {
    bbdesign <- rbind(-1, 1)
  }

  bbdesign <- rbind(bbdesign,
                    0)
  colnames(bbdesign) <- paste0("par", 1:npars)

  # Bounds
  Lows <- rep(lows, length.out = npars)
  upps <- rep(upps, length.out = npars)
  lows <- pmin(Lows, upps)
  upps <- pmax(Lows, upps)
  units <- (upps - lows)
  bbdesign <- t(lows + ((t(bbdesign) + 1) / 2) * units)

  structure(bbdesign,
            npars = npars,
            class = "bbd.start")
}
