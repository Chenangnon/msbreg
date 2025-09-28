#' Auxiliary for Controlling Binomial Model Fitting
#'
#' Auxiliary function for fitting a binomial model.
#' Typically used internally by \link[msbreg]{msbm.fit}
#' for fitting Multistage Binomial (MSB) Models.
#' It may also be used to construct the control
#' argument for \link[msbreg]{msbreg} which wraps
#' \link[msbreg]{msbm.fit}.
#'
#' @param criterion character indicating the objective function(s) to be
#' optimized to find the estimates of model parameters.
#' Currently recognized criteria include:
#'
#' \itemize{
#' \item \code{'ML'} for Maximum Likelihood (ML) estimation,
#' \item \code{'MLJ'} for ML estimation with Jeffreys' prior,
#' \item \code{'MLPJ'} for ML estimation with partial Jeffreys' prior (penalizing only slopes), and
#' \item \code{'MLJIC'} for MLJ estimation with intercept-correction.
#' }
#'
#'
#' Defaults to \code{'MLJ'}. The \code{'MLJIC'} criterion is the
#' \code{'MLJ'} followed with a post-hoc adjustment of intercepts
#' (if any, they are re-estimated by "ML") so that predicted probabilities
#' are unbiased to the first order (see the method \code{FLIC} in
#' \insertCite{puhr2017firth;textual}{msbreg}). The \code{'MLPJ'} criterion
#' is based on the same argument but avoid the post-hoc adjustment by
#' directly excluding intercepts from the information matrix used for
#' Jeffreys prior calculation \insertCite{xu2023penalized}{msbreg}.
#'
#' @param slope.signs NULL, or numeric vector indicating if regression
#' slopes in the MSB model should be forced to have some specific signs.
#' The default \code{slope.signs = 0} (equivalent to
#' \code{slope.signs = NULL}) means no sign forcing. The specification
#' \code{slope.signs = -1} or \code{slope.signs = 1} indicates
#' restriction to the negative or positive orphan real domain,
#' respectively. A scalar input \code{slope.signs} is recycled to
#' the number of slopes in the MSB regression model.
#'
#' @param unit.lambda either \code{NULL}, or a logical value indicating
#' whether a unit asymptote model (\eqn{\lambda = 1}) should be
#' explored as part of the parameter search procedure.
#'
#' The argument \code{unit.lambda} is **ignored if there is no intercept in the**
#' **linear predictor** \eqn{\eta_{\lambda}} of the \eqn{\lambda} model component
#' (in this case, \code{unit.lambda = NULL} means \code{unit.lambda = FALSE},
#' but setting \code{unit.lambda = TRUE} will not have any impact on fitting).
#' This includes when \eqn{\eta_{\lambda}} only has an offset term (then
#' \code{unit.lambda} is ignored).
#'
#' If \eqn{\eta_{\lambda}} includes an intercept, then setting the
#' intercept to \eqn{+\infty} (while all slopes in \eqn{\eta_{\lambda}}, if any,
#' are set to zero) gives the unit asymptote model (with \eqn{\eta_{\lambda} = +\infty}
#' and \eqn{\lambda = 1}).
#' The later is a limiting model which also happens to have a singular
#' information matrix. Thus, numerical optimization can generally only approach it
#' but not reach it (unless there is no offset or predictor, and \eqn{\lambda} is
#' an unknown constant: in this case we can re-parameterize the model directly
#' in terms of \eqn{\lambda \in (0, 1]} with finite bounds).
#'
#' When \eqn{\eta_{\lambda}} includes an intercept, the default behavior
#' (corresponding to \code{unit.lambda = NULL}) depends then on whether the
#' linear predictor also includes offsets or not.
#'
#' \itemize{
#' \item When the MSB model includes offsets in \eqn{\eta_{\lambda}} (in addition
#' to an intercept, and slopes, if any), \code{unit.lambda = NULL} means
#' \code{unit.lambda = FALSE}.
#'
#' This is because when \eqn{\eta_{\lambda} = +\infty}, any offset in the
#' linear predictor \eqn{\eta_{\lambda}} would be ignored.
#' Of course, the user can always set \code{unit.lambda = TRUE} to force
#' the parameter search to explore \eqn{\eta_{\lambda} = +\infty},
#' but beware of the implication that offsets
#' are then ignored during this peculiar search. If the final estimate of
#' \eqn{\lambda} is one, then offsets are not part of the fitted model.
#'
#' \item When the MSB model does not include any offset in \eqn{\eta_{\lambda}}
#' (but has an intercept, and slopes for predictors, if any),
#' \code{unit.lambda = NULL} means \code{unit.lambda = TRUE}.
#' }
#'
#' @param method a character indicating the optimizer to consider
#' for the primary search of parameter estimates.
#' Currently, this can be one of
#'
#' \itemize{
#' \item \code{"BFGS"} (the default) (replaced by \code{"L-BFGS-B"} when
#' \code{slope.signs} is not \code{NULL}, nor \code{0}),
#' \item \code{"CG"};
#' \item \code{"nlminb"} and
#' \item \code{"Nelder-Mead"}.
#' }
#'
#' If \code{method = "CG"} is specified but \code{slope.signs} is not
#' \code{NULL}, then \code{method = "L-BFGS-B"} (which can handle
#' box constraints) is returned instead.
#'
#' @param start character, should a batch of starting parameter values
#' be used? Defaults to \code{start = 'none'} which indicates no batch
#' of starting values. The alternatives are \code{start = 'one'} for
#' one starting value (central point of the search space);
#' \code{start = 'ccd'} for starting values generated using a
#' Central-Composite Design and \code{start = 'bbd'} for
#' starting values generated using a Box-Behnken Design.
#'
#' @param star.points,inf,pars arguments that are all passed to
#' \link[msbreg]{optim.selfstart} when \code{start} is not \code{'none'}
#' or \code{pars} is non-\code{NULL}. In particular, \code{pars} can be
#' used to supply a set of suspected good starting values of model parameters.
#' Note that a \code{pars} input with a number of columns different from the
#' number of model parameters gives an error.
#'
#' @param epsilon positive convergence tolerance, typically used as
#' relative tolerance (e.g. \code{reltol} for \link[stats]{optim}).
#'
#' @param maxit The maximum number of iterations.
#' Defaults to \code{1000}.
#'
#' @param trace logical or non-negative integer, indicating if tracing
#' information on the progress of the fitting procedure should be produced.
#' Higher values may produce more tracing information.
#' Defaults to \code{FALSE}.
#'
#' @param gradient logical, should the gradient of the fitting function be
#' returned as part of the results?
#'
#' @details
#' This function is similar to the \link[stats]{glm.control} function for
#' fitting \link[stats]{glm}'s. The control argument of \link[msbreg]{msbreg}
#' is by default passed to the control argument of \link[msbreg]{msbm.fit}.
#' The latter uses its elements as arguments to \code{msbm.control}
#' which provides defaults and sanity checking.
#'
#' @return A list with components named as the arguments.
#'
#' @export msbm.control
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \link[msbreg]{msbm.fit},
#' the fitting procedure used by \link[msbreg]{msbreg}.
#'
#' @examples
#' library(msbreg)
#' data('test1data', package = 'msbreg')
#'
#' MSBfit0 <- msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                    lambda.formula = ~ 1,
#'                    weights = Total,
#'                    data = test1data)
#'
#' # Reducing 'epsilon' value to 1e-14 (default is 1e-08)
#' MSBfit1 <- msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                    lambda.formula = ~ 1,
#'                    weights = Total,
#'                    data = test1data,
#'                    epsilon = 1e-14)
#'
#' # Slightly lower criterion value for epsilon = 1e-14
#' options(digits = 7)
#' deviance(MSBfit0)
#' deviance(MSBfit1)
#'
#' # For epsilon = 1e-14, we get slightly larger parameter values for stage 1
#' # and slightly smaller parameter values for stage 2
#' summary (MSBfit0)
#' summary (MSBfit1)
#'
#' # Changing fitting criterion to 'ML' produces a larger asymptote probability
#' MSBfit <- msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                    lambda.formula = ~ 1,
#'                    weights = Total,
#'                    data = test1data,
#'                    criterion = 'ML')
#'
#' summary (MSBfit)
#'
msbm.control <- function (criterion = "MLJ",
                          slope.signs = 0,
                          unit.lambda = NULL,
                          method = "BFGS",
                          start = c("none", "one", "bbd", "ccd"),
                          star.points = c("rotatable", "orthogonal", "spherical"),
                          inf = 7, pars = NULL,
                          epsilon = 1e-08,
                          maxit = 1000,
                          trace = FALSE,
                          gradient = FALSE) {
  stopifnot(is.character(criterion))
  criterion <- criterion[1]
  if(!(criterion %in% c("ML", "MLJ", "MLPJ", "MLJIC"))) {
    warning(paste0("currently recognized fitting criteria include:",
                   " 'ML' (Maximum Likelihood)",
                   ", 'MLJ' (ML with Jeffreys invariant prior) ",
                   ", 'MLPJ' (ML with partial Jeffreys prior (only slopes)) ",
                   ", and 'MLJIC' (MLJ with intercept correction)."))
    stop(paste0("fitting criterion ", criterion, " not recognized."))
  }

  if (!is.null(slope.signs)) {
    if (is.character(slope.signs)) {
      if (!all(slope.signs %in% c('negative', 'zero', 'positive')))
        stop("character valued 'slope.signs' must one of 'zero', 'negative' and 'positive'")
      slope.signs[slope.signs == 'zero'] <- '0'
      slope.signs[slope.signs == 'negative'] <- '-1'
      slope.signs[slope.signs == 'positive'] <- '1'
      slope.signs <- as.numeric(slope.signs)
    }
    else if (!all(slope.signs %in% c(-1, 0, 1)))
      stop("numeric valued 'slope.signs' must one of '0', '-1' and '1'")
  }
  else {
    slope.signs <- 0
  }

  stopifnot(is.character(method))
  method <- method[1]
  if (!all(method %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "nlminb"))) {
    stop("value of 'method' not recognized (see '?msbm.control' for a list of available methods)")
  }

  start <- match.arg(start)
  star.points <- match.arg(star.points)
  stopifnot(is.numeric(inf), is.finite(inf))
  stopifnot(is.null(pars) | is.numeric(pars))

  if (any(slope.signs != 0) & identical(method, "CG")) {
    if (trace) {
      warning("method 'CG' replaced by 'L-BFGS-B' to handle box constraints")
    }
    method <- "L-BFGS-B"
  }

  if (!all(c(is.numeric(epsilon), epsilon > 0)))
    stop("value of 'epsilon' must be > 0")
  epsilon <- epsilon[1]

  if (!all(c(is.numeric(maxit), maxit > 0)))
    stop("maximum number of iterations must be > 0")
  maxit <- ceiling(maxit)[1]

  list(criterion = criterion,
       slope.signs = slope.signs,
       unit.lambda = unit.lambda,
       method = method,
       start = start,
       star.points = star.points,
       inf = inf,
       pars = pars,
       epsilon = epsilon,
       maxit = maxit,
       trace = trace,
       gradient = gradient)
}
