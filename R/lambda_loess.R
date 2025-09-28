# Maximum success probability
#
# Semi-parametric estimation of the maximum success probability \eqn{\lambda}
# in a (multiple) multistage binomial (MSB) model. This provides a robust
# estimation of \eqn{\lambda}. The function also returns initial MSB
# regression parameters, but these estimates are not intended as final
# estimates. Currently limited to **at most four** predictors.
#
# @param formula a formula specifying the numeric response and one to **four** numeric predictors.
# @param link same as for \link{msbreg} or \link{glm}.
# @param data,weights,subset,na.action,start same as for \link{glm}, but
# here, the argument \code{data} is mandatory and must be a \code{data.frame}.
# @param Lstart optional starting value for the maximum success probability.
# @param span,enp.target,degree,parametric,drop.square,normalize,family,control.loess arguments passed to \link{loess}.
# @param conf.level confidence level for the returned confidence intervals.
# @param df optional degrees of freedom parameter for confidence interval building.
# @param eps positive small tolerance value to ensure that smoothed probabilities range in \code{(0, 1)}.
# @param control.nls an optional list of control settings. Passed to \link[stats]{nls}.
# @param lower,upper optional bounds on parameters, passed to \link[stats]{nls}.
# @param ... directly supplied \link{loess} control parameters (if \code{control.loess} is not specified).
#
# @details
# This routine is similar to \link{infer.lambda} in
# a multiple binomial regression setting with up to *four* predictors. However,
# \code{lambda.loess} focuses on and robustly estimates only the maximum of
# the success probability. The returned estimates of other parameters should be
# used only as starting values for a more conventional estimation process.
#
# The function simply calls \link{loess} to fit a local polynomial regression to
# the binomial response specified in \code{formula} (analogous to the
# sequential/moving \link{binom.test} used by \link{infer.lambda}).
# Then, the predicted individual mean responses are used as known success
# probabilities in a MRB model framework to estimate the model parameters.
#
# To do:
#
# * consider a more general (more covariates) smoothing function to handle
# more predictors;
#
# * use the loess.as function of \code{fANCOVA} to automate the selection of the
# \code{span} for the case of one or two predictors.
#
# @export lambda.loess
# @importFrom fANCOVA loess.as
#
# @return A list with the following elements:
#
# \item{L}{ a list with components \code{estimate} (point estimate of the maximum success probability),
# \code{se} (the standard error of the point estimate), \code{conf.int} (end points-transform confidence interval),
# and \code{conf.int.t} (Student t approximation based confidence interval).}
# \item{delta}{ a list with components \code{estimate} (point estimate of the link-scale transform of the maximum success probability),
# \code{se} (the standard error of the point estimate).}
# \item{L1, Lhat}{ \link[stats]{nls} fit objects with an additional component \code{summary} (result of applying summary on the \code{nls} object).
# \code{L1} is an initialization fit assuming \code{L = 1}, and \code{Lhat} is the final fit that estimated \code{L}.
# Both fits are on the \code{log} of the smoothed response against the predictors in \code{formula} under the MSB model.}
#
# @seealso
# See \link{infer.lambda} for an empirical estimation method.
#
# @examples
# # Simulate some data: 2 continuous predictors
# set.seed(155)
# dset <- data.frame (x1 = runif(1000, min = -1, max = 4),
#                     x2 = runif(1000, min = -1, max = 4))
#
# xframe <- msbm.frame( ~ x1 | x2,
#                       lambda.formula = ~ 1,
#                       data = dset)
#
# # Set model parameters: true lambda = 0.5
# theta <- c(3, -1.5, 3, -1.5, qlogis(0.5))
#
# dset$y <- as.numeric(simulate(xframe, nsim = 1,
#                               link = msbreg::logit(),
#                               theta = theta))
#
# # Mode posterior likelihood estimate (with Jeffrey's prior)
# msbreg (y ~ x1 | x2,
#         lambda.formula = ~ 1,
#         data = dset)$lambda.values
#
# # Semi-parametric estimate
# lambda.loess (y ~ x1 + x2,
#               data = dset,
#               lower = -3,
#               upper = c(rep(3, 3), -1))$L$estimate
#
# # Empirical estimate
# infer.lambda(y = dset$y, x = cbind(dset$x1, dset$x2))$lambda
#
lambda.loess <- function (formula, # All variables must be in 'data', formula elements are not read from the calling environment.
                          data, # must be a data.frame, and not an environment (unlike 'glm')
                          link = logit,
                          weights,
                          subset,
                          na.action = getOption("na.action"),
                          start = NULL, Lstart = NULL,
                          span = 0.75, enp.target, degree = 2,
                          parametric = FALSE, drop.square = FALSE, normalize = TRUE,
                          family = c("gaussian", "symmetric"),
                          conf.level = 0.95, df = NULL, eps = 1e-8,
                          control.nls = list(),
                          control.loess = loess.control(...),
                          lower = -4, upper = 4, ...) {
  ### Checking
  stopifnot(!missing(formula))
  stopifnot(!missing(data))
  data <- as.data.frame(data)
  if (is.null(Lstart))
    Lstart <- 0.5
  stopifnot(is.numeric(Lstart))
  stopifnot(eps > 0, eps < 0.10001)

  ### Matched call
  mcall <- match.call()

  ### Get Binomial link function
  if (is.character(link))
    link <- get(link, mode = "function", envir = parent.frame())
  if (is.function(link))
    link <- link()
  if (is.null(link$link)) {
    print(link)
    stop("'link' not recognized")
  }

  ### Standard Binomial fit
  if (missing(weights) & missing(subset))
    SMBfit <- glm(formula,
                  family = binomial(link = link$link),
                  data = data,
                  na.action = na.action,
                  start = start)
  else if (missing(weights))
    SMBfit <- glm(formula,
                  family = binomial(link = link$link),
                  data = data,
                  subset = subset,
                  na.action = na.action,
                  start = start)
  else if (missing(subset))
    SMBfit <- glm(formula,
                  family = binomial(link = link$link),
                  data = data,
                  weights = weights,
                  na.action = na.action,
                  start = start)
  else
    SMBfit <- glm(formula,
                  family = binomial(link = link$link),
                  data = data,
                  weights = weights,
                  subset = subset,
                  na.action = na.action,
                  start = start)
  #data <- SMBfit$data

  # Find the names and positions of predictors in formula
  xnames <- names(SMBfit$coefficients)
  xnames <- xnames[xnames %in% colnames(data)]
  #xpositions <- which(colnames(data) %in% xnames)

  ### LOESS (Locally Estimated Scatterplot Smoothing) is a nonparametric method
  # for smoothing a series of data in which no assumptions are made about the
  # underlying structure of the data. LOESS uses local regression to fit a smooth
  # curve through a scatterplot of data.
  loess_mod <- loess(formula,
                     data = data,
                     span = span, enp.target = enp.target, degree = degree,
                     parametric = parametric, drop.square = drop.square,
                     normalize = normalize,
                     family = family,
                     control.loess = control.loess)

  # Ensure the smoothed curve is bounded to (0, 1)
  data$yhat <- pmin(pmax(loess_mod$fitted, eps), 1-eps)
  data$logyhat <- log(data$yhat)

  ### Formula for 'nls' and starting value
  p <- length(xnames)
  # Initial parameter vector for 'nls'
  start <- SMBfit$coefficients

  intercept <- '(Intercept)' %in% names(SMBfit$coefficients)
  if (intercept) {
    rhs <- paste0(sapply(1:p,
                         FUN = function(j) {
                           paste0("link$linkinv(beta_0", xnames[j], " + beta_", xnames[j], " * ", xnames[j], ", log.p = TRUE)")
                         }),
                  collapse = " + ")

    start <- c(rbind(rep(start[1]/p, p),
                     start[-1]))
    names(start) <- c(rbind(paste0('beta_0', xnames),
                            paste0('beta_', xnames)))
  }
  else {
    rhs <- paste0(sapply(1:p,
                         FUN = function(j) {
                           paste0("link$linkinv(beta_", xnames[j], " * ", xnames[j], ", log.p = TRUE)")
                         }),
                  collapse = " + ")
    names(start) <- paste0('beta_', xnames)
  }
  start <- pmax(pmin(start, rep(upper, length.out = length(start))),
                rep(lower, length.out = length(start)))

  ### Assume L = 1
  nlformula <- reformulate(rhs, 'logyhat')

  NL1 <- nls (nlformula,
              start = start,
              data = data,
              lower = lower,
              upper = upper,
              algorithm = "port")
  SumNL1 = summary(NL1)

  ### Assume L in (0, 1)
  rhs <- paste0('(-exp(delta)) + ', rhs)
  nlformula <- reformulate(rhs, 'logyhat')
  start0 <- start
  start <- c(SumNL1$coefficients[,1], log(-log(Lstart[1])))
  names(start) <- c(names(start0), 'delta')
  start <- pmax(pmin(start, rep(upper, length.out = length(start))),
                rep(lower, length.out = length(start)))
  NLhat <- nls (nlformula,
                start = start,
                data = data,
                lower = lower,
                upper = upper,
                algorithm = "port")
  SumNLhat = summary(NLhat)

  # Save summaries of nls fit objects
  NL1$summary <- SumNL1
  NLhat$summary <- SumNLhat

  # L estimate
  delta <- SumNLhat$coefficients[length(start),1]
  se.delta <- SumNLhat$coefficients[length(start),2]
  Lvalue <- exp(-exp(delta))
  se.L <- se.delta * Lvalue * exp(delta) # 'Delta' Method to find se.L

  # residual degrees of freedom
  if (is.null(df))
    df <- NLhat$summary$df[2]
  stopifnot(is.numeric(Lstart))

  #End point-transform confidence interval
  tvalue <- qt((1+conf.level)/2, df = df)
  conf.int.L <- c(delta + se.delta * tvalue, delta - se.delta * tvalue)
  conf.int.L <- exp(-exp(conf.int.L))
  attr(conf.int.L, 'conf.level') <- conf.level

  # Student t approximation based confidence interval
  t.conf.int.L <- c(max(0, Lvalue - se.L * tvalue), min(1, Lvalue + se.L * tvalue))
  attr(t.conf.int.L, 'conf.level') <- conf.level
  out <- list(L = list(estimate = Lvalue,
                       se = se.L,
                       conf.int = conf.int.L,
                       conf.int.t = t.conf.int.L),
              delta = list(estimate = delta,
                           se = se.delta),
              L1 = NL1,
              Lhat = NLhat)

  return(structure(out, class = 'loess.lambda'))

}
