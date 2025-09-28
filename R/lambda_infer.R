#' Maximum Success Probability
#'
#' Empirically estimate the maximum success probability (\eqn{\lambda}) in
#' a Multistage Binomial (MSB) model.
#'
#' @param x real vector or matrix/data.frame of numeric predictor(s).
#' Each predictor represents a MSB model stage.
#'
#' @param y binary or binomial vector of response values (integers), its length must match
#' the number of rows of \code{x}.
#'
#' @param weights an optional vector of *weights* giving the number of trials
#' that resulted in the corresponding number of successes given in \code{y}.
#' Should be \code{NULL} (equivalent to \code{weights = 1}) or a numeric scalar
#' or a numeric vector (same length as \code{y}).
#'
#' @param scale.unit logical, if \code{TRUE} (the default), data columns in \code{x}
#' are scaled to unit variance.
#'
#' @param min.quantile,max.quantile scalars in the open \code{(0, 1)},
#' respectively the minimum and the maximum proportion of the data in \code{x}
#' to be used in the estimation procedure. The minimum \code{min.quantile} must
#' be greater than the positive tolerance value \code{eps}, but less
#' than 0.25, and the maximum \code{max.quantile} must be less than 0.5, but
#' greater than \code{min.quantile}.
#'
#' @param bin.size scalar in the open \code{(0, 1)}, size of a bin (proportion of
#' data points) in the sequence of bins used in the estimation.
#'
#' @param leftmost logical, should we always test against observations in only the leftmost bin?
#'
#' @param alpha scalar in the open \code{(0, 1)}, type I error rate for detecting
#' a significant shift in the success probability of the binary response \code{y}.
#'
#' @param alternative a character string specifying the alternative hypothesis
#' for each test. One of \code{'two.sided'} (default), \code{'greater'} or
#' \code{'less'}. Can be abbreviated (you can specify just the initial letter).
#' Passed to \link{fisher.test} and \link{binom.test}.
#'
#' @param conf.level confidence level of the returned confidence interval for
#' \eqn{\lambda}. Must be a single number between \code{0} and \code{1}.
#' Passed to \link{binom.test}.
#'
#' @param eps positive tolerance \eqn{\epsilon} used to decide when a (non random)
#' quantity is positive.
#'
#' @export infer.lambda
#' @import Rdpack
#'
#' @details
#' This routine uses a sequence of proportion tests (see \link{fisher.test}) to
#' detect the left asymptotic limit of the success probability \eqn{P(Y = y|x)}
#' of a binary response \eqn{Y}, as a function of a vector of predictors \eqn{x}.
#'
#' For one predictor (scalar \eqn{x}), the basic assumption of the estimation
#' procedure is that the observed range \eqn{[x_{min}, x_{max}]} of \eqn{x}
#' is (almost) the true meaningful range of \eqn{x}. In particular,
#' \code{infer.lambda} assumes that the observed \eqn{x} values include
#' (almost) the leftmost possible value of \code{x}. This implies that
#' \eqn{\beta_0 > x_{min}} under the MSB model
#' \eqn{P(Y = 1|x) = \lambda \times h(\beta_0 + \beta_1 x)} where
#' \eqn{\lambda} is the maximum success probability and \eqn{h} is a binomial
#' \code{link} function.
#'
#' If \eqn{P(Y = 1|x)} decreases with x values, then \eqn{\beta_1 \leq 0} and the
#' default returned \eqn{\lambda} value estimates the maximum success probability.
#' Otherwise (i.e. when  \eqn{\beta_1 > 0}), the returned \eqn{\lambda}
#' estimates the minimum success probability. If the maximum is desired in this
#' case, replace \code{x} by \code{-x} in the call to \code{infer.lambda}.
#'
#' For many predictors (when \code{x} is a matrix), each predictor is optionally
#' divided by its standard deviation (if \code{scale.unit = TRUE}, the default) to avoid
#' adding quantities in different units. Then the point of minima (i.e. the first
#' coordinate is the minimum of the first column in \code{x}, ..., the last
#' coordinate is the minimum of the last column in \code{x}) is then taken as
#' origin and the Euclidean distance of each row of \code{x} to this origin is
#' computed. The distance vector \code{d} is finally used in place of \code{x}
#' as in the case of one predictor.
#'
#' Note that the function simply estimates the asymptote \eqn{\lambda} as the
#' proportions of ones in the vicinity of the point of minima. It is the user's
#' responsibility to ensure that the asymptotic behavior happens near this point
#' of minima and not elsewhere. Clearly, when the *upper* asymptote is desired,
#' each column of \code{x} should be such that the proportion of ones in \code{y}
#' decreases with the column values. When the *lower* asymptote is desired,
#' simply replace the appropriate input for the *upper* asymptote (\code{x})
#' by \code{-x}.
#'
#' Also note that the estimation procedure directly implies that \eqn{\hat{\lambda}}
#' will be generally under (upper asymptote) or over (lower asymptote) estimated.
#'
# One could replace the sample estimate \eqn{\hat{\lambda}} of an upper asymptote
# by a corrected estimate such as: \eqn{\frac{4 \times \hat{\lambda} + UB_{\lambda}}{5}}
# where \eqn{UB_{\lambda}} is the upper \eqn{95\%} confidence bound of
# \eqn{\hat{\lambda}} (simple speculation).
#'
#' Confidence interval for the estimated \eqn{\lambda} is obtained using \link{binom.test}.
#'
#' The default argument values \code{min.quantile = bin.size},
#' \code{max.quantile = 3 * bin.size}, and \code{bin.size = 0.1}
#' are meant for small samples and small true \eqn{\lambda} value,
#' typically \eqn{n \leq 100} and \eqn{\lambda \leq 0.5}. The bin size should get
#' smaller as the sample size increases, the same holds for \code{max.quantile}
#' and \code{min.quantile}. Otherwise, estimates of the upper asymptote
#' \eqn{\hat{\lambda}} would be underestimated, and by large.
#'
#' @return A list with elements:
#' \item{lambda}{ point estimate of the asymptotic success probability.}
#' \item{conf.int}{ a vector of two values: lower and upper bounds of the estimate \code{lambda}.
#' The vector has attribute \code{conf.level} which indicates the confidence level of \code{conf.int}.}
#' \item{elbow.x}{ the \code{x} value (vector if \code{x} is a matrix) at which the success probability reaches \code{lambda}.}
#' \item{tests}{ a \code{data.frame} with rows corresponding to quantiles that delimit
#' the bins used for estimation, and columns:
#'    * \code{distance}: distance value (to the point of minima, see details) closing a bin;
#'    * \code{N}: total number of trial in a bin;
#'    * \code{n}: total number of successes among the \code{N} trial in a bin;
#'    * \code{lambda}: estimate from a bin and all left-side data points;
#'    * \code{odds.ratio}: odd ratio for change in the success
#'    probability when moving to a bin;
#'    * \code{p.value}: p value of the test
#'    of unit odd ratio against the alternative specified by \code{alternative}.}
#'
#' \item{call}{ the matched call to the function.}
#'
# @seealso
# \link{lambda.loess} for a semi-parametric estimation of
# \eqn{\lambda} when \eqn{x} is multivariate (up to 4-variate).
#'
#' @examples
#' ## Load 'msbreg'
#' library(msbreg)
#'
#' ## Simulate some data where 'x1' and 'x2' represent some ddg for
#' # a mutant protein folding/binding and 'y' is mutation viability
#' set.seed(155)
#' dset <- data.frame (x1 = runif(1000, min = -1, max = 4),
#'                     x2 = runif(1000, min = -1, max = 4))
#'
#' xframe <- msbm.frame( ~ x1 | x2,
#'                      lambda.formula = ~ 1,
#'                      data = dset)
#'
#' # Set model parameters: true lambda = 0.5
#' theta <- c(3, -1.5, 3, -1.5, qlogis(0.5))
#'
#' dset$y <- as.numeric(simulate(xframe, nsim = 1,
#'                               link = 'logit',
#'                               theta = theta))
#'
#' InfLambda <- infer.lambda(y = dset$y, x = cbind(dset$x1, dset$x2))
#'
#' InfLambda$lambda     # Point estimate
#' InfLambda$conf.int   # Confidence interval
#' InfLambda$tests        # Sequential tests performed
#'
# Can use univariate logistic regressions to find out the sign of beta_1
# correct so that y can be binomial, not just binary
infer.lambda <- function (x,
                          y, weights = NULL,
                          scale.unit = TRUE,
                          min.quantile = bin.size,     # must be > eps and < 0.25
                          max.quantile = 3 * bin.size, # must be <= 0.5
                          bin.size = 0.1,
                          leftmost = TRUE,    # should we always test against the leftmost bin?
                          alpha = 0.2,         # Type I error
                          alternative = c("two.sided", "less", "greater"),
                          conf.level = 0.95,
                          eps = 1e-5) {
  ## Call
  mcall <- match.call()

  ## Check arguments
  if(is.logical(y))
    y <- y+0
  if (is.null(weights))
    weights <- 1
  stopifnot(is.numeric(y), all(y >= 0), all(y == floor(y)))
  stopifnot(is.numeric(weights), all(weights >= 0), all(weights == floor(weights)))
  nobs <- length(y)
  if (is.data.frame(x))
    x <- as.matrix(x)
  stopifnot(is.numeric(x),  NROW(x) == nobs)
  stopifnot(any(y <= weights))
  weights <- rep(weights, length.out = nobs)
  stopifnot(eps > 0)
  stopifnot(min.quantile >= eps, min.quantile < 0.25)
  stopifnot(min.quantile < max.quantile, max.quantile <= 0.5)
  stopifnot(bin.size > 0, bin.size < max.quantile)
  stopifnot((bin.size - eps) <= (max.quantile - min.quantile)) # to have minimum 2 bins

  ## Derive a distance measure from the matrix x
  x0 <- x
  if (NCOL(x) > 1) {
    # Normalize
    if (scale.unit) {
      x <- scale(x, center = FALSE, scale = TRUE)    # divide by standard deviation
    }

    # Distance to the point of minima (min_x1, ..., min_xp)
    x <- t(t(x) - apply(x, MARGIN = 2, FUN = min)) # Each variable starts at zero (min)
    d <- sqrt(rowSums(x^2))
  }
  else {
    if (scale.unit) {
      x <- scale(x, center = FALSE, scale = TRUE)    # divide by standard deviation
    }
    d <- x - min(x)
  }

  ## Define bins of distance values (quantiles)
  bins <- seq(from = min.quantile, to = max.quantile, by = bin.size)
  nbins <- max(length(bins), 2)
  bins[1] <- min.quantile
  bins[nbins] <- max.quantile
  qmin <- quantile(d, prob = bins)

  ## Find number of trials and successes per bin
  n <- sapply(qmin, FUN = function(dmax) {
    pick <- d <= dmax
    c(sum(weights[pick]), sum(y[pick], na.rm = TRUE))
  })
  N <- c(n[1,1], diff(n[1,]))
  n <- c(n[2,1], diff(n[2,]))

  # Run sequential proportion tests (Fisher's Exact Test)
  res.df <- matrix(NA, nrow = nbins, ncol = 5)
  colnames(res.df) <- c('N', 'n', 'lambda', 'odds.ratio', 'p.value')
  rownames(res.df) <- names(N)
  x_argo <- x_arg <- n[1]
  n_argo <- n_arg <- N[1]
  res.df[1,] <- c(N = n_argo, n = x_argo, lambda = x_argo/n_argo, odds.ratio = NA, p.value = NA)
  for (k in 2:nbins) {
    # Bug: was using all observations at a bin and the left most
    #x_arg <- x_arg + n[k]
    #n_arg <- n_arg + N[k]
    #matrix(c(x_argo,          x_arg,
    #         n_argo - x_argo, n_arg - x_arg),
    #       nrow = 2, byrow = TRUE)

    # Corrected to only use the current bin
    Ftests <- fisher.test (matrix(c(x_argo,          n[k],
                                    n_argo - x_argo, N[k] - n[k]),
                                  nrow = 2, byrow = TRUE),
                           alternative = alternative[1])

    if (!leftmost) {
      x_argo <- x_argo + n[k]
      n_argo <- n_argo + N[k]
    }

    res.df[k,] <- c(N[k], n[k],
                    sum(n[1:k])/sum(N[1:k]),
                     Ftests$estimate,
                     Ftests$p.value)
  }
  res.df <- as.data.frame(res.df)

  # Find any significant change point
  pick <- which(res.df$p.value[-1] <= alpha)
  if (length(pick)) {
    pick <- pick[1] # Could think this 'max(1, pick[1] - 1)' is the right code but we removed the first element of 'p.value' when defining 'pick'
  }
  else {
    pick <- nbins
  }

  # Get lambda value
  nfinal <- sum(n[1:pick])
  Nfinal <- sum(N[1:pick])
  lambda <- nfinal / Nfinal # Use all observations on the left side and at the picked bin

  # Set the elbow value
  elbow.d <- qmin[pick]
  whichd <- which(d <= elbow.d)
  if (length(whichd) == 1) {
    take <- whichd
  }
  else {
    take <- which.max(d[whichd])
    take <- whichd[take[1]]
  }

  if (NCOL(x) > 1) {
    elbow.x <- x0[take,]
  }
  else {
    elbow.x <- x0[take]
  }
  attr(elbow.x, 'position') <- take

  # Confidence interval
  conf.int <- binom.test (x = c(nfinal, Nfinal - nfinal),
                       alternative = alternative[1],
                       conf.level = conf.level[1])$conf.int

  res.df <- cbind(distance = qmin, res.df)

  # Adjusting an estimate above the average of 'y'
  ave.y <- mean(y/weights)
  if (lambda <= ave.y) {
    lambda <- ave.y + eps
    conf.int[1] <- min(conf.int[1], ave.y + eps)
  }

  # Adjusting a zero estimate
  if (lambda < eps) {
    lambda <- max(eps, min(1-eps,  2 * sum(n)/sum(N)))
    # conf.int <- pmax(eps, pmin(1-eps,  conf.int))
  }

  return(list(lambda = lambda,
              elbow.x = elbow.x,
              which = whichd,
              conf.int = conf.int,
              tests = res.df,
              distance = d,
              call = mcall))

}

get.lambda <- function (x, y, x.quantile = 0.1) {
  qmin <- quantile(x, prob = x.quantile)
  L <- mean(y[x <= qmin], na.rm = TRUE)
  return(L)
}
