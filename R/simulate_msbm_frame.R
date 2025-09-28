#' Simulate Multistage Binomial Responses
#'
#' Simulate one or more vectors of responses from
#' the distribution corresponding to the response in a
#' Multistage Binomial (MSB) model represented by a \code{"msbm"}
#' or a \code{"msbm.frame"} class object and a vector
#' of model parameters.
#'
#' @param object an object of class \code{"msbm"} or \code{"msbm.frame"},
#' typically returned by the functions \link[msbreg]{msbreg} or
#' \link[msbreg]{msbm.frame}.
#'
#' @param nsim,seed,... same as for the default method (\link[stats]{simulate}).
#'
#' @param theta numeric vector of model parameters. Note that \code{theta}
#' is required for an \code{object} of class \code{"msbm.frame"}, and must
#' be of the same length as the \code{$parnames} component of \code{object}.
#' For an \code{object} of class \code{"msbm"}, \code{theta} is optional,
#' with the default given by the \code{$coefficients} component of
#' \code{object}.
#'
#' @param link a character giving the name of the link function to be used.
#' Alternatively, link can be an object of class \code{"link"},
#' or a function such that \code{link()} returns an object of class \code{"link"}.
#' Defaults to the \link[msbreg]{logit} \code{link} for an \code{object} of
#' class \code{"msbm.frame"}. For an \code{object} of class \code{"msbm"},
#' the default is the \code{$link} component of \code{object}.
#'
#' @aliases simulate.msbm.frame
#' @aliases simulate.msbm
#'
#' @details
#' These are methods of the generic function \link[stats]{simulate} for
#' classes \code{"msbm"} and \code{"msbm.frame"}.
#' The functions provide alternatives to the \link[msbreg]{sim.msbdata} when
#' a \code{"msbm"} or \code{"msbm.frame"} object is already available.
#'
#' The method for class \code{"msbm"} simply extracts the model frame
#' (of class \code{"msbm.frame"}), and if not specified, the vector of
#' estimated model parameters and the \code{link} of the fit in
#' \code{object}, and passes them to the method for \code{"msbm.frame"}.
#'
#' @return a numeric vector when \code{nsim = 1},
#' a numeric \code{"matrix"} class object (\code{nsim} columns)
#' when \code{nsim > 1}. In either case, the returned value has
#' the two attributes:
# \describe{
#' \item{\code{"mu"}}{ the vector of success probabiliies used to generate
#' individual responses in the returned object;}
#' \item{\code{"seed"}}{ if argument \code{seed} is \code{NULL},
#' the attribute is the value of \link{.Random.seed} before the
#' simulation was started; otherwise it is the value of the argument with a
#' \code{"kind"} attribute with value \code{as.list(RNGkind())}.}
# }
#'
#' @seealso \link[stats]{simulate} for the generic and its default.
#' See also \link[msbreg]{sim.msbdata} to simulate a MSB model when
#' \code{object} is not yet available.
#'
#' @examples
#' ##** Simulated example: two-stage logistic model
#' data(test1data)
#' attr(test1data$y, "formula")
#'
#' ##Two-stage logistic model fit with unknown maximum success probability
#' # A multiplicative intercept included
#' MSBres <- msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                    lambda.formula = ~ 1,
#'                    weights = Total,
#'                    data = test1data)
#'
#' summary (MSBres)
#'
#' ## Simulate 1000 replicates of y from the fitted model
#' y <- simulate (MSBres, nsim = 1000)
#'
#' ## Compare the original response with a summary from the 1000 simulations
#' pick <- seq(1, by = 100, to = 1000)
#' cbind(y = test1data$y[pick],
#'       fitted = round(MSBres$fitted.values[pick] * test1data$Total[pick], 3),
#'       sim.mean = rowMeans(y[pick,]),
#'       sim.min = rowStats(y[pick,], stat = min),
#'       sim.median = rowStats(y[pick,], stat = median),
#'       sim.max = rowStats(y[pick,], stat = max))
#'
#' @exportS3Method stats::simulate
simulate.msbm.frame <- function (object, nsim = 1, seed = NULL,
                                 theta, link = "logit", ...) {
  stopifnot(nsim >= 1)

  #* Get Binomial link function
  eval(get.linkfun())

  #* Compute basic quantities and the vector 'mu' of success probabilities
  frame <- object
  eval(basic.mutate())

  size <- object$weights
  if (length(size) > 1) {
    size <- size[validmu]
  }

  #* Save/Set the random generator seed
  eval(setsave.RNGstate())

  ynames <- names(object$y)
  if (is.null(ynames)) {
    ynames <- rownames(object$input.matrix)
  }

  if (nsim == 1) {

    y <- numeric(length = object$dims$nobs) + NA

    y[validmu] <- stats::rbinom (sum(validmu), size = size,
                                 prob = mu[validmu])
    names(y) <- ynames
  }
  else {

    y <- matrix(NA, nrow = object$dims$nobs, ncol = nsim)

    y[validmu,] <- sapply(seq_len(nsim),
                          FUN = function (i) {
                            stats::rbinom (sum(validmu), size = size,
                                           prob = mu[validmu])
                          })

    rownames(y) <- ynames
    colnames(y) <- paste0("sim_", seq_len(nsim))
  }

  return (structure(y,
                    mu = mu,
                    seed = RNGstate))
}

#setMethod ("simulate",
#           signature = "msbm.frame",
#           definition = simulate.msbm.frame)

#' @rdname simulate.msbm.frame
#' @exportS3Method stats::simulate
simulate.msbm <- function (object, nsim = 1, seed = NULL,
                           theta = object$coefficients, link = object$link, ...) {
  stopifnot(nsim >= 1)

  return(simulate.msbm.frame(object = model.frame(object),
                             nsim = nsim, seed = seed,
                             theta = theta,
                             link = link, ...))
}

#setMethod ("simulate",
#           signature = "msbm",
#           definition = simulate.msbm)
