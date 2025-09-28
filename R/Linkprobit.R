# make.me.family; return class me-family

#' Link Objects for Model Fitting
#'
#' Link objects provide a convenient way to specify the details of the links
#' used by model fitting functions such as \link{msbreg}.
#'
#' @param me.family a specification for the measurement error model family.
#' This can be a name/expression, a literal character string, a length-one
#' character vector, or an object of class \code{"me-family"} (such as generated
#' by \link{make.me.family}) provided it is not specified \code{via} one of the
#' standard names such as \code{probit}, \code{sprobit}, \code{robit},
#' \code{logit}, \code{cloglog}, \code{cauchit}, and \code{tcauchit}.
#'
#' Currently, the only available character option is the default
#' \code{me.family = "conjugate"}. This default picks, if it exists, the natural
#' (or conjugate) measurement distribution family for a given link function
#' \code{link}: that is the distribution \code{distr} such that the mixture
#' of the inverse link function (which is a cumulative distribution function)
#' by \code{distr} is of the same family as this inverse link function.
#' This requires the existence of first and second order moments for the
#' distribution corresponding to the link function, a condition which is
#' satisfied for all the link functions given here, except for the
#' \code{cauchit} link (see **Details**).
#'
#' @param df degrees of freedom of the Student t distribution (\link[stats]{TDist}),
#' \eqn{df >2} required to ensure finite variance \eqn{df/(df - 2)}.
#' The infinite value \code{df = Inf} is allowed. The variance is then the limit 1.
#'
#' @param lambda real skewness parameter of the skew normal distribution
#' (\link[msbreg]{snorm}). Note that the default \code{lambda = 0.75} corresponds
#' to a right skewed link with skewness (third standardized moment)
#' approximately equal to \code{0.8173}.
#'
#' @param object the function \code{link} accesses the link of \code{object}
#' which are stored within objects created by modelling functions
#' (e.g., \code{msbreg}).
#'
#' @param ... Additional arguments passed to or from other methods.
#'
#'
#' @details
#' The links given here are for binomial model fittings.
#' If the model include measurement errors, only the conjugate measurement
#' distribution family, when it exists, is currently available for each link.
#' For the \code{probit} link, the conjugate family is the Gaussian family
#' (see \link[stats]{Normal}). This may be explicitly specified by
#' \code{probit(me.family = "norm")}.
#' The \code{sprobit} link which is based on the Sicard skew normal
#' distribution does not have a known conjugate family. Currently, the
#' conjugate family is given by the default \code{sprobit(me.family = "conjugate")}.
#' For the \code{logit} link, the conjugate family is the logistic Bridge
#' distribution (see \link[msbreg]{dbridge}), the explicit specification
#' is \code{logit(me.family = "bridge")}.
#'
#' The \code{robit} link is based on the Student t distribution and thus requires
#' the specification of a degrees of freedom parameter \code{df} which
#' determines both the shape (tail behavior) and the variance of the link
#' \insertCite{liu2004robit,kim2008flexible}{msbreg}.
#' The degrees of freedom parameter must satisfy \code{df} \eqn{>2} to ensure
#' a finite link variance \eqn{\sigma^2 = df/(df - 2)} (\eqn{1 < \sigma^2 < \infty}).
#' So the default \code{df = 2.1} corresponds to a heavy tail distribution with
#' a variance of \code{21}.
#' Note that the argument \code{df} is a fixed parameter for model fitting.
#' The conjugate measurement error family for the \code{robit} link is the
#' Student t family (this assumes that the link and the measurement errors
#' are uncorrelated but dependent, as always under the multivariate t
#' distribution).
#'
#' The \code{cloglog} link is a skewed link with the log-positive stable
#' distribution as its conjugate measurement error family
#' \insertCite{hougaard1986survival}{msbreg}.
#' The explicit specification here is \code{cloglog(me.family = "lstab")}.
#'
#' The \code{cauchit} link does not have a conjugate measurement error family
#' since it has no well defined moment.
#' The \code{tcauchit} link is based on the symmetric truncated Cauchy
#' distribution which retains \eqn{99.99\%} probability mass of the standard
#' Cauchy distribution (the approximate truncation region is
#' \eqn{[-3183, 3183]}).
# qcauchy(.9999)
#' This truncated Cauchy distribution has well defined moments.
#' The \code{tcauchit} link thus has a conjugate measurement error family,
#' with a distribution function corresponding to a mixture of Cauchy distributions.
#'
# In addition to its \code{"conjugate"} family (if any),
# each link accepts the family \code{"norm"} (Gaussian measurement
# errors) which is equivalent to \code{"conjugate"} for the \code{probit} link.
# The \code{logit} link also accepts the family \code{"logis"} (Logistic
# measurement errors) and the family \code{"bridge"} (Bridge measurement
# errors) which is its \code{"conjugate"} family.
# The \code{cloglog} link also accepts \code{"lstab"} which is the same as the
# \code{"conjugate"} family.
#'
#' The function \code{link} is generic.
#' The default method extracts the \code{$link} component of \code{object}.
#' If \code{object} has no \code{$link} component, and is of a class with no
#' \code{link} method, the default method will output \link{NULL}.
#'
#' @usage
#' probit(me.family = "conjugate")
#'
#' sprobit(me.family = "conjugate", lambda = 0.75)
#'
#' robit(me.family = "conjugate", df = 2.1)
# df = 2.87345 gives variance close to pi^2/3
#'
#' logit(me.family = "conjugate")
#'
#' cloglog(me.family = "conjugate")
#'
#' cauchit(me.family = "conjugate")
#'
#' tcauchit(me.family = "conjugate")
#'
#' @return An object of class \code{"link"} (which has a concise \code{print}
#' method). A proper \code{"link"} object is a list with elements:
#'
#' \item{family}{character: the family name of the response (i.e. \code{"binomial"}).}
#' \item{link}{character: the link name.}
#' \item{linkfun}{function: the link.}
#' \item{linkinv}{function: the inverse of the link function.}
#' \item{mefamilyname}{character: the measurement error distribution family name.}
#' \item{mefamily}{a measurement error family, object of class \code{"me-family"}
#' as returned by \link{make.me.family}.}
#' \item{sigma}{numeric: the standard deviation of the distribution
#' corresponding to the link function (\link{NaN} if this is undefined).
#' The returned \code{sigma} element may have some attributes, e.g.,
#' for the \code{robit} link, \code{sigma} has attribute \code{df},
#' i.e. the degree of freedom of the corresponding Student t distribution.}
#' \item{melinkinv}{function: the inverse link function accounting for
#' measurement errors.}
#' \item{mu.eta}{ function: first derivative of the inverse-link function with
#' respect to the linear predictor.}
#' \item{mu.eta.eta}{ function: second derivative of the inverse-link function with
#' respect to the linear predictor.}
#'
#' @note
#' The design was inspired by \code{R}'s \link{family} functions.
#'
#' @export probit
#' @export link
#'
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @importFrom stats dnorm
#'
#' @aliases link logit
#' @aliases cloglog robit
#' @aliases cauchit tcauchit
#' @aliases sprobit
#'
#' @references
#' \insertAllCited{}
#'
#' @import stats
#'
#' @seealso
#' \link{make.me.family}.
#'
probit <- function (me.family = "conjugate") {
  familytemp <- substitute(me.family)
  if (!is.character(familytemp))
    familytemp <- deparse(familytemp)
  if(identical(me.family, "canonical"))
    me.family <- "conjugate"
  okFams <- c("conjugate", "norm", "logis", "bridge")
  link <- "probit"
  if (familytemp %in% okFams[-1]) {
    familystats <- make.me.family (familytemp)
  }
  else if (identical(me.family, "conjugate")) {
    familytemp <- "norm"
    familystats <- make.me.family (familytemp)
  }
  else if (is.character(me.family)) {
    familystats <- make.me.family (me.family)
  }
  else {
    if (inherits(me.family, "me-family")) {
      familystats <- me.family
    }
    else {
      stop(gettextf("me.family \"%s\" not available for %s link; available me.families are %s",
                    me.family, link, paste(sQuote(okFams), collapse = ", ")),
           domain = NA)
    }
  }

  if (is.null(familystats$name))
    familystats$name <- familytemp

  linkfun <- function(mu, lower.tail = TRUE, log.p = FALSE) {
    qnorm(mu, lower.tail = lower.tail, log.p = log.p)
  }

  linkinv <- function(eta, lower.tail = TRUE, log.p = FALSE) {
    thresh <- -qnorm(.Machine$double.eps)
    eta <- pmin(pmax(eta, -thresh), thresh)
    pnorm (eta, lower.tail = lower.tail, log.p = log.p)
  }
  environment(linkfun) <- environment(linkinv) <- asNamespace("stats")

  if (familytemp %in% c("conjugate", "norm")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
      SD <- sqrt(1 + sd^2)
      thresh <- -qnorm(.Machine$double.eps) * SD
      eta <- pmin(pmax(eta, -thresh), thresh)
      pnorm (eta / SD, lower.tail = lower.tail, log.p = log.p)
    }
    environment(melinkinv) <- asNamespace("stats")
  }
  else if (identical(familytemp, "logis")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
      invbit.me (eta, sd = sd, lower.tail = lower.tail, log.p = log.p,
                 family = "logis")
    }
    environment(melinkinv) <- asNamespace("msbreg")
  }
  else if (identical(familytemp, "bridge")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
      invbit.me (eta, sd = sd, lower.tail = lower.tail, log.p = log.p,
                 family = "bridge")
    }
    environment(melinkinv) <- asNamespace("msbreg")
  }
  else {
    stop("not yet implemented for families other than 'norm', 'logis', and 'bridge'")
  }

  mu.eta <- binomial(link = 'probit')$mu.eta
  mu.eta.eta <- function (eta) {
    - eta * stats::dnorm(eta)
  }

  structure(list(family = "binomial",
                 link = link,
                 linkfun = linkfun,
                 linkinv = linkinv,
                 mefamilyname = familystats$name,
                 mefamily = familystats,
                 sigma = 1,
                 melinkinv = melinkinv,
                 mu.eta = mu.eta,
                 mu.eta.eta = mu.eta.eta),
            class = "link")
}

invbit.me <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE,
                             family = "logis") {
  # Use the dimension of eta for the result
  fres <- eta

  # Change to the dimension of sd if required
  if (is.matrix(sd) & !is.matrix(eta))
    fres <- sd

  # Re-cycle vectors if required
  n <- max(length(eta), length(sd))
  eta <- rep(c(eta), length.out = n)
  sd <- rep(c(sd), length.out = n)

  # truncate eta (as in standard binomial GLM) and use pnorm.melogis
  thresh <- -qnorm(.Machine$double.eps) * sqrt(1 + sd^2)
  eta <- pmin(pmax(eta, -thresh), thresh)
  res <- apply (cbind(eta, sd), MARGIN = 1,
                FUN =   switch(family,
                               logis = pnorm.melogis,
                               bridge = pnorm.mebridge))

  if (!lower.tail)
    res <- 1 - res

  if (log.p)
    res <- log(res)

  fres[] <- res

  return(fres)
}

pnorm.melogis <- function(eta, abs.tol = 1e-8) {
  sd <- eta[2]
  eta <- eta[1]

  if (sd <= 0) {
    return(pnorm(eta))
  }

  if (eta == 0) {
    return(0.5)
  }

  invsdlogit <- sqrt(3)/pi
  p <- catch.conditions({
    integrate (f = function(x) {
      pnorm (x, mean = 0, sd = 1,
             lower.tail = TRUE, log.p = FALSE) *
        dlogis(x, location = eta, scale = sd * invsdlogit)
    },
    lower = -Inf, upper = Inf, abs.tol = abs.tol)$value
  })$value

  if (any(class(p) %in% c("simpleError", "error", "condition"))) {
    p <- catch.conditions({
      integrate (f = function(u) {
        pnorm (sd * u, mean = 0, sd = 1,
               lower.tail = TRUE, log.p = FALSE) *
          dlogis(u, location = eta/sd, scale = invsdlogit)
      },
      lower = -Inf, upper = Inf, abs.tol = abs.tol)$value
    })$value

    if (any(class(p) %in% c("simpleError", "error", "condition")) &
        (abs.tol < .Machine$double.eps^0.25)) {
      p <- catch.conditions({
        integrate (f = function(u) {
          pnorm (sd * u, mean = 0, sd = 1,
                 lower.tail = TRUE, log.p = FALSE) *
            dlogis(u, location = eta/sd, scale = invsdlogit)
        },
        lower = -Inf, upper = Inf, abs.tol = .Machine$double.eps^0.25)$value
      })$value

      if (any(class(p) %in% c("simpleError", "error", "condition")))
        return(NA)
    }
  }

  return(p)
}

pnorm.mebridge <- function(eta, abs.tol = 1e-8) {
  sd <- eta[2]
  eta <- eta[1]

  if (sd <= 0) {
    return(pnorm(eta))
  }

  if (eta == 0) {
    return(0.5)
  }

  p <- catch.conditions({
    integrate (f = function(x) {
      pnorm (x, mean = 0, sd = 1,
             lower.tail = TRUE, log.p = FALSE) *
        dbridge (x, location = eta, scale = pi / sqrt(3 * sd^2 + pi^2))
    },
    lower = -Inf, upper = Inf, abs.tol = abs.tol)$value
  })$value

  if (any(class(p) %in% c("simpleError", "error", "condition"))) {
    p <- catch.conditions({
      integrate (f = function(u) {
        pnorm (sd * u, mean = 0, sd = 1,
               lower.tail = TRUE, log.p = FALSE) *
          dbridge (u, location = eta/sd, scale = pi / sqrt(3 + pi^2))
      },
      lower = -Inf, upper = Inf, abs.tol = abs.tol)$value
    })$value

    if (any(class(p) %in% c("simpleError", "error", "condition")) &
        (abs.tol < .Machine$double.eps^0.25)) {
      p <- catch.conditions({
        integrate (f = function(u) {
          pnorm (sd * u, mean = 0, sd = 1,
                 lower.tail = TRUE, log.p = FALSE) *
            dbridge (u, location = eta/sd, scale = pi / sqrt(3 + pi^2))
        },
        lower = -Inf, upper = Inf, abs.tol = .Machine$double.eps^0.25)$value
      })$value

      if (any(class(p) %in% c("simpleError", "error", "condition")))
        return(NA)
    }
  }

  return(p)
}

#' @rdname probit
link <- function(object, ...) {
  UseMethod("link")
}

link.default <- function(object, ...) {
  object$link
}

setGeneric(name = "link",
           def = link.default)
