
#' @export logit
#' @importFrom stats binomial
#' @importFrom stats plogis
#' @importFrom stats qlogis
#' @importFrom stats dlogis
# @rdname probit
logit <- function (me.family = "conjugate") {
  familytemp <- substitute(me.family)
  if (!is.character(familytemp))
    familytemp <- deparse(familytemp)
  if(identical(me.family, "canonical"))
    me.family <- "conjugate"
  okFams <- c("conjugate", "bridge", "norm", "logis")
  link <- "logit"
  if (familytemp %in% okFams[-1]) {
    familystats <- make.me.family (familytemp)
  }
  else if (identical(me.family, "conjugate")) {
    familytemp <- "bridge"
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
    qlogis(mu, location = 0, scale = 1,
           lower.tail = lower.tail, log.p = log.p)
  }

  linkinv <- function (eta, lower.tail = TRUE, log.p = FALSE) {
    plogis(eta, location = 0, scale = 1,
           lower.tail = lower.tail, log.p = log.p)
  }
  environment(linkfun) <- environment(linkinv) <- asNamespace("stats")

  if (familytemp %in% c("conjugate", "bridge")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
      plogis(eta / sqrt(1 + 3 * (sd^2)/(pi^2)), location = 0, scale = 1,
             lower.tail = lower.tail, log.p = log.p)
    }
    environment(melinkinv) <- asNamespace("stats")
  }
  else if (identical(familytemp, "norm")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
      expit.me (eta, sd = sd, lower.tail = lower.tail, log.p = log.p,
                family = "norm")
    }
    environment(melinkinv) <- asNamespace("msbreg")
  }
  else if (identical(familytemp, "logis")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
      expit.me (eta, sd = sd, lower.tail = lower.tail, log.p = log.p,
                family = "logis")
    }
    environment(melinkinv) <- asNamespace("msbreg")
  }
  else {
    stop("not yet implemented for families other than 'norm', 'logis', and 'bridge'")
  }

  mu.eta <- stats::binomial(link = 'logit')$mu.eta
  mu.eta.eta <- function (eta) {
    mu.etavalue <- stats::dlogis(eta)
    muvalue <- stats::plogis(eta)
    mu.etavalue * (1 - 2 * muvalue)
  }

  structure(list(family = "binomial",
                 link = link,
                 linkfun = linkfun,
                 linkinv = linkinv,
                 mefamilyname = familystats$name,
                 mefamily = familystats,
                 sigma = pi/sqrt(3),
                 melinkinv = melinkinv,
                 mu.eta = mu.eta,
                 mu.eta.eta = mu.eta.eta),
            class = "link")
}

expit.me <- function (eta, sd = 1, lower.tail = TRUE,
                      log.p = FALSE, family = "norm") {
  # Use the dimension of eta for the result
  fres <- eta

  # Change to the dimension of sd if required
  if (is.matrix(sd) & !is.matrix(eta))
    fres <- sd

  # Re-cycle vectors if required
  n <- max(length(eta), length(sd))
  eta <- rep(c(eta), length.out = n)
  sd <- rep(c(sd), length.out = n)

  # truncate eta (as in standard binomial GLM) and use plogis.menorm
  res <- apply (cbind(eta, sd), MARGIN = 1,
                FUN = switch(family,
                             norm = plogis.menorm,
                             logis = plogis.melogis))

  if (!lower.tail)
    res <- 1 - res

  if (log.p)
    res <- log(res)

  fres[] <- res

  return(fres)
}

plogis.menorm <- function(eta, abs.tol = 1e-8) {
  sd <- eta[2]
  eta <- eta[1]

  if (sd <= 0) {
    return(plogis(eta))
  }

  if (eta == 0) {
    return(0.5)
  }

  upper <- FALSE
  if (eta > 0) {
    eta <- -eta
    upper <- TRUE
  }

  p <- catch.conditions({
    integrate (f = function(x) {
      stats::plogis (x, location = 0, scale = 1,
             lower.tail = TRUE, log.p = FALSE) *
        stats::dnorm(x, mean = eta, sd = sd)
    },
    lower = -Inf, upper = Inf, abs.tol = abs.tol)$value
  })$value

  if (any(class(p) %in% c("simpleError", "error", "condition"))) {
    p <- catch.conditions({
      integrate (f = function(u) {
        stats::plogis (sd * u, location = 0, scale = 1,
               lower.tail = TRUE, log.p = FALSE) *
          stats::dnorm(u, mean = eta/sd, sd = 1)
      },
      lower = -Inf, upper = Inf, abs.tol = abs.tol)$value
    })$value

    if (any(class(p) %in% c("simpleError", "error", "condition")) &
        (abs.tol < .Machine$double.eps^0.25)) {
      p <- catch.conditions({
        integrate (f = function(u) {
          stats::plogis (sd * u, location = 0, scale = 1,
                  lower.tail = TRUE, log.p = FALSE) *
            stats::dnorm(u, mean = eta/sd, sd = 1)
        },
        lower = -Inf, upper = Inf, abs.tol = .Machine$double.eps^0.25)$value
      })$value
    }

    if (any(class(p) %in% c("simpleError", "error", "condition")))
      return(NA)
  }

  if (upper)
    p <- 1 - p

  return(p)
}

plogis.melogis <- function(eta, abs.tol = 1e-8) {
  sd <- eta[2]
  eta <- eta[1]

  if (sd <= 0) {
    return(plogis(eta))
  }

  if (eta == 0) {
    return(0.5)
  }

  upper <- FALSE
  if (eta > 0) {
    eta <- -eta
    upper <- TRUE
  }

  invsdlogis <- sqrt(3)/pi
  p <- catch.conditions({
    integrate (f = function(x) {
      stats::plogis (x, location = 0, scale = 1,
              lower.tail = TRUE, log.p = FALSE) *
        stats::dlogis(x, location = eta, scale = sd * invsdlogis)
    },
    lower = -Inf, upper = Inf, abs.tol = abs.tol)$value
  })$value

  if (any(class(p) %in% c("simpleError", "error", "condition"))) {
    p <- catch.conditions({
      integrate (f = function(u) {
        stats::plogis (sd * u, location = 0, scale = 1,
                lower.tail = TRUE, log.p = FALSE) *
          stats::dlogis(u, location = eta/sd, scale = invsdlogis)
      },
      lower = -Inf, upper = Inf, abs.tol = abs.tol)$value
    })$value
  }

  if (any(class(p) %in% c("simpleError", "error", "condition")) &
      (abs.tol < .Machine$double.eps^0.25)) {
    p <- catch.conditions({
      integrate (f = function(u) {
        stats::plogis (sd * u, location = 0, scale = 1,
                lower.tail = TRUE, log.p = FALSE) *
          stats::dlogis(u, location = eta/sd, scale = invsdlogis)
      },
      lower = -Inf, upper = Inf, abs.tol = .Machine$double.eps^0.25)$value
    })$value
  }

  if (any(class(p) %in% c("simpleError", "error", "condition")))
    return(NA)

  if (upper)
    p <- 1 - p

  return(p)
}

