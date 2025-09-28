
#' @export robit
#' @importFrom stats binomial
#' @importFrom stats pt
#' @importFrom stats qt
#' @importFrom stats dt
# @rdname probit
# Add a scale parameter (generalized t)
# Consider the prostate cancer data in kim2008flexible for MSB
# Data from d2002biochemical
robit <- function (me.family = "conjugate", df = 2.1) {
  stopifnot(df > 2)
  familytemp <- substitute(me.family)
  if (!is.character(familytemp))
    familytemp <- deparse(familytemp)
  if(identical(me.family, "canonical"))
    me.family <- "conjugate"
  okFams <- c("conjugate", "t", "none")
  link <- "robit"
  if (familytemp %in% okFams[-1]) {
    familystats <- make.me.family (familytemp)
  }
  else if (identical(me.family, "conjugate")) {
    familytemp <- "t"
    familystats <- make.me.family ("t")
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
    qt(mu, df = df, lower.tail = lower.tail, log.p = log.p)
  }

  linkinv <- function (eta, lower.tail = TRUE, log.p = FALSE) {
    thresh <- -qt(.Machine$double.eps, df = df)
    eta <- pmin(pmax(eta, -thresh), thresh)
    pt(eta, df = df, lower.tail = lower.tail, log.p = log.p)
  }

  if (is.finite(df))
    sigma2 <- df / (df - 2)
  else
    sigma2 <- 1
  if (familytemp %in% c("conjugate", "none", "t")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
      thresh <- -qt(.Machine$double.eps, df = df)

      eta <- eta / sqrt(1 + (sd^2) / sigma2)
      eta <- pmin(pmax(eta, -thresh), thresh)
      pt(eta, df = df, lower.tail = lower.tail, log.p = log.p)
    }

    familystats$name <- 'conjugate'
  }
  else {
    stop("not yet implemented for families other than 'conjugate'")
  }

  mu.eta <- function (eta) {
    dt(eta, df = df)
  }

  mu.eta.eta <- function (eta) {
    - eta * ((df+1) / df) * dt(eta, df = df) / (1 + (eta^2)/df)
  }

  structure(list(family = "binomial",
                 link = link,
                 linkfun = linkfun,
                 linkinv = linkinv,
                 mefamilyname = familystats$name,
                 mefamily = familystats,
                 sigma = structure(sqrt(sigma2), df = df),
                 melinkinv = melinkinv,
                 mu.eta = mu.eta,
                 mu.eta.eta = mu.eta.eta),
            class = "link")
}

