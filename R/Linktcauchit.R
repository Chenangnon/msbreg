
# Use a truncated cauchy distribution within +- 3183.099
#' @export tcauchit
#' @importFrom stats binomial
#' @importFrom stats pcauchy
#' @importFrom stats qcauchy
#' @importFrom stats dcauchy
# @rdname probit
tcauchit <- function (me.family = "conjugate") {
  thresh <- stats::qcauchy(.9999)
  familytemp <- substitute(me.family)
  if (!is.character(familytemp))
    familytemp <- deparse(familytemp)
  if(identical(me.family, "canonical"))
    me.family <- "conjugate"
  okFams <- c("conjugate", "none")
  link <- "tcauchit"
  if (familytemp %in% okFams[-1]) {
    familystats <- make.me.family ('none')
  }
  else if (identical(me.family, "conjugate")) {
    familytemp <- "conjugate"
    familystats <- make.me.family ('none') ######################### CORRECT !!!
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

  fname <- familystats$name
  if (is.null(fname))
    fname <- familystats$name <- familytemp

  linkfun <- function(mu, lower.tail = TRUE, log.p = FALSE) {
    mu <- pmax(pmin(mu, 1 - .Machine$double.eps),
               .Machine$double.eps)
    qtcauchy(mu, location = 0, scale = 1, tau = thresh,
             lower.tail = lower.tail, log.p = log.p)
  }

  linkinv <- function (eta, lower.tail = TRUE, log.p = FALSE) {
    eta <- pmin(pmax(eta, -thresh + .Machine$double.eps),
                thresh - .Machine$double.eps)
    ptcauchy(eta, location = 0, scale = 1, tau = thresh,
             lower.tail = lower.tail, log.p = log.p)
  }

  if (familytemp %in% c("conjugate", "none")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {

      eta <- eta / sqrt(1 + (sd^2)/(thresh / atan(thresh) - 1))

      eta <- pmin(pmax(eta, -thresh + .Machine$double.eps),
                  thresh - .Machine$double.eps)

      ptcauchy(eta, location = 0, scale = 1, tau = thresh,
               lower.tail = lower.tail, log.p = log.p)
    }

    fname <- familystats$name <- "conjugate"
  }
  else {
    stop("not yet implemented for families other than 'conjugate' and 'none'")
  }

  mu.eta <- function (eta) {
    pmax(dtcauchy(eta, location = 0, scale = 1, tau = thresh), .Machine$double.eps)
  }

  mu.eta.eta <- function (eta) {
    pmax(2 * dtcauchy(eta, location = 0, scale = 1, tau = thresh) * eta / (1 + eta^2), .Machine$double.eps)
  }

  structure(list(family = "binomial",
                 link = link,
                 linkfun = linkfun,
                 linkinv = linkinv,
                 mefamilyname = fname,
                 mefamily = familystats,
                 sigma = sqrt(thresh / atan(thresh) - 1),
                 melinkinv = melinkinv,
                 mu.eta = mu.eta,
                 mu.eta.eta = mu.eta.eta),
            class = "link")
}

