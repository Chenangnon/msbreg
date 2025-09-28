
# Use a truncated cauchy distribution within +- 3744
#' @export cauchit
#' @importFrom stats binomial
#' @importFrom stats pcauchy
#' @importFrom stats qcauchy
#' @importFrom stats dcauchy
# @rdname probit
# sums of independent Cauchy random variables will be Cauchy.
cauchit <- function (me.family = "conjugate") {
  familytemp <- substitute(me.family)
  if (!is.character(familytemp))
    familytemp <- deparse(familytemp)
  if(identical(me.family, "canonical"))
    me.family <- "conjugate"
  okFams <- c("conjugate", "none")
  link <- "cauchit"
  if (familytemp %in% okFams[-1]) {
    familystats <- make.me.family ("none")
  }
  else if (identical(me.family, "conjugate")) {
    familytemp <- "none"
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

  if (is.null(familystats$name))
    familystats$name <- familytemp

  linkfun <- function(mu, lower.tail = TRUE, log.p = FALSE) {
    qcauchy(mu, location = 0, scale = 1,
           lower.tail = lower.tail, log.p = log.p)
  }

  linkinv <- function (eta, lower.tail = TRUE, log.p = FALSE) {
    thresh <- -qcauchy(.Machine$double.eps)
    eta <- pmin(pmax(eta, -thresh), thresh)

    pcauchy(eta, location = 0, scale = 1,
            lower.tail = lower.tail, log.p = log.p)
  }
  environment(linkfun) <- environment(linkinv) <- asNamespace("stats")

  if (familytemp %in% c("conjugate", "none")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
      thresh <- -qcauchy(.Machine$double.eps)
      eta <- pmin(pmax(eta, -thresh), thresh)

      pcauchy(eta, location = 0, scale = 1,
              lower.tail = lower.tail, log.p = log.p)
    }
    environment(melinkinv) <- asNamespace("msbreg")
  }
  else {
    stop("not yet implemented for families other than 'conjugate' and 'none'")
  }

  mu.eta <- function (eta) {
    pmax(dcauchy(eta, location = 0, scale = 1), .Machine$double.eps)
  }

  mu.eta.eta <- function (eta) {
    pmax(2 * dcauchy(eta, location = 0, scale = 1) * eta / (1 + eta^2), .Machine$double.eps)
  }

  structure(list(family = "binomial",
                 link = link,
                 linkfun = linkfun,
                 linkinv = linkinv,
                 mefamilyname = familystats$name,
                 mefamily = familystats,
                 sigma = NaN,
                 melinkinv = melinkinv,
                 mu.eta = mu.eta,
                 mu.eta.eta = mu.eta.eta),
            class = "link")
}

