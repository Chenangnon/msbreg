
#' @export cloglog
#' @importFrom stats binomial
# @rdname probit
cloglog <- function (me.family = "conjugate") {
  familytemp <- substitute(me.family)
  if (!is.character(familytemp))
    familytemp <- deparse(familytemp)
  if(identical(me.family, "canonical"))
    me.family <- "conjugate"
  okFams <- c("conjugate", "lstab")
  link <- "cloglog"
  if (familytemp %in% okFams[-1]) {
    familystats <- make.me.family ('none')
  }
  else if (identical(me.family, "conjugate")) {
    familytemp <- "lstab" # log-positive stable
    familystats <- make.me.family ('none')
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
    if (log.p) {
      mu <- exp(mu)
    }
    log(-log(if(lower.tail) 1 - mu else mu))
  }

  linkinv <- function (eta, lower.tail = TRUE, log.p = FALSE) {

    pval <- pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)

    pval <- if(lower.tail) pval else 1 - pval

    if(log.p) log(pval) else pval
  }
  environment(linkfun) <- environment(linkinv) <- asNamespace("msbreg")

  if (familytemp %in% c("conjugate", "lstab")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {

      pval <- pmax(pmin(-expm1(-exp(eta / sqrt(1 + 6 * (sd^2)/(pi^2)))), 1 - .Machine$double.eps), .Machine$double.eps)

      pval <- if(lower.tail) pval else 1 - pval

      if(log.p) log(pval) else pval
    }
    environment(melinkinv) <- asNamespace("stats")

    fname <- familystats$name <- "lstab"
  }
  else {
    stop("not yet implemented for families other than 'lstab'")
  }

  mu.eta <- stats::binomial(link = 'cloglog')$mu.eta
  mu.eta.eta <- function (eta) {
    eta <- pmin(eta, 700)
    -expm1 (eta) * pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
  }

  structure(list(family = "binomial",
                 link = link,
                 linkfun = linkfun,
                 linkinv = linkinv,
                 mefamilyname = fname,
                 mefamily = familystats,
                 sigma = pi/sqrt(6),
                 melinkinv = melinkinv,
                 mu.eta = mu.eta,
                 mu.eta.eta = mu.eta.eta),
            class = "link")
}
