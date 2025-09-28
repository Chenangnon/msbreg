
#' @export sprobit
#' @importFrom stats binomial
# @rdname probit
# Consider the prostate cancer data in kim2008flexible for MSB
# Data from d2002biochemical
sprobit <- function (me.family = "conjugate", lambda = 0.75) {
  familytemp <- substitute(me.family)
  if (!is.character(familytemp))
    familytemp <- deparse(familytemp)
  if(identical(me.family, "canonical"))
    me.family <- "conjugate"
  okFams <- c("conjugate", "none")
  link <- "sprobit"
  if (familytemp %in% okFams[-1]) {
    familystats <- make.me.family ('none')
  }
  else if (identical(me.family, "conjugate")) {
    familytemp <- "conjugate"
    familystats <- make.me.family ("none")
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
    qsnorm(mu, lambda = lambda, lower.tail = lower.tail, log.p = log.p)
  }

  linkinv <- function (eta, lower.tail = TRUE, log.p = FALSE) {
    thresh <- -qsnorm(.Machine$double.eps, lambda = lambda)
    eta <- pmin(pmax(eta, -thresh), thresh)
    psnorm(eta, lambda = lambda, lower.tail = lower.tail, log.p = log.p)
  }

  if (familytemp %in% c("conjugate", "none")) {
    melinkinv <- function (eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
      thresh <- -qsnorm(.Machine$double.eps, lambda = lambda)

      eta <- eta / sqrt(1 + (sd^2))
      eta <- pmin(pmax(eta, -thresh), thresh)
      psnorm(eta, lambda = lambda, lower.tail = lower.tail, log.p = log.p)
    }

    familystats$name <- 'conjugate'
  }
  else {
    stop("not yet implemented for families other than 'conjugate'")
  }

  mu.eta <- function (eta) {
    dsnorm(eta, lambda = lambda)
  }

  #alpha <- (tanh(-lambda[!zeroscale]) + 1)/2 # Would be faster? No idea
  alpha <- 1 / (1 + exp(2*lambda))

  # Scale omega
  omega <- 1 / sqrt(alpha * (1 - alpha) + (1 - 2/pi) * (1 - 2 * alpha)^2)

  # Location
  mlocation <- - sqrt(2/pi) * (1 - 2 * alpha) * omega

  mu.eta.eta <- function (eta) {
    w <- alpha * omega
    w[eta > mlocation] <- (1 - alpha) * omega
    ((mlocation - eta) / (w * omega)) * dsnorm(eta, lambda = lambda)
  }

  structure(list(family = "binomial",
                 link = link,
                 linkfun = linkfun,
                 linkinv = linkinv,
                 mefamilyname = familystats$name,
                 mefamily = familystats,
                 sigma = structure(1, lambda = lambda),
                 melinkinv = melinkinv,
                 mu.eta = mu.eta,
                 mu.eta.eta = mu.eta.eta),
            class = "link")
}
