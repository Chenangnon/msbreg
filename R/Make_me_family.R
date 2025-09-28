#'
#' Distribution Functions for Measurement Error Families
#'
#' This function is used with \link{link} functions in \link{msbreg}.
#' Given the name of a family, it returns a probability density function (pdf),
#' a cumulative distribution function (cdf), and a quantile (inverse cdf)
#' function (quf).
#'
#' @param me.family character; currently one of \code{"norm"}, \code{"t"},
#' \code{"logis"}, and \code{"bridge"}, corresponding to the Gaussian, Student t.
#' Logistic and Bridge distribution families, respectively; or \code{"none"}
#' (no particular distribution).
#'
#' @export make.me.family
#'
#' @return An object of class \code{me-family}, a list with components:
#' \item{pdf}{probability density function of the family.}
#' \item{cdf}{cumulative distribution function of the family.}
#' \item{quf}{quantile (inverse cdf) function of the family.}
#' \item{name}{a name to be used for the family.}
#'
#' When \code{me.family = 'none'}, the list elements \code{pdf},
#' \code{cdf}, and \code{qdf} are \code{NULL}.
#'
#' Otherwise, \code{pdf}, \code{cdf}, and \code{qdf} are
#' functions. The arguments of these three functions are
#' \item{eta}{numeric; vector of quantiles,}
#' \item{p}{numeric; vector of probabilities,}
#' \item{sd}{numeric; vector of standard deviations,}
#' \item{log,log.p}{logical; should the \code{log} of density/probability values
#' be returned? Defaults to \code{FALSE},}
#' \item{lower.tail}{logical, should the probability in the lower tail be
#' returned? Defaults to \code{TRUE}.}
#'
#' For the Student t family, the three functions have each the
#' additional argument
#' \item{df}{numeric; vector of degrees of freedom,
#' which defaults to \code{2.1}.}
#'
#' @seealso \link[msbreg]{link}, \link[msbreg]{msbreg}.
#'
make.me.family <- function(me.family) {
  switch(me.family,
         none = {
           return(structure(list(pdf = NULL, cdf = NULL, quf = NULL,
                                 name = me.family),
                            class = "me-family"))
         },
         norm = {
           pdf <- function(eta, sd = 1, log = FALSE) {
             dnorm (eta, sd = sd, log = log)
           }
           cdf <- function(eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
             thresh <- -qnorm(.Machine$double.eps)
             eta <- pmin(pmax(eta, -thresh), thresh)
             pnorm (eta, sd = sd, lower.tail = lower.tail, log.p = log.p)
           }
           quf <- function(p, sd = 1, lower.tail = TRUE, log.p = FALSE) {
             qnorm (p, sd = sd, lower.tail = lower.tail, log.p = log.p)
           }

           environment(pdf) <- environment(cdf) <- environment(quf) <- asNamespace("stats")
         },
         t = {
           pdf <- function(eta, df = 2.1, sd = 1, log = FALSE) {
             out <- dt (eta/sd, df = df, log = log)
             if (log)
               out - logb(sd)
             else
               out / sd
           }
           cdf <- function(eta, df = 2.1, sd = 1, lower.tail = TRUE, log.p = FALSE) {
             thresh <- -qnorm(.Machine$double.eps)
             eta <- pmin(pmax(eta/sd, -thresh), thresh)
             pt (eta, df = df, lower.tail = lower.tail, log.p = log.p)
           }
           quf <- function(p, df = 2.1, sd = 1, lower.tail = TRUE, log.p = FALSE) {
             out <- qt (p, df = df, lower.tail = lower.tail, log.p = log.p)
             out * sd
           }

           environment(pdf) <- environment(cdf) <- environment(quf) <- asNamespace("stats")
         },
         logis = {
           pdf <- function(eta, sd = 1, log = FALSE) {
             dlogis (eta, scale = sd * sqrt(3) / pi, log = log)
           }
           cdf <- function(eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
             plogis (eta, scale = sd * sqrt(3) / pi, lower.tail = lower.tail, log.p = log.p)
           }
           quf <- function(p, sd = 1, lower.tail = TRUE, log.p = FALSE) {
             qlogis (p, scale = sd * sqrt(3) / pi, lower.tail = lower.tail, log.p = log.p)
           }

           environment(pdf) <- environment(cdf) <- environment(quf) <- asNamespace("stats")
         },
         bridge = {
           pdf <- function(eta, sd = 1, log = FALSE) {
             dbridge (eta, scale = pi / sqrt(3*sd^2 + pi^2), log = log) # scale = phi
           }
           cdf <- function(eta, sd = 1, lower.tail = TRUE, log.p = FALSE) {
             pbridge (eta, scale = pi / sqrt(3*sd^2 + pi^2), lower.tail = lower.tail, log.p = log.p)
           }
           quf <- function(p, sd = 1, lower.tail = TRUE, log.p = FALSE) {
             qbridge (p, scale = pi / sqrt(3*sd^2 + pi^2), lower.tail = lower.tail, log.p = log.p)
           }

           environment(pdf) <- environment(cdf) <- environment(quf) <- asNamespace("msbreg")
         },
         {
           stop("not implemented for families other than 'norm', 'bridge' and 'logis'")
         })

  structure(list(pdf = pdf, cdf = cdf, quf = quf,
                 name = me.family),
            class = "me-family")
}
