
#' @rdname MCmsbm
sim.boot <- function (statistics, theta, data, nsim,
                      parallel, ncpus, cl, stat.names = NULL, ...) {

  if (!is.function(statistics))
    statistics <- get(statistics, mode = "function",
                      envir = parent.frame())

  raw <- boot::boot (data = data, statistic = statistics,
                     R = nsim, sim = "parametric", stype = "i", ran.gen = function(d, p) d,
                     mle = theta, parallel = parallel, ncpus = ncpus, cl = cl)

  outmat <- raw$t
  raw$t <- NULL

  if (NCOL(outmat) == length(stat.names))
    colnames(outmat) <- stat.names

  list(estimates = outmat,
       raw = raw)
}
