# Taken from secr::setNumThreads of Murray Efford (2025)

#' @importFrom RcppParallel defaultNumThreads
#' @importFrom RcppParallel setThreadOptions
#'
setNumThreads <- function (ncores, ...) {
  current <- as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", ""))
  if (missing(ncores))
    ncores <- NULL
  if (is.na(current)) {
    if (is.null(ncores)) {
      ncores <- min(RcppParallel::defaultNumThreads(), 2)
    }
  }
  else if (is.null(ncores) || (ncores == current)) {
    return(current)
  }

  if (ncores > RcppParallel::defaultNumThreads()) {
    ncore0 <- ncores
    ncores <- RcppParallel::defaultNumThreads()
    warning(paste0("requested ncores (", ncore0, ") exceeds number available: using the maximum ", ncores))
  }
  if (ncores < 1) {
    ncores <- 1
    warning("specified ncores < 1: using 1 core")
  }
  ncores <- min(ncores, RcppParallel::defaultNumThreads())
  RcppParallel::setThreadOptions(ncores, ...)
  return(as.integer(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", "")))
}
