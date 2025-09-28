#'
#' Row and Column Statistics
#'
#' Form row and column statistics for numeric matrix (or data frames).
#'
#' @param x an array containing numeric, integer or logical values,
#' or a numeric data frame.
#'
#' @param stat a function, or a character naming a function
#' that takes a vector input and returns a vector of statistics.
#' The default is \link{median}.
#'
#' If the function \code{stat} takes additional inputs, they should
#' be supplied as named arguments through \code{...}. None of such
#' arguments can however be named as \code{x} or \code{stat}.
#'
#' @param ... additional arguments passed to the function \code{stat}.
#' A common argument is \code{na.rm}, a logical indicating if missing
#' values (including \code{NaN}) should be omitted from the calculations.
#'
#' @usage
#'
#' colStats (x, stat = median, ...)
#'
#' rowStats (x, stat = median, ...)
#'
#' @details
#' The functions are convenience routines, built as analogous of
#' \link{colMeans} and \link{rowMeans} for any function appropriate
#' \code{stat}.
#'
#' These routines are strictly equivalent to the use of \link{apply}
#' with \code{FUN = stat} and appropriate margins.
#' As a result, the argument \code{simplify} of \link{apply} can be
#' passed through \code{...} to these routines.
#'
#' Arguments in \code{...} cannot have the same name as any of the
#' other arguments, or arguments of  \link{apply} such as \code{MARGIN}
#' or \code{FUN}. Care may also be needed to avoid partial matching
#' (to for instance \code{stat}, \code{MARGIN} or \code{FUN}).
#'
#' @return The output from \link{apply}.
#'
#' This depends on the presence of \code{simplify} in \code{...},
#' and the number of dimensions of \code{x}.
#' See \link{apply} for further details.
#'
#' @export colStats
#' @export rowStats
#' @aliases rowStats
#'
#' @seealso \link[msbreg]{colVars} and \link[msbreg]{rowVars}.
#'
#' @examples
#' ##** Compute row and column sums for a matrix:
#' x <- cbind(x1 = 3, x2 = c(4:1, 2:5))
#' dimnames(x)[[1]] <- letters[1:8]
#'
#' ##* Scalar statistic
#' require(msbreg)
#' rowStats (x, stat = min)
#' colStats (x, stat = max)
#' colStats (x, stat = sd)
#'
#' ##* Vector statistic
#' # Quantiles
#' rowStats (x, stat = quantile,
#'           probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))
#'
#' # Self defined function stat
#' unistats <- function(x) {
#'   c(mean = mean(x), min = min(x), med = median(x), max = max(x),
#'     sd = sd(x), IQR = IQR(x))
#' }
#'
#' colStats (x, stat = unistats)
#'
#' rowStats (x, stat = unistats)
#'
colStats <- function (x, stat = median, ...) {
  if (is.data.frame(x))
    x <- as.matrix(x)
  x <- cbind(x)
  if (!is.array(x) || length(dn <- dim(x)) < 2L)
    stop("'x' must be an array of at least two dimensions")

  if (is.character(family))
    stat <- get(stat, mode = "function",
                envir = parent.frame())

  apply(X = x, MARGIN = 2, FUN = stat, ...)
}

rowStats <- function (x, stat = median, ...) {
  if (is.data.frame(x))
    x <- as.matrix(x)
  x <- cbind(x)
  if (!is.array(x) || length(dn <- dim(x)) < 2L)
    stop("'x' must be an array of at least two dimensions")

  if (is.character(family))
    stat <- get(stat, mode = "function",
                envir = parent.frame())

  apply(X = x, MARGIN = 1, FUN = stat, ...)
}
