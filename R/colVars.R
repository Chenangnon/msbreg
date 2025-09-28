#'
#' Row and Column Statistics
#'
#' Form row and column variance and standard deviation
#' for numeric arrays (or data frames).
#'
#' @param x an array of two or more dimensions, containing numeric,
#' complex, integer or logical values, or a numeric data frame.
#'
#' @param na.rm logical. Should missing values (including \code{NaN})
#' be omitted from the calculations?
#'
#' @param dims integer: Which dimensions are regarded as *rows* or *columns*
#' to operate over. See \link{colMeans}.
#'
#' @usage
#'
#' colVars(x, na.rm = FALSE, dims = 1L)
#'
#' colSds(x, na.rm = FALSE, dims = 1L)
#'
#' rowVars(x, na.rm = FALSE, dims = 1L)
#'
#' rowSds(x, na.rm = FALSE, dims = 1L)
#'
#' @details
#' These are convenience routines,
#' built as analogous of \link{colMeans} and \link{rowMeans}
#' for variance (\link{var}) and standard deviation (\link{sd}).
#' These functions are equivalent to use of \link{apply} with
#' \code{FUN = var} or \code{FUN = sd} with appropriate margins,
#' but are a lot faster. See \link{colMeans} for details
#' on \code{NaN} and \code{NA} handling.
#'
#' For a complex input, variance or standard deviation is
#' computed for the real and the imaginary parts separately
#' and assembled into a complex result.
#'
#' @export colVars
#' @export colSds
#' @export rowVars
#' @export rowSds
#'
#' @aliases colSds
#' @aliases rowVars
#' @aliases rowSds
#'
#' @return A numeric or complex array of suitable size,
#' or a vector if the result is one-dimensional.
#' The dimnames (or names for a vector result)
#' are taken from the original array \code{x}.
#'
#' @seealso \link{colMeans}.
#'
#' @examples
#' ##* Compute row and column variances for a matrix:
#' x <- cbind(x1 = 3, x2 = c(4:1, 2:5))
#' rowVars(x)
#' colVars(x)
#'
#' ##* an array
#' dim(UCBAdmissions)
#' rowMeans(UCBAdmissions)
#' rowVars(UCBAdmissions)
#'
#' rowMeans(UCBAdmissions, dims = 2)
#' rowVars(UCBAdmissions, dims = 2)
#'
#' colVars(UCBAdmissions)
#' colVars(UCBAdmissions, dims = 2)
#'
#' ##* complex case
#' x <- cbind(x1 = 3 + 2i, x2 = c(4:1, 2:5) - 5i)
#' x[3, ] <- NA; x[4, 2] <- NA
#' rowVars(x)
#' rowVars(x, na.rm = TRUE)
#'
#' colVars(x)
#' colVars(x, na.rm = TRUE)
#'
#' rowSds(x)
#' rowSds(x, na.rm = TRUE)
#'
#' colSds(x)
#' colSds(x, na.rm = TRUE)
#'
colVars <- function (x, na.rm = FALSE, dims = 1L) {
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.complex(x))
    return(colVars (Re(x), na.rm = na.rm, dims = dims) + (0+1i) *
             colVars (Im(x), na.rm = na.rm, dims = dims))

  xmeans <- colMeans(x, na.rm = na.rm, dims = dims)
  x2means <- colMeans(x^2, na.rm = na.rm, dims = dims)
  dn <- dim(x)
  n <- prod(dn[id <- seq_len(dims)])

  xvar <- (x2means - xmeans^2) * n / (n - 1)

  xvar
}

colSds <- function (x, na.rm = FALSE, dims = 1L) {
  if (is.data.frame(x))
    x <- as.matrix(x)

  if (is.complex(x)) {
    xsd <- sqrt(colVars (Re(x), na.rm = na.rm, dims = dims)) + (0+1i) *
      sqrt(colVars (Im(x), na.rm = na.rm, dims = dims))

    return(xsd)
  }

  xvar <- colVars(x, na.rm = na.rm, dims = dims)

  sqrt(xvar)
}

rowVars <- function (x, na.rm = FALSE, dims = 1L) {
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.complex(x))
    return(rowVars (Re(x), na.rm = na.rm, dims = dims) + (0+1i) *
             rowVars (Im(x), na.rm = na.rm, dims = dims))

  xmeans <- rowMeans(x, na.rm = na.rm, dims = dims)
  x2means <- rowMeans(x^2, na.rm = na.rm, dims = dims)
  dn <- dim(x)
  p <- prod(dn[-(id <- seq_len(dims))])

  xvar <- (x2means - xmeans^2) * p / (p - 1)

  xvar
}

rowSds <- function (x, na.rm = FALSE, dims = 1L) {
  if (is.data.frame(x))
    x <- as.matrix(x)

  if (is.complex(x)) {
    xsd <- sqrt(rowVars (Re(x), na.rm = na.rm, dims = dims)) + (0+1i) *
      sqrt(rowVars (Im(x), na.rm = na.rm, dims = dims))

    return(xsd)
  }

  xvar <- rowVars(x, na.rm = na.rm, dims = dims)
  sqrt(xvar)
}
