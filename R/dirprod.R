#'
#' Direct Product of Matrices
#'
#' The function \code{dirprod} forms the direct (Kronecker) product
#' \eqn{\mathbf{D} = \mathbf{A} \otimes \mathbf{B}} of a \eqn{m \times n} matrix
#' \eqn{\mathbf{A}} and a \eqn{p \times q} matrix \eqn{\mathbf{B}}.
#' The matrix \eqn{\mathbf{D}} is of size \eqn{mp \times nq} and is defined as the
#' block matrix with \eqn{m} row-blocks and \eqn{n} column-blocks where the
#' \eqn{ij}th block is \eqn{a_{ij} \mathbf{B}}.
#' The function allows an iteration of the direct product operator
#' \eqn{\otimes} through \eqn{\mathbf{A} \otimes \mathbf{B} \otimes \mathbf{C} = }
#' \eqn{(\mathbf{A} \otimes \mathbf{B}) \otimes \mathbf{C}} =
#' \eqn{\mathbf{A} \otimes (\mathbf{B} \otimes \mathbf{C})}.
#' The function \code{dirpower} provides a shorthand for iterated direct products
#' of the same matrix:
#' \eqn{\mathbf{A}^{\otimes n} = \mathbf{A} \otimes \mathbf{A} \otimes \cdots \otimes \mathbf{A}}
#' where the matrix \eqn{\mathbf{A}} appears \eqn{n} times.
#'
#' @param A,B numeric matrices. A numeric vector is interpreted as a one-column
#' matrix (and coerced to a matrix using \link{as.matrix}).
#'
#' @param ... additional numeric matrices.
#'
#' @param n integer representing the power in the desired direct matrix power.
#'
#' @details
#' The direct product (also known as *Kronecker product*, *tensor product* or
#' *left direct product*) is very useful in several linear algebra applications,
#' including but not limited to the Sylvester equation
#' \eqn{\mathbf{A} \mathbf{X} + \mathbf{X} \mathbf{B} = \mathbf{C}} and the Lyapunov
#' equation \eqn{\mathbf{X} \mathbf{A} + \mathbf{A}^{*} \mathbf{X} = \mathbf{H}}
#' where \eqn{\mathbf{X}} is the unknown \insertCite{schacke2004kronecker}{msbreg}.
#' Of particular interest to the \code{msbreg} package, the use of the direct
#' product greatly simplifies matrix differentiation and the related derivation
#' of information matrix \insertCite{neudecker1969some,magnus1985matrix}{msbreg}.
#'
#' Some basic properties of the direct product \insertCite{schacke2004kronecker}{msbreg}:
#'
#' \itemize{
#'  \item associative: \eqn{(\mathbf{A} \otimes \mathbf{B}) \otimes \mathbf{C} = \mathbf{A} \otimes (\mathbf{B} \otimes \mathbf{C})};
#'  \item right-distributive: \eqn{(\mathbf{A} + \mathbf{B}) \otimes \mathbf{C} = \mathbf{A} \otimes \mathbf{C} + \mathbf{B} \otimes \mathbf{C}};
#' \item left-distributive: \eqn{\mathbf{A} \otimes (\mathbf{B} + \mathbf{C}) = \mathbf{A} \otimes \mathbf{B} + \mathbf{A} \otimes \mathbf{C}};
#' \item mixed product: \eqn{(\mathbf{A} \otimes \mathbf{B}) (\mathbf{C} \otimes \mathbf{D} = (\mathbf{A} \mathbf{C}) \otimes (\mathbf{B} \mathbf{D})}
#' if \eqn{\mathbf{A} \mathbf{C}} and \eqn{\mathbf{B} \mathbf{D}} can be computed given individual matrix dimensions.
#' }
#'
#' The call \code{dirpower (A, n)} is a simple shortcut for
#' \code{dirprod (A, A, ..., A)} involving \code{n} times the matrix \code{A}.
#' By definition, \code{dirpower (A, 1) = A}.
#' The convention \code{dirpower (A, 0) = 1} is used in light of
#' \eqn{\mathbf{A}^{\otimes n} \otimes \mathbf{A}^{\otimes m} = \mathbf{A}^{\otimes (n + m)}}.
#'
#' @usage
#' dirprod (A, B, ...)
#' dirpower (A, n)
#'
#' @return An \eqn{m \cdot p \times n \cdot q} matrix.
#'
#' @note
#' The function \code{dirprod} mainly reproduces the function
#' \link[dae]{mat.dirprod} in Chris Brien's package \code{dae}
#' \insertCite{brien2024dae}{msbreg}. The implementation here
#' allows for the direct product of more than two matrices.
#'
#' @export dirprod
#' @export dirpower
#' @aliases dirpower
#'
#' @references
#' \insertAllCited{}
#'
#' @examples
#' library(msbreg)
#'
#' A <- matrix(1:4, nrow = 2, ncol = 2)
#' A
#'
#' B <- matrix(c(1/3, 1/2, 1), nrow = 1, ncol = 3)
#' B
#'
#' # Basic direct product
#' dirprod (A, B)
#'
#' # The direct product
#' dirprod (A, B, c(0.5, 2))
#'
#' # is equivalent to the iterated direct products:
#'
#' dirprod (dirprod (A, B), c(0.5, 2))
#'
#' # and
#'
#' dirprod (A, dirprod (B, c(0.5, 2)))
#'
#' # The direct power
#' dirpower (A, 3)
#'
#' # is equivalent to the direct product:
#' dirprod (A, A, A)
#'

# Code from dae::mat.dirprod of Chris Brien
dirprod <- function (A, B, ...) {
  if (...length() > 0) {
    return(dirprod (dirprod (A, B), ...))
  }

  stopifnot(is.numeric(A), is.numeric(B))

  if(is.vector(A))
    A <- as.matrix(A)
  if(is.vector(B))
    B <- as.matrix(B)

  rA <- nrow(A)
  cA <- ncol(A)
  rB <- nrow(B)
  cB <- ncol(B)
  Aexp <- A[rep(1:rA, each = rB), rep(1:cA, each = cB)]
  Bexp <- eval(parse(text = paste("cbind(", paste(rep("B", cA),
                                                  collapse = ","), ")")))
  Bexp <- eval(parse(text = paste("rbind(", paste(rep("Bexp", rA),
                                                  collapse = ","), ")")))
  Aexp * Bexp
}

# Shortcut for 'dirprod (A, A, ..., A)' with n 'A' terms
dirpower <- function (A, n) {
  stopifnot(is.numeric(n))
  n <- n[1]
  stopifnot(floor(n) == n, is.finite(n))

  if (n == 0) {
    return(matrix(1, nrow = 1, ncol = 1))
  }

  if (n == 1)
    return(A)

  eval(parse(text = paste("dirprod (", paste(rep("A", n),
                                          collapse = ", "), ")")))
}
