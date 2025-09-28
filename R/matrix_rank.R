# Function taken from package 'matrixcalc'
matrix.rank <- function (x, method = c("qr", "chol"), tol = 1e-07) {
  method <- method[1]
  if (method == "chol") {
    if (!is.square.matrix(x))
      stop("argument x is not a square matrix: required for method 'chol'")
    ans <- catch.conditions({
      chol(x, pivot = TRUE)
    })$value

    if (any(class(ans) %in% c("simpleError", "error",
                                 "condition", "try-error"))) {
      D <- svd(x)$d
      ans <- sum(abs(D) > tol)
    }
    else {
      ans <- attr(ans, "rank")
    }

  }
  else {
    ans <- catch.conditions({
      qr(x)
    })$value

    if (any(class(ans) %in% c("simpleError", "error",
                              "condition", "try-error"))) {
      D <- svd(x)$d
      ans <- sum(abs(D) > tol)
    }
    else {
      ans <- ans$rank
    }
  }
  return(ans)
}
