
# Used by 'msbm_varcov'
safe_solve <- function (a, b, ...,
                        qr = is.qr(a),
                        chol = FALSE) {
  if (qr) {
    if (is.qr(a)) {
      out <- catch.conditions({
        solve.qr(a, b, ...)
      })
    }
    else {
      out <- catch.conditions({
        qr.solve(a, b, ...)
      })
    }
  }
  else if (chol) {
    out <- catch.conditions({
      chol2inv(x=a, ...)
    })
  }
  else {
    out <- catch.conditions({
      solve(a, b, ...)
    })
  }

  return(structure(out$value, warning = out$warning))
}
