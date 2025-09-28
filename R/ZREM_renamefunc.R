
# A function to rename columns of simulation results from MCmsbm
renamefunc <- function(x) {
  hasit <- grep('(Intercept).lambda', x = x, fixed = TRUE)
  hasit <- length(hasit)

  if (hasit) {
    x0 <- strsplit(x, split = '(Intercept).lambda', fixed = TRUE)[[1]]
    out <- paste0(x0[1], 'delta')
    if (length(x0) > 1) {
      out <- paste0(out, x0[2])
    }

    return(out)
  }

  hasit <- grep('.x', x = x, fixed = TRUE)
  hasit <- length(hasit)
  if (hasit) {
    x0 <- strsplit(x, split = '.x', fixed = TRUE)[[1]]
    out <- paste0(x0[1], '.beta')
    if (length(x0) > 1) {
      out <- paste0(out, x0[2])
    }

    return(out)
  }

  hasit <- grep('(Intercept)', x = x, fixed = TRUE)
  hasit <- length(hasit)
  hasdot <- grep('(Intercept).', x = x, fixed = TRUE)
  hasdot <- length(hasdot)

  if (!hasit & !hasdot)
    return(x)

  if (hasdot) {
    x0 <- strsplit(x, split = '(Intercept).', fixed = TRUE)[[1]]
    out <- paste0(x0[1], 'beta0')
    if (length(x0) > 1) {
      out <- paste0(out, x0[2])
    }

    return(out)
  }

  x0 <- strsplit(x, split = '(Intercept)', fixed = TRUE)[[1]]
  out <- paste0(x0[1], 'beta0')
  if (length(x0) > 1) {
    out <- paste0(out, x0[2])
  }
  else {
    out <- paste0(out, '1')
  }

  out
}


picksefunc <- function(x) {
  hasse <- grep('.se.', x = x, fixed = TRUE)
  hasse <- length(hasse)

  hasse > 0
}

pickconvfunc <- function(x) {
  hasse <- grep('.converged', x = x, fixed = TRUE)
  hasse <- length(hasse)

  hasse > 0
}
