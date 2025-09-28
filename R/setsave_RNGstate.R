
setsave.RNGstate <- function() {
  expression({
    if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      runif(1)

    seed <- get0 ("seed", ifnotfound = NULL)

    if (is.null(seed))
      RNGstate <- get(".Random.seed", envir = .GlobalEnv)
    else {
      R.seed <- get(".Random.seed", envir = .GlobalEnv)
      set.seed(seed)
      RNGstate <- structure(seed, kind = as.list(RNGkind()))
      on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
    }
  })
}
