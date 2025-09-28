
# @rdname test.separation
test.sep <- function (y, x,
                      weights = 1,# CURRENTLY ONLY USED when binary = TRUE for logit-Poisson transformation
                      offset = 0, # CURRENTLY NOT USED BY THE IMPLEMENTED ALGORITHM
                      binary = TRUE,
                      epsilon = 1e-08,
                      method = 'iterative_rectifier',
                      maxit = 500,
                      ...) {
  stopifnot(epsilon > 0, epsilon < 1)
  method <- method[1]
  stopifnot(is.finite(maxit))
  maxit <- max(maxit, 1)

  stopifnot(all(!is.na(cbind(y, weights, offset, x))))

  # Need adjustment to cope with offset > 0
  offset <- 0 # Setting offset = 0 for now ... !!!

  if(!(tolower(method) %in% c("iterative rectifier",
                              "iterative_rectifier",
                              "iterative-rectifier", "ir"))) {
    stop("only method = 'iterative_rectifier' is currently available")
  }

  nobs <- NROW(y)
  stopifnot(NROW(x) == nobs,
            length(weights) %in% c(1, nobs),
            length(epsilon) %in% c(1, nobs))


  # Use the logit-Poisson transformation if required
  if (binary) {
    y <- c(y, weights - y)
    x <- rbind(cbind(x, diag(nobs)), # Check this zero!!!
               cbind(0 * x, diag(nobs)))
    #    if (length(weights) > 1) {
    #      weights <- c(weights, weights)
    #    }
  }

  # Initialize the iterative rectifier algorithm
  y <- (y > 0) + 0 # y / weights # Not really useful
  zeroy <- y == 0
  K <- ceiling(sum(zeroy) / (epsilon^2)) # + 1
  W <- y + 1
  W[!zeroy] <- K
  U <- zeroy + 0

  # Initial step
  Reg <- catch.conditions({
    lm.wfit(x = x, y = U, w = W, offset = offset)
  })$value
  if (any(class(Reg) %in% c("simpleError", "error",
                            "condition", "try-error"))) {
    stop(paste0("the algorithm failled at iteration ", 0, ": ",
                Reg,
                "; check the inputs!"))
  }
  Uhat <- Reg$fitted.values
  zerahat <- abs(Uhat) < epsilon
  if (any(zerahat)) {
    Uhat[zerahat] <- rep(0, sum(zerahat))
  }
  # Rectify
  U <- pmax(Uhat, 0)

  # Iterate
  for (iter in 1:maxit) {
    Reg <- catch.conditions({
      lm.wfit(x = x, y = U, w = W, offset = offset)
    })$value
    if (any(class(Reg) %in% c("simpleError", "error",
                              "condition", "try-error"))) {
      stop(paste0("the algorithm failled at iteration ", iter, ": ",
                  Reg,
                  "; check the inputs!"))
    }
    Uhat <- Reg$fitted.values
    zerahat <- abs(Uhat) < epsilon
    if (any(zerahat)) {
      Uhat[zerahat] <- rep(0, sum(zerahat))
    }

    conv <- all(Uhat >= 0)
    if (conv) {
      break
    }

    # Rectify
    U <- pmax(Uhat, 0)
  }
  ABSresids <- abs(U - Uhat)

  # End
  iseparated <- abs(Uhat[1:nobs]) > epsilon
  separable <- any(iseparated)
  separated <- which(iseparated)

  out <- separable
  if (length(separated))
    attr(out, "separated") <- separated
  attr(out, "convergence") <- c(converged = conv,
                                niter = iter)
  attr(out, "resid") <- c(min = min(ABSresids[1:nobs]),
                          max = max(ABSresids[1:nobs]))
  attr(out, "epsilon") <- epsilon

  return(out)
}



















oldtest.sep <- function (y, x,
                      weights = 1,# CURRENTLY ONLY USED when binary = TRUE for logit-Poisson transformation
                      offset = 0, # CURRENTLY NOT USED BY THE IMPLEMENTED ALGORITHM
                      binary = TRUE,
                      epsilon = 1e-08,
                      method = 'iterative_rectifier',
                      maxit = 500,
                      ...) {
  stopifnot(epsilon > 0, epsilon < 1)
  method <- method[1]
  stopifnot(is.finite(maxit))

  stopifnot(all(!is.na(cbind(y, weights, offset, x))))

  # Need adjustment to cope with offset > 0
  offset <- 0 # Setting offset = 0 for now ... !!!

  if(!(tolower(method) %in% c("iterative rectifier",
                              "iterative_rectifier",
                              "iterative-rectifier", "ir"))) {
    stop("only method = 'iterative_rectifier' is currently available")
  }

  nobs <- NROW(y)
  stopifnot(NROW(x) == nobs,
            length(weights) %in% c(1, nobs),
            length(epsilon) %in% c(1, nobs))


  # Use the logit-Poisson transformation if required
  if (binary) {
    y <- c(y, weights - y)
    x <- rbind(cbind(x, diag(nobs)),
               cbind(0 * x, diag(nobs)))
    #    if (length(weights) > 1) {
    #      weights <- c(weights, weights)
    #    }
  }

  # Initialize the iterative rectifier algorithm
  y <- (y > 0) + 0 # y / weights # Not really useful
  U <- y - 1
  K <- ceiling(sum(U * U) / (epsilon^2)) + 1
  W <- y + 1
  zeroy <- y == 0
  W[!zeroy] <- K

  # Step 0
  Reg <- catch.conditions({
    lm.wfit(x = x, y = U, w = W, offset = offset)
  })$value
  if (any(class(Reg) %in% c("simpleError", "error",
                            "condition", "try-error"))) {
    stop(paste0("the algorithm failled at initialization: ",
                Reg,
                "; check the inputs!"))
  }
  Uhat <- Reg$fitted.values
  conv <- max(abs(U - Uhat)) < epsilon
  U[zeroy] <- pmin(Uhat[zeroy], 0)
  iter <- 1

  # Iterate
  while ((iter < maxit) & !conv) {
    Reg <- catch.conditions({
      lm.wfit(x = x, y = U, w = W, offset = offset)
    })$value
    if (any(class(Reg) %in% c("simpleError", "error",
                              "condition", "try-error"))) {
      stop(paste0("the algorithm failled after ", iter, " iterations: ",
                  Reg,
                  "; check the inputs!"))
    }
    Uhat <- Reg$fitted.values
    conv <- max(abs(U - Uhat)) < epsilon
    U[zeroy] <- pmin(Uhat[zeroy], 0)
    iter <- iter + 1
  }
  ABSresids <- abs(U - Uhat)

  # End
  iseparated <- abs(Uhat[1:nobs]) > epsilon
  separable <- any(iseparated)
  separated <- which(iseparated)

  out <- separable
  if (length(separated))
    attr(out, "separated") <- separated
  attr(out, "convergence") <- c(converged = conv,
                                niter = iter)
  attr(out, "resid") <- c(min = min(ABSresids[1:nobs]),
                          max = max(ABSresids[1:nobs]))
  attr(out, "epsilon") <- epsilon

  return(out)
}
