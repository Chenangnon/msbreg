#' @rdname identifiability.msbm
#' @exportS3Method msbreg::identifiability
identifiability.msbm.frame <- function (object, continuous.levels = 4, ...) {
  if (inherits(object, "msbm"))
    object <- model.frame (object)

  stopifnot(inherits(object, what = 'msbm.frame'))
  dims <- object$dims
  y <- object$y
  ntrials <- object$weights

  #* Check inner-identifiability for each component
  # Non-singularity (full rank) + Non-separation (overlap)
  ffun <- function(j) {
    n_colsj <- as.vector(dims$pj[j] + attr(object$stage.dictionary[[j]], "intercept"))
    if (n_colsj <= 1) {
      return(structure(c(identifiable = 1,
                         singular = 0,
                         rank.all = n_colsj,
                         rank.successes = n_colsj,
                         rank.failures = n_colsj,
                         separable = 0,
                         converged = 1,
                         niter = 0),
                       separated = integer(0L),
                       aliased = integer(0L),
                       resid = c(0, 0),
                       epsilon = 0))
    }

    Offsetj <- if (attr(object$stage.dictionary[[j]], "offset"))
      object$offset.matrix[, attr(object$stage.dictionary[[j]], "offset.index")]
    else 0

    Xj <- object$input.matrix[, object$stage.dictionary[[j]], drop = FALSE]
    if (attr(object$stage.dictionary[[j]], "intercept")) {
      Xj <- cbind(1, Xj)
    }
    Idj <- .identifiability.x.y (x = Xj, y = y,
                                 weights = ntrials,
                                 offset = Offsetj,
                                 method = "iterative_rectifier")
    separated <- unlist(attr(Idj, "separated"))
    if (length(separated) == 0)
      separated <- integer(0L)
    aliased <- unlist(attr(Idj, "aliased"))
    if (length(aliased) == 0)
      aliased <- integer(0L)

    return(structure(c(identifiable = Idj[1],
                       attr(Idj, "summary")),
                     separated = separated,
                     aliased = aliased,
                     resid = attr(Idj, "resid"),
                     epsilon = attr(Idj, "epsilon")))
  }
  iner_id <- lapply(1:dims$q, FUN = ffun)
  iner_separated <- lapply(iner_id, FUN = function(x) {attr(x, 'separated')})
  iner_aliased <- lapply(iner_id, FUN = function(x) {attr(x, 'aliased')})
  iner_resid <- sapply(iner_id, FUN = function(x) {attr(x, 'resid')})
  if (dims$q == 1)
    iner_resid <- cbind(iner_resid)
  iner_epsilon <- sapply(iner_id, FUN = function(x) {attr(x, 'epsilon')})
  iner_id <- do.call("cbind", iner_id)
  colnames(iner_id) <- colnames(iner_resid) <- names(iner_separated) <-
    names(iner_aliased) <- names(object$stage.dictionary)

  if (dims$d - attr(object$alpha.dictionary, "intercept")) {
    iner_id.alpha <- .identifiability.x.y (x = cbind(if (attr(object$alpha.dictionary, "intercept")) 1,
                                                     object$input.matrix[, object$alpha.dictionary]),
                                           y = y,
                                           weights = ntrials,
                                           offset = if (attr(object$alpha.dictionary, "offset"))
                                             object$offset.matrix[, attr(object$alpha.dictionary, "offset.index")]
                                           else 0,
                                           method = "iterative_rectifier")
    iner_separated$alpha <- unlist(attr(iner_id.alpha, "separated"))
    if (is.null(iner_separated$alpha))
      iner_separated$alpha <- integer(0L)
    iner_aliased$alpha <- attr(iner_id.alpha, "aliased")
    if (is.null(iner_aliased$alpha))
      iner_aliased$alpha <- integer(0L)
    iner_resid <- cbind(iner_resid, alpha = attr(iner_id.alpha, "resid"))
    iner_epsilon <- c(iner_epsilon, attr(iner_id.alpha, "epsilon"))
    iner_id.alpha <- c(identifiable = iner_id.alpha[1],
                       attr(iner_id.alpha, "summary"))
  }
  else if (attr(object$alpha.dictionary, "intercept")) {
    iner_separated$alpha <- integer(0L)
    iner_aliased$alpha <- integer(0L)
    iner_resid <- cbind(iner_resid, alpha = 0)
    iner_epsilon <- c(iner_epsilon, 0)
    iner_id.alpha <- c(identifiable = 1,
                       singular = 0,
                       rank.all = dims$d,
                       rank.successes = dims$d,
                       rank.failures = dims$d,
                       separable = 0,
                       convergence = 1,
                       niter = 0)
  }
  else {
    iner_id.alpha <- NULL
  }

  if (dims$r - attr(object$lambda.dictionary, "intercept")) {
    iner_id.lambda <- .identifiability.x.y (x = cbind(if (attr(object$lambda.dictionary, "intercept")) 1,
                                                      object$input.matrix[, object$lambda.dictionary]),
                                            y = y,
                                            weights = ntrials,
                                            offset = if (attr(object$lambda.dictionary, "offset"))
                                              object$offset.matrix[, attr(object$lambda.dictionary, "offset.index")]
                                            else 0,
                                            method = "iterative_rectifier")
    iner_separated$lambda <- unlist(attr(iner_id.lambda, "separated"))
    if (is.null(iner_separated$lambda))
      iner_separated$lambda <- integer(0L)
    iner_aliased$lambda <- attr(iner_id.lambda, "aliased")
    if (is.null(iner_aliased$lambda))
      iner_aliased$lambda <- integer(0L)
    iner_resid <- cbind(iner_resid, lambda = attr(iner_id.lambda, "resid"))
    iner_epsilon <- c(iner_epsilon, attr(iner_id.lambda, "epsilon"))
    iner_id.lambda <- c(identifiable = iner_id.lambda[1],
                        attr(iner_id.lambda, "summary"))
  }
  else if (attr(object$lambda.dictionary, "intercept")) {
    iner_separated$lambda <- integer(0L)
    iner_aliased$lambda <- integer(0L)
    iner_resid <- cbind(iner_resid, lambda = 0)
    iner_epsilon <- c(iner_epsilon, 0)
    iner_id.lambda <- c(identifiable = 1,
                        singular = 0,
                        rank.all = dims$r,
                        rank.successes = dims$r,
                        rank.failures = dims$r,
                        separable = 0,
                        convergence = 1,
                        niter = 0)
  }
  else {
    iner_id.lambda <- NULL
  }

  iner_id <- cbind(iner_id,
                   alpha = iner_id.alpha,
                   lambda = iner_id.lambda)
  iner.identifiability <- all(iner_id[1,] > 0)
  iner_id <- iner_id[-1,]

  #* Chech joint-identifiability across stages
  covs_labels <- object$auxy$covs_labels
  stage.covs <- object$tables$stage.covs
  # Build a matrix of offsets for all components
  Offsets <- sapply(object$stage.dictionary,
                    FUN = function (stg) {
                      offsetj <- attr(stg, "offset")
                      if (!offsetj)
                        return(rep(0, length.out = dims$nobs))
                      return(object$offset.matrix[,attr(stg, "offset.index")])
                    })
  Offsets <- cbind(Offsets,
                   alpha = if (attr(object$alpha.dictionary, "offset"))
                     rep(attr(object$alpha.dictionary, "offset.index"), length.out = dims$nobs)
                   else
                     rep(0, length.out = dims$nobs),
                   lambda = if (attr(object$lambda.dictionary, "offset"))
                     rep(attr(object$lambda.dictionary, "offset.index"), length.out = dims$nobs)
                   else
                     rep(0, length.out = dims$nobs))

  # Build a binary indicator matrix for all components and covariates
  component.covs <- rbind(stage.covs,
                          alpha = (covs_labels %in% attr(object$alpha.dictionary, "covs.names")) + 0,
                          lambda = (covs_labels %in% attr(object$lambda.dictionary, "covs.names")) + 0)

  # Find covariates that appear in only one component
  unique.covs <- which(colSums(component.covs) == 1)

  ffun <- function (j) {
    # Vars that are in stage j
    Vj <- which(stage.covs[j,] == 1)

    # Find those only in stage j
    Ij <- which(Vj %in% unique.covs)

    # At least one Var appears only in stage j?
    idt <- length(Ij) > 0
    cov.continuous <- 0
    if (idt) {
      # If yes, check if each is numeric and "continuous" [here: means at least
      # four or ('continuous.levels') distinct values]
      # Ij in an index for variables in stage j,
      # Associate it with Vj to find index for variables in all stages

      # Check that variables in Ij are numeric and "continuous"
      num.covs_labels <- object$auxy$num.covs_labels
      all.contrasts <- colnames(object$input.matrix)
      cov.continuous <- sapply(Vj[Ij], FUN = function(k) {
        # If the variable is not numeric
        if (!(covs_labels[k] %in% num.covs_labels)) {
          return(2) # At most binary column!
        }

        pos <- which(all.contrasts == covs_labels[k])
        uv <- unique(object$input.matrix[,pos])

        return(length(uv))
      })

      # Count the number of numeric and "continuous" variables
      cov.continuous <- sum(cov.continuous > continuous.levels)
    }

    # If Stage j has an offset specific to the stage,
    # we are okay (provided the same offset is not used anywhere else!).
    offset.continuous <- length(unique(Offsets[,j])) > continuous.levels

    return(c(cov.continuous, offset.continuous))
  }
  joint_stages <- sapply(1:dims$q, FUN = ffun)
  if (dims$q == 1) {
    joint_stages <- cbind(joint_stages)
  }
  colnames(joint_stages) <- names(object$stage.dictionary)
  rownames(joint_stages) <- c("continuous.predictor",
                              "continuous.offset")

  #joint.identifiability <- all(colSums(joint_stages) > 0)
  joint.identifiability <- sum(colSums(joint_stages) > 0) >= (dims$q - 1)

#  if (!joint.identifiability) {
#    joint.identifiability <- sum(colSums(joint_stages) > 0) >= (dims$q - 1)
#    # Check if only one stage lacks a continuous predictor/offset
#    if (joint.identifiability) {
#      if (length(object$alpha.dictionary) & length(object$lambda.dictionary)) {
#        al <- !all(object$alpha.dictionary == object$lambda.dictionary)
#      }
#      else if (!length(object$alpha.dictionary) & !length(object$lambda.dictionary)) {
#
#      }
#
#    }
#  }

  resid <- t(iner_resid)
  colnames(resid) <- c("min", "max")
  out <- structure(joint.identifiability & iner.identifiability,
                   summary  = iner_id,
                   continuous = joint_stages,
                   separated = if (length(unlist(iner_separated))) iner_separated,
                   aliased = if (length(unlist(iner_aliased))) iner_aliased,
                   resid = resid,
                   epsilon = max(iner_epsilon),
                   class = 'id.analysis')
  return(out)
}

#setMethod("identifiability",
#          signature(object = "msbm.frame"),
#          definition = identifiability.msbm.frame)
