
toget_verbosewarning.frame <- function() {
  expression({
    if (verbose) {
      verbosewarning <- function (wmg, call., ...) {

        warning (paste0('In msbm.frame(...): ', wmg), call. = FALSE, ...)
      }
    }
    else {
      verbosewarning <- function (wmg, call., ...) {
        invisible()
      }
    }
  })
}

toget_main.frame <- function () {
  expression({
    mainframe <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights"),
               names(mainframe), 0L)
    mainframe <- mainframe[c(1L, m)]
    mainframe$na.action <- quote(na.pass)
    mainframe$drop.unused.levels <- drop.unused.levels
    mainframe[[1L]] <- quote(stats::model.frame)

    ########### Extract model frame for each stage ##########
    #* Stage frame
    chr_formula <- as.character(formula)

    xcharacters <- xcharacters0 <- strsplit(x = tail(chr_formula, 1),
                                            split = "|", fixed = TRUE)[[1]]
    if (length(chr_formula) == 3) {
      # Ensure that the response is not included in the left hand of formula when the '.' notation is used
      # Solution: add ‘- the response’ to each stage formula (when there is a response)
      # Useful when there is only one stage
      # For many stages, useful if stage has most variables, except a few.
      xcharacters <- sapply(xcharacters, FUN = function(x) {
        paste0(x, ' - ', chr_formula[2])
      })
    }

    stage.frame <- as.list(xcharacters)
    stage.envir <- parent.frame()
    adddatabeforeeval <- 0

    j <- 1
    formulaj <- as.formula(paste0('~', xcharacters[j]))
    environment(formulaj) <- environment(formula)
    stage.frame[[j]] <- mainframe
    stage.frame[[j]]$formula <- formulaj
    stage.frame_j1 <- catch.conditions({
      eval(stage.frame[[j]], envir = stage.envir)
    })$value

    if (any(class(stage.frame_j1) %in% c("simpleError", "error",
                                         "condition", "try-error"))) {
      # Review the 'envir' argument
      if (!any(c(is.list(data), is.environment(data)))) {
        if (!is.null(mcall$data))
          stage.envir <- eval(mcall$data, envir = parent.frame())
      }

      stage.frame_j1 <- catch.conditions({
        eval(stage.frame[[j]], envir = stage.envir)
      })$value

      if (any(class(stage.frame_j1) %in% c("simpleError", "error",
                                           "condition", "try-error"))) {

        stage.frame[[j]]$data <- data
        stage.frame_j1 <- eval(stage.frame[[j]], envir = stage.envir)
        adddatabeforeeval <- 1
      }
    }
    stage.frame[[j]] <- stage.frame_j1

    if (length(xcharacters) > 1) {
      for (j in seq_along(xcharacters)[-1]) {
        formulaj <- as.formula(paste0('~', xcharacters[j]))
        environment(formulaj) <- environment(formula)
        stage.frame[[j]] <- mainframe
        stage.frame[[j]]$formula <- formulaj
        if (adddatabeforeeval) {
          stage.frame[[j]]$data <- data
        }
        stage.frame[[j]] <- eval(stage.frame[[j]], envir = stage.envir)
      }
    }
    xcharacters <- xcharacters0

    #* Stage number of terms
    tj <- sapply(stage.frame,
                 FUN = function (x) {
                   length(attr(attr(x, 'terms'), 'term.labels'))
                 })

    #* Stage intercepts
    intercepts <- sapply(stage.frame,
                         FUN = function (x) {
                           attr(attr(x, 'terms'),
                                'intercept')
                         })
    #* Stage offsets
    ffun <- function (stg) {
      offsetj <- as.vector(stats::model.offset(stg))
      offsetj
    }
    offset.list <- lapply(stage.frame, FUN = ffun)
    stage.has.offset <- sapply(offset.list,
                               FUN = function (x) {
                                 # any(abs(x) > 0)
                                 !is.null(x)
                               })

    #* Drop empty stages
    # If there is no predictor in a stage (tj = 0), we must have a non zero offset.
    # If a stage only has an intercept, with no offset (i.e. offset = 0), remove it.
    # But we keep stages that only have an intercept, but have an offset as predictor
    # The latter allows to fit a model with offset-only stages.
    # This feature is intended for re-estimating stage intercepts while keeping
    # slopes fix to some already computed estimates.
    # The non constant part of each linear predictor can be calculated and supplied
    # as offset, and then re-estimate the slopes so that the mean of predicted
    # probabilities is equal to the mean of the binary responses.

    #* Remove empty stages
    keep <- (tj + stage.has.offset) > 0

    q <- sum(keep)
    if (q == 0) {
      verbosewarning ("no model term found in 'formula': at least one stage is required to fit a mutistage model")
      stop ("no stage detected: use 'glm' with 'family = binomial()' to fit a constant binomial glm")
    }
    if (!all(keep)) {
      # Remove empty stages
      stage.frame <- stage.frame[keep]
      tj <- tj[keep]
      intercepts <- intercepts[keep]
      offset.list <- offset.list[keep]
      stage.has.offset <- stage.has.offset[keep]
    }
    names(stage.frame) <- names(tj) <- names(intercepts) <-
      paste0('stage.', 1:q)

    #* Total number of terms and intercepts
    pt <- sum(tj)
    pi <- sum(intercepts)

    # Vector of covariate names (labels)
    stage.covs.name <- lapply(stage.frame,
                              FUN = function(stg) attr(stg, 'names'))
    all_covs.name <- unique(unlist(stage.covs.name))

    # Identify numeric covariates
    num_covs.name <- unique(unlist(
      lapply(stage.frame,
             FUN = function(stg) {
               attr(stg, 'names')[attr(attr(stg, 'terms'),
                                       'dataClasses') == "numeric"]
             })
    ))

    # Filter out columns added for weights or offsets
    if (length(all_covs.name)) {
      keepcovnames <- sapply(all_covs.name,
                             FUN = function(x) {
                               osplit <- strsplit (x = x, split = 'offset(', fixed = TRUE)[[1]]
                               wsplit <- strsplit (x = x, split = '(weights)', fixed = TRUE)[[1]]
                               (length(osplit) == 1) & (length(wsplit) == 1)
                             })
      all_covs.name <- all_covs.name[keepcovnames]
    }
    n_covs <- length(all_covs.name)
    stage.covs.name <- lapply(stage.covs.name,
                              FUN = function (x) {
                                if (length(x)) {
                                  x <- x[x %in% all_covs.name]
                                }
                                return(x)
                              })
    num_covs.name <- num_covs.name[num_covs.name %in% all_covs.name]
  })
}

toget_main.binary.tables <- function() {
  expression({
    # Binary table of predictors (rows) and stages (columns)
    if (n_covs == 0) {
      stage.covs_table <- matrix(, nrow = q, ncol = 0)
      dimnames(stage.covs_table) <- list(names(stage.frame),
                                         all_covs.name)
    }
    else if (n_covs == 1 && q == 1) {
      stage.covs_table <- as.matrix(1)
      dimnames(stage.covs_table) <- list(names(stage.frame),
                                         all_covs.name)
    }
    else {
      stage.covs_table <- sapply(stage.covs.name,
                                 FUN = function(stg.nm) {
                                   all_covs.name %in% stg.nm
                                 })
      stage.covs_table <- stage.covs_table + 0
      if (n_covs == 1) {
        stage.covs_table <- cbind(stage.covs_table)
      }
      else {
        stage.covs_table <- t(stage.covs_table)
      }

      colnames(stage.covs_table) <- all_covs.name
    }

    # Vector of model terms (labels)
    stage.terms.labels <- lapply(stage.frame,
                                 FUN = function(stg) {
                                   attr(attr(stg, 'terms'), 'term.labels')
                                 })
    all.terms.labels <- unique(unlist(stage.terms.labels))
    n_terms <- length(all.terms.labels)

    # Binary table of stages (rows) and terms (columns)
    if (n_terms == 0) {
      stage.terms_table <- matrix(, nrow = q, ncol = 0)
      dimnames(stage.terms_table) <- list(names(stage.frame),
                                         all.terms.labels)
    }
    else if (n_terms == 1 && q == 1) {
      stage.terms_table <- as.matrix(1)
      dimnames(stage.terms_table) <- list(names(stage.frame),
                                          all.terms.labels)
    }
    else {
      stage.terms_table <- sapply(stage.terms.labels,
                                  FUN = function(stgtms) {
                                    all.terms.labels %in% stgtms
                                  })
      if (n_terms == 1) {
        stage.terms_table <- rbind(stage.terms_table)
      }

      rownames(stage.terms_table) <- all.terms.labels
      stage.terms_table <- t(stage.terms_table) + 0
    }

    # Binary table of predictors (rows) and terms (columns)
    if (n_covs == 0 | n_terms == 0) {
      covs.terms_table <- matrix(, nrow = n_covs, ncol = n_terms)
    }
    else if (n_covs == 1 && n_terms == 1) {
      covs.terms_table <- as.matrix(1)
    }
    else {
      ffun <- function (k) {
        # Find the first stage where the term appears
        stg <- which(stage.terms_table[,k] == 1)[1]

        # Find the column position of the term 'trm' in stage 'stg'
        j <- which(all.terms.labels[k] == attr(attr(stage.frame[[stg]], 'terms'),
                                               'term.labels'))

        # Extract the vector of term factors in stage 'stg' corresponding to term 'trm'
        trm.fact <- attr(attr(stage.frame[[stg]], 'terms'),
                         'factors')[,j]

        # Extract predictors involved in term 'trm'
        preds <- rownames(attr(attr(stage.frame[[stg]], 'terms'),
                               'factors'))[trm.fact == 1]

        # Return a logical vector indicating predictors present in term 'trm'
        all_covs.name %in% preds
      }
      covs.terms_table <- sapply(1:n_terms, FUN = ffun)
      if (n_covs == 1) {
        covs.terms_table <- rbind(covs.terms_table)
      }
      covs.terms_table <- covs.terms_table + 0
    }
    dimnames(covs.terms_table) <- list(all_covs.name,
                                       all.terms.labels)

    # Vector of term orders
    if (n_terms == 0) {
      all.terms.orders <- numeric(0L)
    }
    else {
      ffun <- function (k) {
        # Find the first stage where the term appears
        stg <- which(stage.terms_table[,k] == 1)[1]

        # Find the column position of the term 'trm' in stage 'stg'
        j <- which(all.terms.labels[k] == attr(attr(stage.frame[[stg]], 'terms'),
                                               'term.labels'))
        # Term order
        trm.order <- attr(attr(stage.frame[[stg]], 'terms'), 'order')[j]

        # Return the term order
        trm.order
      }
      all.terms.orders <- sapply(1:n_terms, FUN = ffun)
    }
  })
}

toget_main.stage.matrix <- function () {
  expression({
    #* Stage model matrix
    input.matrix <- lapply(stage.frame,
                           FUN = stats::model.matrix,
                           na.action = stats::na.pass,
                           data = data)
    stage.term.order <- lapply(stage.frame, FUN = function(x) attr(attr(x, "terms"), "order"))
    names(stage.term.order) <- names(stage.frame)

    #* Total number of observations
    nobs <- sapply(input.matrix, FUN = NROW)
    if(!all(nobs == nobs[1])) {
      verbosewarning("extracted design matrices do not have the same number of rows!")
      stop("should not get here: something went wrong!")
    }
    nobs <- as.numeric(nobs[1])

    #* Check the length of offset at each stage
    ffun <- function (j) {
      offsetj <- offset.list[[j]]
      if (!(length(offsetj) %in% c(0, nobs))) {
        stop(gettextf("number of offsets is %d at stage %d, should equal %d (number of observations) or 0",
                      length(offsetj), j, nobs), domain = NA)
      }
      invisible(NULL)
    }
    sapply(1:q, FUN = ffun)

    #* Stage number of contrasts
    pj <- sapply(input.matrix,
                 FUN = function (stg) { NCOL(stg) }) - intercepts
    p <- sum(pj) + sum(intercepts)

    # Vector of contrasts (labels)
    stage.contrasts.labels <- lapply (input.matrix,
                                      FUN = function(stg) {
                                        contr <- colnames(stg)[attr(stg, 'assign') > 0]
                                        return(contr)
                                      })
    all.contrasts.labels <- unique(unlist(stage.contrasts.labels))
    n_contrasts <- length(all.contrasts.labels)

    # Check that each high order contrast has all involved first order contrasts in the model too
    ffun <- function (contr_label) {
      contrs <- strsplit (x = contr_label, split = ':', fixed = TRUE)[[1]]
      if (length(contrs) > 1) {
        all(contrs %in% all.contrasts.labels)
      }
      else {
        return(TRUE)
      }
    }
    check_contr <- sapply(all.contrasts.labels, FUN = ffun)
    if (!all(check_contr)) {
      check_contr <- which(!check_contr)
      verbosewarning(paste0(length(check_contr), ' higher order contrast(s) without the constituent first order contrasts in the model: ',
                            paste(all.contrasts.labels[check_contr], collapse = ', '), '.'))
      stop("any term in a higher order contrast must be a separate model term")
    }

    # Binary table of stages (rows) and contrasts (columns)
    if (n_contrasts == 0) {
      stage.contrasts_table <- matrix(, nrow = q, ncol = 0)
      dimnames(stage.contrasts_table) <- list(names(stage.frame),
                                              NULL)
    }
    else if (n_contrasts == 1 && q == 1) {
      stage.contrasts_table <- as.matrix(1)
      dimnames(stage.contrasts_table) <- list(names(stage.frame),
                                              all.contrasts.labels)
    }
    else {
      stage.contrasts_table <- sapply(stage.contrasts.labels,
                                      FUN = function(stgcontrs) {
                                        all.contrasts.labels %in% stgcontrs
                                      })
      if (n_contrasts == 1) {
        stage.contrasts_table <- rbind(stage.contrasts_table)
      }

      rownames(stage.contrasts_table) <- all.contrasts.labels
      stage.contrasts_table <- t(stage.contrasts_table) + 0
    }

    # Binary table of terms (rows) and contrasts (columns)
    if (n_covs == 0 | n_contrasts == 0) {
      term.contrasts_table <- matrix(, nrow = q, ncol = 0)
    }
    else if (n_covs == 1 && n_contrasts == 1) {
      term.contrasts_table <- as.matrix(1)
    }
    else {
      ffun <- function (k) {
        # Find the first stage where the contrast appears
        stg <- which(stage.contrasts_table[,k] == 1)[1]

        # Find the column position of the contrast in stage 'stg'
        j <- which(all.contrasts.labels[k] == stage.contrasts.labels[[stg]])

        # Find the position of the term from which the contrast derives
        ctr_assign <- attr(input.matrix[[stg]], 'assign')
        ctr_assign <- ctr_assign[ctr_assign > 0]
        trm_pos <- ctr_assign[j]

        # Find the term from which the contrast derives
        trm_j <- attr(attr(stage.frame[[stg]], 'terms'),
                      'term.labels')[trm_pos]

        # Return a logical vector indicating terms present in contrast
        all.terms.labels %in% trm_j
      }
      term.contrasts_table <- sapply(1:n_contrasts, FUN = ffun)
      if (n_terms == 1) {
        term.contrasts_table <- rbind(term.contrasts_table)
      }
      term.contrasts_table <- term.contrasts_table + 0
    }
    dimnames(term.contrasts_table) <- list(all.terms.labels,
                                           all.contrasts.labels)

    #* Full design matrix
    design.matrix <- input.matrix[[1]][, attr(input.matrix[[1]], 'assign') > 0, drop = FALSE]
    included <- colnames(design.matrix)
    if (q > 1) {
      for (j in 2:q) {
        not_in <- !(stage.contrasts.labels[[j]] %in% included)
        if (any (not_in)) {
          matj <- input.matrix[[j]][, attr(input.matrix[[j]], 'assign') > 0, drop = FALSE]
          design.matrix <- cbind(design.matrix, matj[, not_in])
          included <- c(included, stage.contrasts.labels[[j]][not_in])
        }
        if (all(all.contrasts.labels %in% included)) {
          break
        }
      }
      if (!identical(all.contrasts.labels, included)) {
        stop("should not get here: something went wrong!")
      }
      colnames(design.matrix) <- included
    }

    #* Full offset matrix
    offset.matrix <- do.call("cbind", offset.list)

    #* Index for extracting the contrasts of a given stage from the full design matrix
    #* Add term attributes to the index
    all.contrasts.labels <- colnames(design.matrix)
    ffun <- function (j) {
      if (length (stage.contrasts.labels[[j]])) {
        index <- sapply(stage.contrasts.labels[[j]],
                        FUN = function(ctrjk) {
                          which(ctrjk == all.contrasts.labels)
                        })

        stgj.trms <- sapply(stage.terms.labels[[j]],
                            FUN = function(lbl) {
                              which(lbl == all.terms.labels)
                            })
        assign <- attr(input.matrix[[j]], "assign")
        assign <- assign[assign > 0]
        assign <- stgj.trms[assign]
      }
      else {
        index <- assign <- integer(0)
      }

      mt <- attr(stage.frame[[j]], "terms")
      orderj <- attr(mt, "order")
      interceptj <- attr(mt, "intercept")

      attr(index, 'order') <- orderj
      attr(index, 'assign') <- assign
      attr(index, 'intercept') <- interceptj
      attr(index, 'offset') <- stage.has.offset[j] + 0
      attr(index, 'offset.index') <- if (stage.has.offset[j]) sum(stage.has.offset[1:j])
      attr(index, 'only.offset') <- (length(index) + interceptj) == 0

      return(index)
    }
    stage.dictionary <- lapply(1:q, FUN = ffun)
    names(stage.dictionary) <- names(stage.frame)
  })
}

toget_alphalambda.comps <- function () {
  expression({
    ########### Extract model frames for alpha #############
    # A function to extract matrix and attributes
    ffun <- function(aframe, design.matrix, offset.matrix) {
      dictionary <- numeric(0L)

      if (!is.empty.model(aframe)) {
        aoffset <- stats::model.offset(aframe)
        amatrix <- stats::model.matrix (aframe,
                                        na.action = stats::na.pass,
                                        data = data)
        ncols <- NCOL(amatrix)

        if (NROW(amatrix) != nobs) {
          if (NROW(amatrix) == 1) {
            amatrix <- matrix(c(amatrix), nrow = 1, ncol = ncols, byrow = TRUE) # Useless!!!
          }
          else {
            verbosewarning("extracted design matrices do not have the same number of rows!")
            stop("should not get here: something went wrong!")
          }
        }

        if(all(attr(amatrix, "assign") == 0)) {
          amatrix <- 1
          ncols <- 1
          attr(amatrix, 'intercept') <-
            attr(dictionary, 'intercept') <- 1
        }
        else {

          attr(amatrix, 'intercept') <-
            attr(attr(aframe, "terms"), "intercept")

          if (ncols > 1 & attr(amatrix, 'intercept')) {
            aattr <- attributes(amatrix)
            aattr$assign <- aattr$assign[-1]
            amatrix <- amatrix[,-1, drop = FALSE]
            aattr$dim <- attr(amatrix, "dim")
            aattr$dimnames <- attr(amatrix, "dimnames")
            attributes(amatrix) <- aattr
          }

          dictionary <- (NCOL(design.matrix) + 1):(NCOL(design.matrix) + ncols - attr(amatrix, 'intercept'))
          names(dictionary) <- colnames(amatrix)

          attr(dictionary, 'order') <- attr(amatrix, 'order')
          attr(dictionary, 'assign') <- attr(amatrix, 'assign')

          adups <- names(dictionary) %in% colnames(design.matrix)
          if (any(adups)) {
            addmat <- amatrix
            addmat <- if (!all(adups)) addmat[adups,]
            design.matrix <- cbind(design.matrix, addmat)

            dattr <- attributes(dictionary)
            dictionary <- sapply(names(dictionary),
                                 FUN = function(x) {
                                   which(x == colnames(design.matrix))[1]
                                 })
            attributes(dictionary) <- dattr
          }
          else {
            design.matrix <- cbind(design.matrix, amatrix)
          }
        }

        if (ncols - attr(amatrix, 'intercept') > 0) {
          alpha_covs.name <- unlist(attr(aframe, 'names'))
          keepcovnames <- sapply(alpha_covs.name,
                                 FUN = function(x) {
                                   osplit <- strsplit (x = x, split = 'offset(', fixed = TRUE)[[1]]
                                   wsplit <- strsplit (x = x, split = '(weights)', fixed = TRUE)[[1]]
                                   (length(osplit) == 1) & (length(wsplit) == 1)
                                 })
          alpha_covs.name <- alpha_covs.name[keepcovnames]
          attr(amatrix, 'covs.names') <- alpha_covs.name
        }
        else {
          attr(amatrix, 'covs.names') <- character(0)
        }
      }
      else {
        ncols <- 0
        aoffset <- stats::model.offset(aframe)
        amatrix <- matrix(, 0L, 0L)
        attr(amatrix, 'intercept') <- 0
        attr(amatrix, 'covs.names') <- character(0)
      }

      attr(dictionary, 'covs.names') <- attr(amatrix, 'covs.names')
      attr(dictionary, 'intercept') <- attr(amatrix, 'intercept')
      attr(dictionary, 'offset') <- attr(amatrix, 'offset') <- !is.null(aoffset)
      if (attr(dictionary, 'offset')) {
        offset.matrix <- cbind(offset.matrix, alpha = aoffset)
        attr(dictionary, 'offset.index') <-
          attr(amatrix, 'offset.index') <- NCOL(offset.matrix)
      }

      return(list(frame = aframe,
                  design.matrix = design.matrix,
                  offset.matrix = offset.matrix,
                  matrix = amatrix,
                  offset = aoffset,
                  dictionary = dictionary,
                  ncols = ncols))
    }

    alpha.frame <- mainframe
    alpha.frame$weights <- NULL
    alpha.frame$formula <- alpha.formula
    alpha.frame <- eval(alpha.frame, envir = environment(formula), enclos = parent.frame())
    alpha.frame <- ffun (alpha.frame, design.matrix, offset.matrix)
    design.matrix <- alpha.frame$design.matrix
    offset.matrix <- alpha.frame$offset.matrix
    alpha.matrix <- alpha.frame$matrix
    alpha.offset <- alpha.frame$offset
    alpha.dictionary <- alpha.frame$dictionary
    d <- alpha.frame$ncols
    alpha.frame <- alpha.frame$frame

    ########### Extract model frame for lambda #############
    lambda.frame <- mainframe
    lambda.frame$weights <- NULL
    lambda.frame$formula <- lambda.formula
    lambda.frame <- eval(lambda.frame, envir = environment(formula), enclos = parent.frame())
    lambda.frame <- ffun (lambda.frame, design.matrix, offset.matrix)
    design.matrix <- lambda.frame$design.matrix
    offset.matrix <- lambda.frame$offset.matrix
    lambda.matrix <- lambda.frame$matrix
    lambda.offset <- lambda.frame$offset
    lambda.dictionary <- lambda.frame$dictionary
    r <- lambda.frame$ncols
    lambda.frame <- lambda.frame$frame

    ########### Extract model frame for kappa #############
    kappa.frame <- mainframe
    kappa.frame$weights <- NULL
    kappa.frame$formula <- kappa.formula
    kappa.frame <- eval(kappa.frame, envir = environment(formula), enclos = parent.frame())
    kappa.frame <- ffun (kappa.frame, design.matrix, offset.matrix)
    design.matrix <- kappa.frame$design.matrix
    offset.matrix <- kappa.frame$offset.matrix
    kappa.matrix <- kappa.frame$matrix
    kappa.offset <- kappa.frame$offset
    kappa.dictionary <- kappa.frame$dictionary
    s <- kappa.frame$ncols
    kappa.frame <- kappa.frame$frame
  })
}

toget_sample.weights <- function() {
  expression({
    ########### Extract model frame for 'sample.weights' #############
    if (missing(sample.weights)) {
      sample.weights <- 1
    }
    else {
      char_sample.weights <- deparse(substitute(sample.weights))
      if (identical(char_sample.weights, "NULL")) {
        sample.weights <- 1
      }
      else {
        sample.weights_formula <- as.formula(paste0('~ 0 + ', char_sample.weights))
        environment(sample.weights_formula) <- environment(formula)
        sweightsframe <- mainframe
        sweightsframe$formula <- sample.weights_formula
        sweightsframe <- eval(sweightsframe, envir = environment(formula), enclos = parent.frame())

        if (is.empty.model(sweightsframe)) {
          sample.weights <- 1
        }
        else {
          sweights.matrix <- stats::model.matrix (sweightsframe,
                                                  na.action = stats::na.pass,
                                                  data = data)
          nsweights <- NCOL(sweights.matrix)
          if (nsweights == 0) {
            sample.weights <- 1
          }
          else if (nsweights == 1) {
            if (NROW(sweights.matrix) != nobs) {
              verbosewarning("column extracted from argument 'sample.weights' do not have the same number of rows as other design matrices!")
              stop("should not get here: something went wrong!")
            }
            sample.weights <- c(sweights.matrix)

            stopifnot(all(stats::na.omit(sample.weights) >= 0))
          }
          else {
            verbosewarning("extracted design matrices have more than on column!")
            stop("argument 'sample.weights' specified as a formula should have only a left side with only one column!")
          }
        }
      }
    }
  })
}

toget_inputstage.me <- function() {
  expression({
    ########### Extract model frame for measurement errors #############
    # By default, no measurement error
    includes.me <- FALSE
    names.input.me <- cov.with.me <- input.me.formula <-
      terms.me_list <- me.contrasts_table <- contrasts.has.me <-
      contrasts.with.me_labels <- contrasts.sd.matrix <- NULL
    me.frame <- list()
    n_me.cols <- 0

    if (length(input.me)) {
      # TO DO: allow higher order terms with user specified standard deviations
      # If not given, then we can have their standard deviations deduced from the included first order terms.
      # Currently, we only deduce ME for higher order terms from the included first order terms, using Taylor expansion arround the mean

      stopifnot(is.list(input.me))
      names.input.me <- names(input.me)
      if(length(names.input.me) != length(input.me)) {
        stop("all elements of 'input.me' must be named (see examples in '?msbm.frame')")
      }

      cname.design.matrix <- colnames(design.matrix)
      if (!all(names.input.me %in% cname.design.matrix)) {
        stop("all predictor named in 'input.me' must be a term included in 'formula'")
      }

      # Build an indicator of predictors with specified me
      cov.has.me <- sapply(all_covs.name,
                           FUN = function (x) x %in% names.input.me)
      all_allowed <- all_covs.name[cov.has.me] %in% num_covs.name
      if (!all(all_allowed)) {
        verbosewarning(paste0("measurement errors specified for non-numerical predictor(s):",
                              paste0(all_covs.name[cov.has.me][!all_allowed], collapse = ', ')))
        stop("measurement errors are only allowed for numerical predictors")
      }

      # Vector of labels of predictors with specified ME
      includes.me <- TRUE
      cov.with.me <- all_covs.name[cov.has.me]
      n_me.cols <- length(cov.with.me)

      # Reorder 'input.me' as 'cov.with.me'
      me.pos <- sapply(cov.with.me,
                       FUN = function (x) which(x == names.input.me))
      input.me <- input.me[me.pos]
      names.input.me <- names.input.me[me.pos]

      # Build a formula with all ME columns
      me.rhs <- sapply(input.me,
                       FUN = function(frml) {
                         chfrml <- as.character(frml)
                         chfrml[length(chfrml)]
                       })
      input.me.formula <- paste0(c("~ -1", me.rhs), collapse = '+')
      input.me.formula <- stats::as.formula(input.me.formula)
      environment(input.me.formula) <- environment(formula)

      # Build the standard deviation analogous to 'design.matrix'
      me.frame <- mainframe
      me.frame$formula <- input.me.formula
      #me.frame$weights <- NULL
      me.frame <- eval(me.frame, envir = environment(formula), enclos = parent.frame())
      # me.matrix <- me.frame
      me.matrix <- stats::model.matrix (me.frame,
                                        na.action = stats::na.pass,
                                        data = data)

      if(NCOL(me.matrix) > n_me.cols) {
        verbosewarning("only one column of measurement error standard deviations is allowed for one predictor!")
        stop("many columns of measurement error standard deviations specified for predictor(s)!")
      }
      if(NROW(me.matrix) != nobs) {
        verbosewarning("extracted matrix of measurement error standard deviations (from 'input.me') do not have the same number of rows as design matrices (from 'formula')!")
        stop("should not get here: something went wrong!")
      }
      if (any(stats::na.omit(me.matrix) < 0)) {
        # Change this warning to an error ?
        verbosewarning("extracted matrix of measurement error standard deviations (from 'input.me') have negative values: using the absolute values")
        me.matrix <- abs(me.matrix)
      }

      # List of terms per ME
      ffun <- function (me.cov) {
        # Find terms involving the covariate
        trms <- covs.terms_table[which(all_covs.name == me.cov),]
        trms <- all.terms.labels[which(trms == 1)]
        return(trms)
      }
      terms.me_list <- sapply(names.input.me, FUN = ffun)

      # Binary table of ME (rows) contrasts (columns)
      ffun <- function (trms) {
        # Find contrasts involving the covariate
        contrs <- colSums(term.contrasts_table[which(all.terms.labels %in% trms),, drop = FALSE])
        return(contrs)
      }
      me.contrasts_table <- sapply(terms.me_list, FUN = ffun)
      if (n_contrasts == 1) {
        me.contrasts_table <- cbind(me.contrasts_table)
      }
      else {
        me.contrasts_table <- t(me.contrasts_table)
      }
      contrasts.has.me <- colSums(me.contrasts_table) > 0
      contrasts.with.me_labels <- all.contrasts.labels[contrasts.has.me]

      # Matrix of standard errors for each contrast with ME
      ffun <- function (contr_label) {
        # Find the model term involved in this contrast
        trmk <- term.contrasts_table[, which(all.contrasts.labels == contr_label)]
        trmk <- which(trmk == 1)
        trm.label <- all.terms.labels[trmk]
        trm.order <- all.terms.orders[trmk]

        # Find all predictors involved in this term
        covsk <- all_covs.name[which(covs.terms_table[, trmk] == 1)]
        covk.has.me <- covsk %in% names.input.me
        nbcov <- length(covsk)
        nbcovk <- sum(covk.has.me)

        # Compute standard deviations
        if (nbcov == 1) {
          if (trm.order == 1) {
            me.sd <- me.matrix[, which(names.input.me == covsk)]
          }
          else {
            me.mu <- design.matrix[, which(cname.design.matrix %in% covsk)]^2
            me.sd <- me.matrix[, which(names.input.me %in% covsk)]^2
            if (trm.order == 2) {
              me.sd <- sqrt(pmax(0, me.sd * (4 * me.mu - me.sd)))
            }
            else if (trm.order == 3) {
              me.sd <- sqrt(pmax(0, 9 * me.sd * me.mu * (me.mu - me.sd)))
            }
            else {
              if (trm.order > 4) {
                verbosewarning("measurement error standard deviation for term of order higher than 4 derived assuming order = 4")
              }
              me.sd <- sqrt(pmax(0, 4 * (me.sd^2) * me.mu * (4 * me.mu - 9 * me.sd)))
            }
          }
        }
        else if (nbcov == 2) {
          contrs <- strsplit(x = contr_label, split = ':', fixed = TRUE)[[1]]
          me.mu <- design.matrix[, which(cname.design.matrix %in% contrs)]^2
          me.sd <- me.matrix[, which(names.input.me %in% covsk)]^2
          if (nbcovk == 1) {
            if (which(covk.has.me) == 1) {
              me.sd <- cbind(me.sd, 0)
            }
            else {
              me.sd <- cbind(0, me.sd)
            }
          }
          if (trm.order == 2) {
            me.sd <- sqrt(pmax(0, me.sd[, 1] * (me.sd[, 2] + me.mu[, 2]) + me.mu[, 1] * me.sd[, 2]))
          }
          else if (trm.order >= 3) {
            stop("measurement error not implemented for predictors involved in third or higher order terms")
          }
        }
        else {
          stop("measurement error not implemented for predictors involved in third order terms")
        }

        return(me.sd)
      }

      contrasts.sd.matrix <- sapply(contrasts.with.me_labels, FUN = ffun)
      rownames(contrasts.sd.matrix) <- rownames(design.matrix)

      # Build a ME dictionary for MSB model stages
      ffun <- function (stg.dict) {
        contrs <- all.contrasts.labels[stg.dict]
        stg_contr.has.me <- contrs %in% contrasts.with.me_labels
        stage.has.me <- any(stg_contr.has.me)
        contr.with.me <- me.index <- NULL
        if (stage.has.me) {
          contr.with.me <- which(stg_contr.has.me)
          me.index <- sapply(contrs[stg_contr.has.me], FUN = function (contrsj) {
            which(contrasts.with.me_labels == contrsj)
          })
          names(contr.with.me) <- names(me.index)
        }

        list(has.me = stage.has.me,
             contr.with.me = contr.with.me,
             me.index = me.index)
      }
      me.dictionary <- lapply(stage.dictionary, FUN = ffun)
    }
    else {
      me.dictionary <- lapply(1:q, FUN = function (j) {
        list(has.me = FALSE,
             contr.with.me = character(0),
             me.index = character(0))
      })
    }

    # Stage-wise measurement error (offsets)
    stage.has.me.offset <- logical(q)
    stage.me.formula <- stage.me.rhs <- NULL
    stage.me.frame <- me.offset.index <-
      vector(mode = "list", length = q)
    sd.offset.matrix <- NULL

    if (length(stage.me)) {
      stopifnot(methods::is(stage.me, class2 = "formula"))

      #* Stage me frame
      chr_meformula <- as.character(stage.me)
      if (length(chr_meformula) == 3) # Remove left hand side if present
        chr_meformula <- chr_meformula[-2]
      mecharacters <- strsplit (x = chr_meformula[2],
                                split = '|', fixed = TRUE)[[1]]
      stopifnot(length(mecharacters) <= length(keep)) # same length as the number of stages in "formula"?

      # Remove ME offset for dropped stages
      if (any(!keep)) {
        mecharacters <- mecharacters[keep[1:length(mecharacters)]]
      }

      # Model frames for ME offsets
      me.offsetj <- vector(mode = "list", length = q)
      for (j in seq_along(mecharacters)) {
        if(length(mecharacters[j]) > 0 & mecharacters[j] != "") {
          meformulaj <- as.formula(paste0('~ -1 + offset(', mecharacters[j], ')'))
          environment(meformulaj) <- environment(formula)
          stage.me.frame[[j]] <- mainframe
          stage.me.frame[[j]]$weights <- NULL
          stage.me.frame[[j]]$offset <- NULL
          stage.me.frame[[j]]$formula <- meformulaj
          stage.me.frame[[j]] <-
            eval(stage.me.frame[[j]], envir = environment(formula), enclos = parent.frame())
          stage.has.me.offset[j] <-
            attr(attr(stage.me.frame[[j]], "terms"), "offset") > 0

          if (stage.has.me.offset[j]) {
            me.offsetj[[j]] <- stats::model.offset (stage.me.frame[[j]])
            stage.has.me.offset[j] <- !is.null(me.offsetj[[j]])

            if (stage.has.me.offset[j]) {
              stage.has.me.offset[j] <- any(me.offsetj[[j]] > 0) |
                length(me.offsetj[[j]]) > 1

              if (stage.has.me.offset[j]) {
                stage.me.rhs <- c(stage.me.rhs, mecharacters[j])
              }
            }
          }
        }
      }

      # Join all ME offset formulas
      stage.me.formula <- paste0(c("~ -1", stage.me.rhs), collapse = '+')
      stage.me.formula <- stats::as.formula(stage.me.formula)

      # Check that each stage with me offset includes a numeric predictor or a non-null offset
      if (any(stage.has.me.offset)) {
        can.have.me.offset <- sapply(1:q, FUN = function(j) {
          !all(attr(stage.dictionary[[j]], 'offset') == 0) |
            any(stage.covs.name[[j]] %in% num_covs.name)
        })

        allright <- (can.have.me.offset + 0) >= (stage.has.me.offset + 0)
        if (!all(allright)) {
          verbosewarning(paste0("measurement errors specified for stage(s) with no numerical predictor or offset: stage ",
                                paste0(which(!allright), collapse = ', ')))
          stop("measurement errors are only allowed for stages with numerical predictors/offsets")
        }
      }

      # Build the standard deviation offset analogous to 'design.matrix'
      sd.offset.matrix <- NULL
      if (any(stage.has.me.offset)) {
        for (j in which(stage.has.me.offset)) {
          sd.offset.matrix <- cbind(sd.offset.matrix,
                                    rep(me.offsetj[[j]], length.out = nobs))
        }

        # Change this warning to an error ?
        if (any(stats::na.omit(sd.offset.matrix) < 0)) {
          verbosewarning("extracted matrix of measurement error standard deviation offsets (from 'stage.me') have negative values: using the absolute values")
          sd.offset.matrix <- abs(sd.offset.matrix)
        }
      }

      colnames(sd.offset.matrix) <- stage.me.rhs

      me.offset.index <- cumsum (stage.has.me.offset)
      me.offset.index[!stage.has.me.offset] <- 0
    }
    includes.me.offset <- any(stage.has.me.offset)

    # Insert 'me.dictionary' and me offsets into 'stage.dictionary'
    stage.dictionary <- lapply(1:q,
                               FUN = function (j) {
                                 outj <- stage.dictionary[[j]]
                                 attr(outj, 'has.me') <- me.dictionary[[j]]$has.me
                                 attr(outj, 'contr.with.me') <- me.dictionary[[j]]$contr.with.me
                                 attr(outj, 'me.index') <- me.dictionary[[j]]$me.index

                                 attr(outj, 'has.me.offset') <- stage.has.me.offset[[j]]

                                 if (attr(outj, 'has.me.offset'))
                                   attr(outj, 'me.offset.index') <- me.offset.index[[j]]

                                 return(outj)
                               })
    names(stage.dictionary) <- names(stage.frame)
  })
}

toget_offset.tables <- function() {
  expression({
    # Binary table of stages (rows) and offsets, if any (columns)
    noffs <- sum(stage.has.offset)
    if (noffs) {
      stage.offsets_table <- sapply(stage.dictionary,
                                    FUN = function(stg) {
                                      out <- numeric(noffs)
                                      if (!attr(stg, "offset")) {
                                        return(out)
                                      }
                                      out[attr(stg, "offset.index")] <- 1

                                      out
                                    })
      if (noffs == 1)
        stage.offsets_table <- cbind(stage.offsets_table)
      else
        stage.offsets_table <- t(stage.offsets_table)

      offnm <- lapply(stage.frame[stage.has.offset],
                      FUN = function(x) {
                        attr(x, "names")[attr(attr(x, "terms"), 'offset')]
                      })

      colnames(stage.offsets_table) <- offnm
      rownames(stage.offsets_table) <- rownames(stage.terms_table)

      if (attr(alpha.dictionary, 'offset')) {
        offnm <- c(offnm,
                   attr(alpha.frame, "names")[attr(attr(alpha.frame, "terms"), "offset")])
      }

      if (attr(lambda.dictionary, 'offset')) {
        offnm <- c(offnm,
                   attr(lambda.frame, "names")[attr(attr(lambda.frame, "terms"), "offset")])
      }

      if (attr(kappa.dictionary, 'offset')) {
        offnm <- c(offnm,
                   attr(kappa.frame, "names")[attr(attr(kappa.frame, "terms"), "offset")])
      }

      colnames(offset.matrix) <- offnm
    }
    else
      stage.offsets_table <- NULL

    # Binary table of stages (rows) and me.offsets, if any (columns)
    nmeoffs <- sum(stage.has.me.offset)
    if (nmeoffs) {
      stage.me.offsets_table <- sapply(stage.dictionary,
                                       FUN = function(stg) {
                                         out <- numeric(nmeoffs)
                                         if (!attr(stg, "has.me.offset")) {
                                           return(out)
                                         }
                                         out[attr(stg, "me.offset.index")] <- 1

                                         out
                                       })
      if (nmeoffs == 1)
        stage.me.offsets_table <- cbind(stage.me.offsets_table)
      else
        stage.me.offsets_table <- t(stage.me.offsets_table)

      colnames(stage.me.offsets_table) <- colnames(sd.offset.matrix)
      rownames(stage.me.offsets_table) <- rownames(stage.terms_table)
    }
    else
      stage.me.offsets_table <- NULL
  })
}

toget_response.frame <- function () {
  expression({
    ########### Extract model frame for the response #############
    mainformula <- paste(c(xcharacters,
                           if (d - attr(alpha.matrix, 'intercept')) as.character(alpha.formula)[2],
                           if (r - attr(lambda.matrix, 'intercept')) as.character(lambda.formula)[2],
                           if (s - attr(kappa.matrix, 'intercept')) as.character(kappa.formula)[2],
                           if (includes.me) as.character(input.me.formula)[2],
                           if (includes.me.offset) as.character(stage.me.formula)[2]),
                         collapse = '+')
    mainformula <- paste0(c(if (length(chr_formula) == 3) chr_formula[2],
                            ' ~ ', mainformula),
                          collapse = '')
    mainformula <- stats::as.formula(mainformula)
    environment(mainformula) <- environment(formula)
    mainframe <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"),
               names(mainframe), 0L)
    mainframe <- mainframe[c(1L, m)]
    mainframe$drop.unused.levels <- drop.unused.levels
    if (missing(na.action)) {
      mainframe$na.action <- quote(stats::na.omit)
    }
    mainframe$formula <- mainformula
    mainframe[[1L]] <- quote(stats::model.frame)
    mainframe <- eval(mainframe, envir = environment(formula), enclos = parent.frame())
    y <- model.response(mainframe, "any")
    yname <- colnames(mainframe)[1]
    ycolname <- get.ycolname (yname)

    if (attr(attr(mainframe, "terms"), "response") == 0) {
      y <- rep(NA, nobs - length(attr(mainframe, 'na.action')))
      nm <- rownames(mainframe)
      if (!is.null(nm))
        names(y) <- nm
      ycolname <- NA
    }
    else if (length(dim(y)) == 1L) {
      nm <- rownames(y)
      dim(y) <- NULL
      if (!is.null(nm))
        names(y) <- nm
    }
    weights <- as.vector(model.weights(mainframe))
    if (!is.null(weights) && !is.numeric(weights))
      stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
      stop("negative weights not allowed")
    if (is.null(weights)) {
      weights <- 1
    }

    if (attr(attr(mainframe, "terms"), "response") == 0) {
      ntrials <- weights
    }
    else if (NCOL(y) == 1) {
      if (is.factor(y))
        y <- y != levels(y)[1L]
      n <- rep.int(1, nobs)
      y[weights == 0] <- 0
      if (any(na.omit(y) < 0 | na.omit(y) > 1))
        stop("y values must be 0 <= y <= 1")
      my <- y * weights
      y <- round(my)
      if (any(abs(na.omit(my - y)) > 0.001))
        warning(gettextf("non-integer #successes in a %s glm!",
                         "binomial"), domain = NA)
      ntrials <- weights
    }
    else if (NCOL(y) == 2) {
      if (any(abs(y - round(y)) > 0.001))
        warning(gettextf("non-integer counts in a %s glm!",
                         "binomial"), domain = NA)
      y1 <- y[, 1L]
      n <- y1 + y[, 2L]
      y <- y1
      if (any(n0 <- n == 0))
        y[n0] <- 0
      y <- y * weights
      ntrials <- weights * n
    }
    else stop(gettextf("for the '%s' family, y must be a vector of 0 and 1's\nor a 2 column matrix where col 1 is #successes and col 2 is #failures",
                       "binomial"), domain = NA)
    if (length(ntrials) > 1)
      names(ntrials) <- names(y)
  })
}

toget_na.action.done <- function () {
  expression({
    # Apply na.action (if any) to other model components
    applied.na.action <- attr(mainframe, 'na.action')

    if (length(applied.na.action)) {
      nobs <- length(y)
      if (length(sample.weights) > 1) {
        sample.weights <- sample.weights[-applied.na.action]
        names(sample.weights) <- names(y)
      }

      design.matrix <- design.matrix[-applied.na.action, , drop = FALSE]

      if (NCOL(offset.matrix)) {
        offset.matrix <- offset.matrix[-applied.na.action, , drop = FALSE]
        rownames(offset.matrix) <- names(y)
      }

      if (includes.me) {
        contrasts.sd.matrix <- contrasts.sd.matrix[-applied.na.action, , drop = FALSE]
        rownames(contrasts.sd.matrix) <- names(y)
      }

      if (includes.me.offset) {
        sd.offset.matrix <- sd.offset.matrix[-applied.na.action, , drop = FALSE]
        rownames(sd.offset.matrix) <- names(y)
      }
    }

    #* Restore the 'na.action' in Global environment
    options(na.action = Glob.na.action)

  })
}

toget_all_par_bnames <- function () {
  expression({
    #* beta parameters
    all.bnames <- unlist(stage.contrasts.labels)
    if (q == 1) {
      betas.nm <- c(if (attr(stage.dictionary[[1]], "intercept")) "(Intercept)",
                    stage.contrasts.labels[[1]])
    }
    else {
      ffun <- function(j) {
        betaj.nm <- if (attr(stage.dictionary[[j]], "intercept")) paste0("(Intercept).", j)
        bnm <- stage.contrasts.labels[[j]]
        if (length(bnm)) {
          bnm <- sapply(bnm, FUN = function(nm) {
            if (sum(all.bnames %in% nm) > 1) {
              paste0(nm, '.stage.', j)
            }
            else {
              nm
            }
          })
          betaj.nm <- c(betaj.nm, bnm)
        }

        betaj.nm
      }
      betas.nm <- lapply(1:q, FUN = ffun)
    }

    gamma.nm <- delta.nm <- nu.nm <- NULL
    #* gamma parameters
    if (d > 0) {
      gamma.nm <- if (attr(alpha.matrix, 'intercept')) "(Intercept).alpha"

      if (d - attr(alpha.matrix, "intercept") > 0) {
        anm <- colnames(alpha.matrix)
        if (length(anm) != NCOL(alpha.matrix)) {
          stop("should not get here: something went wrong!")
        }

        anm <- sapply(anm, FUN = function(nm) {
          if (any(all.bnames %in% nm)) {
            paste0(nm, '.', 'alpha')
          }
          else {
            nm
          }
        })
        gamma.nm <- c(gamma.nm, anm)
      }
    }

    #* delta parameters
    if (r > 0) {
      delta.nm <- if (attr(lambda.matrix, 'intercept')) "(Intercept).lambda"
      if (r - attr(lambda.matrix, "intercept") > 0) {
        lnm <- colnames(lambda.matrix)
        if (length(lnm) != NCOL(lambda.matrix)) {
          stop("should not get here: something went wrong!")
        }

        lnm <- sapply(lnm, FUN = function(nm) {
          if (any(c(all.bnames, colnames(lambda.matrix)) %in% nm)) {
            paste0(nm, '.', 'lambda')
          }
          else {
            nm
          }
        })
        delta.nm <- c(delta.nm, lnm)
      }
    }

    #* nu parameters
    if (s > 0) {
      nu.nm <- if (attr(kappa.matrix, 'intercept')) "(Intercept).kappa"
      if (s - attr(kappa.matrix, "intercept") > 0) {
        knm <- colnames(kappa.matrix)
        if (length(knm) != NCOL(kappa.matrix)) {
          stop("should not get here: something went wrong!")
        }

        knm <- sapply(knm, FUN = function(nm) {
          if (any(c(all.bnames, colnames(kappa.matrix)) %in% nm)) {
            paste0(nm, '.', 'kappa')
          }
          else {
            nm
          }
        })
        nu.nm <- c(nu.nm, knm)
      }
    }

    # full names
    parnames <- c(unlist(betas.nm), gamma.nm, delta.nm, nu.nm)
    names(parnames) <- NULL
    parnames <- make.unique(parnames, sep = '.')
  })
}

toget_par_indices <- function() {
  expression({
    # Indicator of regression slopes in the vector of model parameters
    endpj <- as.numeric(cumsum(pj + intercepts))
    startpj <- c(0, endpj[-q]) + 1
    ffun <- function(j) {
      if (pj[j] == 0) {
        return(NULL)
      }

      id <- startpj[j]:endpj[j]
      if (intercepts[j]) {
        id <- id[-1]
      }

      return(id)
    }
    slope.id <- lapply(1:q, FUN = ffun)
    slope.id <- unlist(slope.id)

    # Indicator of regression intercepts in the vector of model parameters
    int.id <- sapply(1:q, FUN = function(j) {
      return(if (intercepts[j]) startpj[j] else 0)
    })
    eta.id <- which(int.id > 0)
    int.id <- int.id[eta.id]

    if (attr(alpha.dictionary, 'intercept')) {
      int.id <- c(int.id, endpj[q] + 1)
      eta.id <- c(eta.id, q + 1)
    }
    if (attr(lambda.dictionary, 'intercept')) {
      int.id <- c(int.id, endpj[q] +
                    (attr(alpha.dictionary, 'intercept'))
                  + 1)
      eta.id <- c(eta.id, q + 2)
    }
    if (attr(kappa.dictionary, 'intercept')) {
      int.id <- c(int.id, endpj[q] + (attr(alpha.dictionary, 'intercept')) +
                  (attr(lambda.dictionary, 'intercept')) + 1)
      eta.id <- c(eta.id, q + 3)
    }
    parindex <- list(slopes = slope.id,
                      intercepts = int.id,
                      eta.intercepts = eta.id)
  })
}

get.ycolname <- function (yname) {
  ycolname <- strsplit (x = yname,
                        split = ',', fixed = TRUE)[[1]][1]

  if (identical(ycolname, yname)) {
    ycolname <- strsplit (x = yname,
                          split = ')', fixed = TRUE)[[1]][1]

    if (identical(ycolname, yname)) {

      ycolname <- strsplit (x = yname,
                            split = '/', fixed = TRUE)[[1]][1]
      return(ycolname)
    }
  }

  ycolname <- strsplit (x = ycolname,
                        split = 'cbind(', fixed = TRUE)[[1]][2]

  return(ycolname)
}
