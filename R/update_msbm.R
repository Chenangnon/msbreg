
#' Update and Re-evaluate a MSB Model Call
#'
#' Update a Multistage-Binomial (MSB) model fit.
#' \code{update.msbm} is a method for \link[stats]{update}.
#' It allows to update any (or all) of the components of MSB
#' model fit, or the fitting method.
#'
#' @param object an existing \code{"msbm"} fit object, typically returned
#' by \link[msbreg]{msbreg}.
#'
#' @param action character, a string indicating the updating action to
#' perform. The available choices include
#' \describe{
#' \item{\code{"model.frame"}}{ update the model frame related arguments
#' in the call that generated \code{object};}
#'
#' \item{\code{"apply.IC"}}{ alter the model fit object by adjusting the fit
#' to include intercept-correction (if not yet applied);}
#'
#' \item{\code{"cancel.IC"}}{ alter the model fit object by adjusting the fit
#' to exclude intercept-correction (if already applied).}
#' }
#'
#' When \code{action = "apply.IC"} or \code{action = "cancel.IC"}, all other
#' arguments are ignored, and irrespective of the argument \code{evaluate},
#' the output is always a model fit (and never an un-evaluated \code{call}).
#'
#' Cases \code{action = "apply.IC"} and \code{action = "cancel.IC"}
#' are IN-DEVELOPMENT. NOT READY FOR USAGE.
#'
#' All the remaining arguments are used only when \code{action = "model.frame"}.
#'
#' @param formula.,input.me.,stage.me.,alpha.formula.,lambda.formula.,kappa.formula. inputs
#' indicating adjustments to the corresponding argument of the original fit in
#' \code{object} (for instance, \code{formula.} is used to update the argument
#' \code{formula} of \link[msbreg]{msbreg}, \code{input.me.} is used for
#' \code{input.me} when provided, ...).
#'
#' @param ... further argument to \link[msbreg]{msbreg}. This can include one
#' or many of \code{link}, \code{data}, \code{weights}, \code{sample.weights},
#' \code{subset}, \code{na.action}, \code{start}, \code{control}, \code{y},
#' \code{frame}, and \code{method}; or any named argument of
#' \link[msbreg]{msbm.control}.
#'
#' @param evaluate logical, if \code{TRUE} evaluate the new \code{call} else
#' return the \code{call}.
#'
#' @details
#' The function mostly relies on \link[stats]{update.formula}: it works by
#' updating formulas in the \code{$call} component of \code{object}, and
#' identifying changes with respect to specified new input arguments.
#'
#' When \code{evaluate = TRUE}, the function might *fail* if any of the formulas
#' \code{formula.}, \code{alpha.formula.}, or \code{lambda.formula.} is updated
#' and \code{start} *is not accordingly adjusted* in the call to \code{update}.
#'
#' @return
#' The updated call when \code{action = "model.frame"} and
#' \code{evaluate = FALSE}. Otherwise, the fitted object.
#'
#' @exportS3Method stats::update
# @export update.msbm
#'
#' @seealso \link[msbreg]{msbreg} for fitting MSB models.
#'
#' @importFrom stats as.formula
#'
#' @examples
#' ##** Simulated two-stage logistic model example
#' data(test1data)
#' attr(test1data$y, "formula")
#'
#' # Estimating a minimum success probability, assumed in [0, 1)
#' MSBres <- msbreg (y/Total ~ x1 + offset(off1) | x2,
#'                    alpha.formula = ~ 1,
#'                    weights = Total,
#'                    data = test1data)
#' summary (MSBres)
#'
#' # The minimum success probability estimate
#' # is essentially zero (1.614e-05).
#' # Drop the related estimated intercept
#' MSBres1 <- update (MSBres, alpha.formula. = ~ 0)
#' summary (MSBres1)
#'
#' # Better fit by the AIC metric
#'
update.msbm <- function(object,
                        action = c("model.frame", "apply.IC", "cancel.IC"),
                        formula. = NULL,
                        input.me. = NULL,
                        stage.me. = NULL,
                        alpha.formula. = NULL,
                        lambda.formula. = NULL,
                        kappa.formula. = NULL,
                        ...,
                        evaluate = TRUE) {
  action <- match.arg(action)
  switch(action[1],
         apply.IC = {
           if (!object$IC)
             object <- apply.IC.msbm(object, ...)
           return(object)
         },
         cancel.IC = {
           if (object$IC)
             object <- cancel.IC.msbm(object, ...)
           return(object)
         },
         {
           eval(update_msbmbase())
           return(object)
         }
         )
}
#setMethod ("update",
#           signature = "msbm",
#           definition = update.msbm)

update_msbmbase <- function () {
  expression({
    call <- getCall(object)
    if (is.null(call))
      stop("need a 'msbm' class object with a call component")

    # Should the start argument be adjusted
    #adjuststart <- FALSE

    #* Main model formula
    if (!missing(formula.)) {
      call$formula <- update(formula(object), formula.)

      #adjuststart <- TRUE
    }

    #* Input measurement errors
    if (!missing(input.me.) & is.list(input.me.)) {
      newinputme <- names(input.me.)
      if (!is.null(newinputme)) {
        noldme <- length (call$input.me)
        if (noldme < 2) {
          call$input.me <- input.me.
        }
        else {
          oldinputme <- names(call$input.me)[-1]
          mnames <- match(newinputme, oldinputme)
          existing <- !is.na(mnames)
          if (any(existing)) {
            for (a in newinputme[existing])
              call$input.me[[a]] <- update(as.formula(call$input.me[[a]]),
                                           input.me.[[a]])
          }

          notexisting <- is.na(mnames)
          if (any(notexisting)) {
            for (a in newinputme[notexisting])
              call$input.me[[a]] <- as.formula(input.me.[[a]])
          }
        }
      }
    }

    #* Stage measurement errors (offsets)
    if (!missing(stage.me.))
      call$stage.me <- if(length(call$stage.me))
        update(as.formula(call$stage.me), stage.me.)
    else
      as.formula(stage.me.)

    #* Alpha component
    if (!missing(alpha.formula.)) {
      call$alpha.formula <- if (length(call$alpha.formula))
        update(as.formula(call$alpha.formula), alpha.formula.)
      else
        as.formula(alpha.formula.)

      #adjuststart <- TRUE
    }

    #* Lambda component
    if (!missing(lambda.formula.)) {
      call$lambda.formula <- if (length(call$alpha.formula))
        update(as.formula(call$lambda.formula), lambda.formula.)
      else
        lambda.formula.

      #adjuststart <- TRUE
    }

    #* kappa component
    if (!missing(kappa.formula.)) {
      call$kappa.formula <- if (length(call$alpha.formula))
        update(as.formula(call$kappa.formula), kappa.formula.)
      else
        kappa.formula.

      #adjuststart <- TRUE
    }

    #* Additional arguments to 'msbreg'
    extras <- match.call(expand.dots = FALSE)$...
    if (length(extras)) {
      existing <- !is.na(match(names(extras), names(call)))
      for (a in names(extras)[existing])
        call[[a]] <- extras[[a]]
      if (any(!existing)) {
        call <- c(as.list(call), extras[!existing])
        call <- as.call(call)
      }

      #if (adjuststart & !('start' %in% names(extras)[existing])) {
        # Drop the start argument in the original fit to avoid mismatch of the length parameter vector
        # call$start <- NULL # (I finally decided to revise 'msbm.fit' to deal with such mismatch)
      #}
    }

    object <- if (evaluate)
      eval(call, parent.frame())
    else call
  })
}
