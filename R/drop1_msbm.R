
#' @exportS3Method stats::drop1
drop1.msbm <- function(object, scope,
                       component = c("stage", "alpha", "lambda"),
                       target = c("component", "predictor"),
                       test = c("none", "Rao", "LRT", "Chisq", "F"),
                       k = 2, ...) {
  component <- match.arg(component)
  target <- match.arg(target)
  test <- match.arg(test)

  if (identical(component, "alpha")) {

  }
  else if  (identical(component, "lambda")) {

  }


}
