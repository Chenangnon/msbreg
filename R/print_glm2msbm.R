
# @exportS3Method print glm2msbm
#' @exportS3Method stats::print
print.glm2msbm <- function(x, ...) {
  print.msbm (x, ...)
}

#setMethod ("print",
#           signature = "glm2msbm",
#           definition = print.glm2msbm)


# @exportS3Method summary glm2msbm
#' @exportS3Method stats::summary
summary.glm2msbm <- function(object, ...) {
  summary.msbm (object, ...)
}

#setMethod ("summary",
#           signature = "glm2msbm",
#           definition = summary.glm2msbm)
