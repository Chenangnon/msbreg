# @rdname score.test
# @export scoretable
scoretable <- function (object, ..., sort = TRUE, dmax = 10) {
  if ("score.test" %in% class(object))
    object <- list(object)
  if (length(match.call(expand.dots = F)$...) == 0)
    score.list <- object
  else score.list <- c(object, list(...))

  if (any(sapply(score.list, function(x) !("score.test" %in%
                                           class(x)))))
    stop("arguments must be 'score.test' objects")

}
