#' Endometrial Cancer Data
#'
#' Histology grade and risk factors for 79 cases of endometrial cancer.
#' See \insertCite{heinze2002solution;textual}{msbreg} for a presentation of the study.
#'
#' @docType data
#'
#' @usage data(endometrial)
#'
#' @format An object of class \code{"data.frame"} with 79 rows and 4 variables:
#'
#' \describe{
#' \item{\code{NV}}{neovasculization with coding 0 for absent and 1 for present,}
#'
#' \item{\code{PI}}{pulsality index of arteria uterina,}
#'
#' \item{\code{EH}}{endometrium height,}
#'
#' \item{\code{HG}}{histology grade with coding 0 for low grade and 1 for high grade.}
#' }
#'
#' @keywords datasets
#'
#' @source Downloaded from \url{https://users.stat.ufl.edu/~aa/glm/data/}.
#'
#' @seealso See examples in \link[msbreg]{test.separation}.
#'
#' @references
#' \insertAllCited{}
#'
# \donttest{head(endometrial)}
"endometrial"

#endometrial <- as.data.frame(data.table::fread('https://users.stat.ufl.edu/~aa/glm/data/Endometrial.dat',
#                                               header = TRUE))
# save(endometrial, file = "data/endometrial.rda", version = 2,
#      compress = 'xz', compression_level = 9)
