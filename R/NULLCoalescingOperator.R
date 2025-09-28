
#
# Reproducing the null coalescing operator
# %||% from ?base::Control
# Because it is not available in some versions of R
# Specifically, it is not available in R 4.2.3
# Since I want 'msbreg' to work from R 4.0.0 and
# latter versions, I include the code here
`%||%` <- function(x, y) if (is.null(x)) y else x
