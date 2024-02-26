#' Title
#'
#' @param X
#' @param FUN
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calVurveApply <- function(X, FUN, ...) {
  parent <- environment(FUN)
  e1 <- new.env(parent = parent)
  multiassign(names(pData(X)), pData(X), envir = e1)
  environment(FUN) <- e1
  res <- apply(X = exprs(X), MARGIN = 1, FUN, ...)
  names(res) <- row.names(X)
  res
}
