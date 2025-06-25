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
  e <- new.env(parent = environment(FUN))
  
  pdata <- pData(X)
  for(n in names(pdata)) {
    assign(n, pdata[[n]], envir = e)
  }
  
  environment(FUN) <- e
  res <- apply(X = exprs(X), MARGIN = 1, FUN, ...)
  names(res) <- row.names(X)
  res
}
