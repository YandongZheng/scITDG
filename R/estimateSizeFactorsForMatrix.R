
#' Title
#'
#' @param ITDGds
#'
#' @return
#' @export
#'
#' @examples
estimateSizeFactorsForMatrix <- function(ITDGds) {
  CM <- exprs(ITDGds)
  CM <- round(CM)
  cell_total <- apply(CM, 2, sum)
  sfs <- cell_total / exp(mean(log(cell_total)))
  sfs[is.na(sfs)] <- 1
  ITDGds$Size_Factor <- sfs
  ITDGds
}
