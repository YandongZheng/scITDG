#' Title
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
AddTimePointData <- function(object) {
  length.out <- nrow(pData(object))
  time.num <- seq(0.1, 100, length.out = length.out)
  pData(object)$Time_num <- time.num
  return(object)
}
