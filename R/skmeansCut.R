# spherical k-means
#' Title
#'
#' @param x
#' @param k
#' @param method
#'
#' @return
#' @export
#'
#' @examples
skmeansCut <- function(x, k, method = "pclust") {
  list(cluster = skmeans::skmeans(x, k, method = "pclust")$cluster)
}
