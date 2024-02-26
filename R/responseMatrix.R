#' Title
#'
#' @param models
#' @param newdata
#' @param response_type
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
responseMatrix <- function(models, newdata = NULL, response_type = "response", cores = 1) {
  res_list <- mclapply(models, function(x) {
    if (is.null(x)) { NA } else {
      if (x@family@vfamily %in% c("negbinomial", "negbinomial.size")) {
        predict(x, newdata = newdata, type = response_type)
      } else if (x@family@vfamily %in% c("uninormal")) {
        predict(x, newdata = newdata, type = response_type)
      }
      else {
        10^predict(x, newdata = newdata, type = response_type)
      }
    }
  }, mc.cores = cores)

  res_list_lengths <- lapply(res_list[is.na(res_list) == FALSE],
                             length)
  stopifnot(length(unique(res_list_lengths)) == 1)
  num_na_fits <- length(res_list[is.na(res_list)])
  if (num_na_fits > 0) {
    na_matrix <- matrix(rep(rep(NA, res_list_lengths[[1]]),
                            num_na_fits), nrow = num_na_fits)
    row.names(na_matrix) <- names(res_list[is.na(res_list)])
    non_na_matrix <- Matrix::t(do.call(cbind, lapply(res_list[is.na(res_list) ==
                                                                FALSE], unlist)))
    row.names(non_na_matrix) <- names(res_list[is.na(res_list) ==
                                                 FALSE])
    res_matrix <- rbind(non_na_matrix, na_matrix)
    res_matrix <- res_matrix[names(res_list), ]
  }
  else {
    res_matrix <- Matrix::t(do.call(cbind, lapply(res_list, unlist)))
    row.names(res_matrix) <- names(res_list[is.na(res_list) ==
                                              FALSE])
  }
  res_matrix
}
