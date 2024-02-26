#' Title
#'
#' @param ITDGds
#' @param exp.cur
#' @param scale.max
#' @param scale.min
#'
#' @return
#' @export
#'
#' @examples
ScaleExpCur <- function(ITDGds, exp.cur, scale.max = 3, scale.min = -3) {
  exp.cur <- exp.cur[!apply(exp.cur, 1, sum) == 0, ]
  fitInfo <- ITDGds@dispFitInfo$blind
  if (is.null(fitInfo)) {
    stop(
      print("Error: No dispersion model named 'blind'. You must call estimateSizeFactorsForMatrix before calling this function.")
    )
  }
  coefs <- attr(fitInfo$disp_func, "coefficients" )
  if (is.null(exp.cur)) {
    stop(
      print("Error: No fit smooth spline curves named 'exp.cur'. You must call getGeneExpCur before calling this function.")
    )
  }
  ## NOTE: this converts to a dense matrix
  exp.cur <- log((1 + coefs["extraPois"] + 2 * coefs["asymptDisp"] * exp.cur +
                    2 * sqrt(coefs["asymptDisp"] * exp.cur *
                               ( 1 + coefs["extraPois"] + coefs["asymptDisp"] * exp.cur))) /
                   (4 * coefs["asymptDisp"])) / log(2)
  exp.cur <- exp.cur[!apply(exp.cur, 1, sd) == 0, ]
  exp.cur <- Matrix::t(scale(Matrix::t(exp.cur), center = TRUE))
  exp.cur <- exp.cur[is.na(row.names(exp.cur)) == FALSE, ]
  exp.cur[is.nan(exp.cur)] = 0
  exp.cur[exp.cur > scale.max] = scale.max
  exp.cur[exp.cur < scale.min] = scale.min
  return(exp.cur)
}
