#' Title
#'
#' @param cds
#' @param new_data
#' @param relative_expr
#'
#' @return
#' @export
#'
#' @examples
getGeneExpCur <- function(cds, new_data, relative_expr = T) {
  expressionFamily <- cds@expressionFamily
  expression_curve_matrix <- calVurveApply(X = cds,
                                           FUN = getFitModel,
                                           trend_formula = "~sm.ns(Time_num, df = 3)",
                                           expressionFamily = expressionFamily,
                                           disp_func = cds@dispFitInfo[['blind']]$disp_func,
                                           relative_expr = relative_expr,
                                           new_data = new_data)
  expression_curve_matrix <- as.matrix(do.call(rbind, expression_curve_matrix))
  row.names(expression_curve_matrix) <- row.names(fData(cds))
  return(expression_curve_matrix)
}
