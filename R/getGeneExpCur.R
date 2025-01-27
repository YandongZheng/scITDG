#' Title
#'
#' @param itdg_cds
#' @param new_data
#' @param relative_expr
#'
#' @return
#' @export
#'
#' @examples
getGeneExpCur <- function(itdg_cds, new_data, relative_expr = T) {
  expressionFamily <- itdg_cds@expressionFamily
  expression_curve_matrix <- calVurveApply(X = itdg_cds,
                                           FUN = getFitModel,
                                           trend_formula = "~sm.ns(Time_num, df = 3)",
                                           expressionFamily = expressionFamily,
                                           disp_func = itdg_cds@dispFitInfo[['blind']]$disp_func,
                                           relative_expr = relative_expr,
                                           new_data = new_data)
  expression_curve_matrix <- as.matrix(do.call(rbind, expression_curve_matrix))
  row.names(expression_curve_matrix) <- row.names(fData(itdg_cds))
  return(expression_curve_matrix)
}
