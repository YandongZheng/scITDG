#' Title
#'
#' @param x
#' @param trend_formula
#' @param expressionFamily
#' @param relative_expr
#' @param new_data
#' @param disp_func
#'
#' @return
#' @export
#'
#' @examples
getFitModel <- function(x, trend_formula, expressionFamily,
                        relative_expr, new_data, disp_func) {
  environment(calFitModel) <- environment()
  environment(responseMatrix) <- environment()
  
  model_fits <- calFitModel(x,
                            modelFormulaStr = trend_formula,
                            expressionFamily = expressionFamily,
                            relative_expr = relative_expr,
                            disp_func = disp_func)
  
  if (is.null(model_fits)) {
    expression_curve <- as.data.frame(matrix(rep(NA, nrow(new_data)), nrow = 1))
  } else {
    expression_curve <- as.data.frame(responseMatrix(list(model_fits),
                                                     new_data,
                                                     response_type = "response"))
  }
  
  colnames(expression_curve) <- row.names(new_data)
  expression_curve
}
