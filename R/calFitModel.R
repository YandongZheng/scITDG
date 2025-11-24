
#' Title
#'
#' @param x
#' @param modelFormulaStr
#' @param expressionFamily
#' @param relative_expr
#' @param disp_func
#' @param verbose
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calFitModel <- function(x, modelFormulaStr, expressionFamily,
                        relative_expr, disp_func = NULL, verbose=FALSE,
                        ...){
  modelFormulaStr <- paste0("f_expression", modelFormulaStr)
  orig_x <- x
  if (expressionFamily@vfamily %in% c("negbinomial.size")) {
    if (relative_expr) {
      x <- x / Size_Factor
    }
    f_expression <- round(x)

    if (is.null(disp_func) == FALSE) {
      expr_hint <- mean(round(orig_x))
      if (expr_hint > 0 && is.null(expr_hint) == FALSE) {
        disp_guess_fit <- disp_func(expr_hint)
        disp_guess <- disp_guess_fit
      } else {
        disp_guess <- NULL
      }

      if (is.null(disp_guess) == FALSE && disp_guess > 0 && is.na(disp_guess) == FALSE) {
        size_guess <- 1 / disp_guess
        expressionFamily <- negbinomial.size(size = 1 / disp_guess, ...)
      }
    }
  } else {
    f_expression <- x
  }
  tryCatch({
    FM_fit <- VGAM::vglm(as.formula(modelFormulaStr),
                         family = expressionFamily, epsilon = 1e-4)
    FM_fit
  }, error = function(e) {
    print(e);
    backup_expression_family <- NULL

    if (expressionFamily@vfamily %in% c("negbinomial.size")) {
      expr_hint <- max(round(orig_x))
      if (expr_hint > 0 && is.null(expr_hint) == FALSE) {
        disp_guess_fit <- disp_func(expr_hint)
        disp_guess <- disp_guess_fit
      } else {
        disp_guess <- NULL
      }
      backup_expression_family <- negbinomial()
    } else {
      backup_expression_family <- NULL
    }

    if (is.null(backup_expression_family) == FALSE) {
      test_res <- tryCatch({
        FM_fit <- VGAM::vglm(as.formula(modelFormulaStr),
                             family = backup_expression_family,
                             epsilon = 1e-1, checkwz = TRUE)
        FM_fit
      },
      error = function(e) {
        NULL
      })
      test_res
    } else {
      NULL
    }
  })
}

