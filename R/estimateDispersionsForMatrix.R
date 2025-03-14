#' Title
#'
#' @param ITDGds
#' @param modelFormulaStr
#' @param relative_expr
#' @param min_cells_detected
#' @param remove_outliers
#'
#' @return
#' @export
#'
#' @examples
estimateDispersionsForMatrix <- function(ITDGds,
                                         modelFormulaStr = "~ 1",
                                         relative_expr = TRUE,
                                         min_cells_detected = 1,
                                         remove_outliers = TRUE) {

  dispModelName = "blind"
  stopifnot(is(ITDGds, "scITDGDataSet"))


  if (any(is.na(ITDGds$Size_Factor)))
    stop("NAs found in size factors. Have you called 'estimateSizeFactors'?")

  ITDGds@dispFitInfo = new.env(hash = TRUE)

  mu <- NA
  model_terms <- unlist(lapply(stringr::str_split(modelFormulaStr, "~|\\+|\\*"), stringr::str_trim))
  model_terms <- model_terms[model_terms != ""]
  progress_opts <- options()$dplyr.show_progress
  options(dplyr.show_progress = T)

  if (length(model_terms) > 1 || (length(model_terms) == 1 && model_terms[1] != "1")) {
    ITDGds_pdata <- dplyr::group_by(dplyr::select(tibble::rownames_to_column(pData(ITDGds)),
                                                  "rowname", .dots = model_terms),
                                    .dots = model_terms)
    disp_table <- as.data.frame(ITDGds_pdata %>% dplyr::do(
      disp_calc_helper_NB(ITDGds[,.$rowname], ITDGds@expressionFamily, min_cells_detected)))
  } else {
    ITDGds_pdata <- dplyr::group_by(dplyr::select(tibble::rownames_to_column(pData(ITDGds)),
                                                  "rowname"))
    disp_table <- as.data.frame(ITDGds_pdata %>% dplyr::do(
      disp_calc_helper_NB(ITDGds[,.$rowname], ITDGds@expressionFamily, min_cells_detected)))
  }


  if (!is.list(disp_table))
    stop("Parametric dispersion fitting failed, please set a different lowerDetectionLimit")
  disp_table <- subset(disp_table, is.na(mu) == FALSE)
  res <- parametricDispersionFit(disp_table)
  fit <- res[[1]]
  coefs <- res[[2]]

  if (remove_outliers) {
    CD <- cooks.distance(fit)
    cooksCutoff <- 4 / nrow(disp_table)
    outliers <- union(names(CD[CD > cooksCutoff]), setdiff(row.names(disp_table), names(CD)))
    res <- parametricDispersionFit(disp_table[row.names(disp_table) %in% outliers == FALSE,])
    fit <- res[[1]]
    coefs <- res[[2]]
  }
  names(coefs) <- c( "asymptDisp", "extraPois" )

  ans <- function(q) {
    coefs[1] + coefs[2] / q
  }
  attr(ans, "coefficients" ) <- coefs

  dfi <- list(disp_table = disp_table, disp_func = ans)
  ITDGds@dispFitInfo[[dispModelName]] <- dfi
  ITDGds
}
