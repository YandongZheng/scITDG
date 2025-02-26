#' Title
#'
#' @param ITDGds
#' @param expressionFamily
#' @param min_cells_detected
#'
#' @return
#' @export
#'
#' @examples
disp_calc_helper_NB <- function(ITDGds, expressionFamily, min_cells_detected) {
  rounded <- round(exprs(ITDGds))
  nzGenes <- Matrix::rowSums(rounded > ITDGds@lowerDetectionLimit)
  nzGenes <- names(nzGenes[nzGenes > min_cells_detected])

  x <- t(t(rounded[nzGenes,]) / pData(ITDGds[nzGenes,])$Size_Factor)
  xim <- mean(1 / pData(ITDGds[nzGenes,])$Size_Factor)
  f_expression_mean <- Matrix::rowMeans(x)

  # For NB: Var(Y)=mu*(1+mu/k)
  f_expression_var <- Matrix::rowMeans((x - f_expression_mean)^2)
  disp_guess_meth_moments <- f_expression_var - xim * f_expression_mean
  disp_guess_meth_moments <- disp_guess_meth_moments / (f_expression_mean^2) #fix the calculation of k

  res <- data.frame(mu = as.vector(f_expression_mean),
                    disp = as.vector(disp_guess_meth_moments))
  res[res$mu == 0]$mu = NA
  res[res$mu == 0]$disp = NA
  res$disp[res$disp < 0] <- 0
  res <- cbind(gene_id = row.names(fData(ITDGds[nzGenes,])), res)
  res
}
