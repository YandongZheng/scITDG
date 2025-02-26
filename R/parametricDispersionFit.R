#' Title
#'
#' @param disp_table
#' @param initial_coefs
#'
#' @return
#' @export
#'
#' @examples
parametricDispersionFit <- function(disp_table, initial_coefs = c(1e-6, 1)) {
  coefs <- initial_coefs
  iter <- 0
  while (TRUE) {
    residuals <- disp_table$disp / (coefs[1] + coefs[2] / disp_table$mu )
    good <- disp_table[which( (residuals > initial_coefs[1]) & (residuals < 10000) ),]
    fit <- glm(disp ~ I(1/mu), data = good, family = Gamma(link = "identity"), start = coefs)
    oldcoefs <- coefs
    coefs <- coefficients(fit)

    if (coefs[1] < initial_coefs[1]) {
      coefs[1] <- initial_coefs[1]
    }

    if (coefs[2] < 0) {
      stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
    }

    if (sum(log(coefs / oldcoefs )^2) < initial_coefs[1])
      break
    iter <- iter + 1

    if (iter > 10 ) {
      warning( "Dispersion fit did not converge." )
      break
    }
  }

  if ( !all(coefs > 0)) {
    stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
  }
  list(fit, coefs)
}
