#' Title
#'
#' @param X
#' @param FUN
#' @param ...
#'
#' @return
#' @export
#'
#' @examples

calVurveApply <- function(X, FUN, ...) {
  # 获取 phenoData
  pdata <- pData(X)
  
  # 创建数据列表
  data_list <- as.list(pdata)
  
  # 创建新环境（父环境为调用环境）
  e <- new.env(parent = environment(FUN))
  
  # 将数据添加到新环境
  list2env(data_list, envir = e)
  
  # 设置函数环境
  environment(FUN) <- e
  
  # 应用函数
  res <- apply(X = exprs(X), MARGIN = 1, function(gene_expr) {
    # 将基因表达数据添加到环境
    e$gene_expr <- gene_expr
    FUN(gene_expr, ...)
  })
  
  names(res) <- rownames(X)
  return(res)
}
