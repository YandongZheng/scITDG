#' Title
#'
#' @param object
#' @param celltype.list
#' @param p2p.deg.use
#' @param logfc.threshold
#' @param pvalue.threshold
#' @param celltype.use
#' @param time.points.use
#' @param time.points
#' @param sample.ncell.use
#' @param seed.use
#' @param n.cores
#' @param save.wd
#'
#' @return
#' @export
#'
#' @examples
CalculatescITDG <- function(object,
                            celltype.list,
                            p2p.deg.use,
                            logfc.threshold = 0.25,
                            pvalue.threshold = 0.05,
                            celltype.use,
                            time.points.use,
                            time.points,
                            sample.ncell.use = 200,
                            seed.use = NULL,
                            n.cores = 1,
                            save.wd = NULL) {

  if (is.null(save.wd)) save.wd <- "./"
  if (!dir.exists(save.wd)) dir.create(save.wd, recursive = TRUE)

  GeneExpCur.wd <- file.path(save.wd, "GeneExpCur")
  dir.create(GeneExpCur.wd, showWarnings = FALSE, recursive = TRUE)


  # Filter DEGs based on thresholds --------------------------------------------
  p2p.deg.list <- p2p.deg.use %>%
    dplyr::filter(abs(avg_log2FC) >= logfc.threshold,
                  p_val_adj <= pvalue.threshold) %>%
    dplyr::select(gene, cell.type) %>%
    dplyr::distinct()


  original_plan <- future::plan()
  on.exit(future::plan(original_plan), add = TRUE)

  if (n.cores > 1) {
    future::plan(future::multisession, workers = n.cores)
  }
  options(future.globals.maxSize = 5*1024^3)


  # Parallel analysis of all cell types-----------------------------------------
  exp_cur_files <- future_map(
    celltype.list,
    function(x) {
      result <- ProcessscITDG(
        p2p.deg.list = p2p.deg.list,
        cal.celltype = x,
        object = object,
        celltype.use = celltype.use,
        time.points.use = time.points.use,
        time.points = time.points,
        sample.ncell.use = sample.ncell.use,
        seed.use = seed.use,
        GeneExpCur.wd = GeneExpCur.wd
      )
      return(result)
    },
    .options = furrr_options(
      globals = c("ProcessscITDG",
                  "p2p.deg.list", "object", "celltype.use",
                  "time.points.use", "time.points", "sample.ncell.use",
                  "GeneExpCur.wd", "tissue"),
      packages = c("scITDG", "dplyr", "Biobase", "Seurat", "SeuratObject",
                   "magrittr", "VGAM", "parallel")  # 确保依赖包加载
    )
  ) %>% purrr::discard(is.null)



  # Expression curves were incorporated for all cell types----------------------
  if (length(exp_cur_files) == 0) stop("No valid cell types processed")

  # Combine expression curves from all cell types
  exp.cur.all <- exp_cur_files %>%
    purrr::map(~ {
      exp_cur <- readRDS(.x)  # 加载文件
      exp_cur
    }) %>%
    purrr::reduce(rbind)  # 使用reduce函数逐行合并


  # Calculate the distance matrix-----------------------------------------------
  mat <- as.matrix(exp.cur.all)
  row.dist <- stats::as.dist((1 - cor(t(mat), method = "spearman")) / 2)
  row.dist[is.na(row.dist)] <- 1
  saveRDS(row.dist, file.path(save.wd, "row.dist.rds"))

  return(list(
    exp.cur = exp.cur.all,
    row.dist = row.dist,
    object.gene = rownames(object),
    time.points = time.points,
    sample.ncell.use = sample.ncell.use
  ))
}
