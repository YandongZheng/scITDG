#' Title
#'
#' @param object
#' @param celltype.use
#' @param celltype
#' @param time.points.use
#' @param time.points
#' @param cellnum.min
#' @param logfc.threshold
#' @param test.use
#'
#' @return
#' @export
#'
#' @examples
getP2PDEGs <- function(object, celltype.use = NULL, celltype = NULL,
                       time.points.use = NULL, time.points = NULL,
                       cellnum.min = 10, logfc.threshold = 0.25, test.use = "wilcox") {
  if (is.null(x = celltype.use)) {
    stop("Please specify celltype.use!")
  }
  if (is.null(x = time.points.use)) {
    stop("Please specify time.points.use!")
  }
  if (is.null(x = celltype)) {
    stop("Please specify celltype!")
  }
  if (is.null(x = time.points)) {
    stop("Please specify time.points!")
  }
  if (length(x = time.points) < 2) {
    stop("At least 2 ident must be specified in `time.points`")
  }
  Seurat::Idents(object = object) <- celltype.use
  cells <- which(object@active.ident == celltype)
  object.sub <- object[, cells]
  Seurat::Idents(object = object.sub) <- time.points.use

  filter.time.points <- which(table(object.sub@active.ident) >= cellnum.min) %>% names
  num.timepoints <- length(filter.time.points)

  if (num.timepoints > 1) {
    p2pdeg.com <- data.frame()
    for (i in 1:(num.timepoints - 1)) {
      time1 <- filter.time.points[i]
      for (j in (i + 1):num.timepoints) {
        time2 <- filter.time.points[j]
        p2pdeg <- object.sub %>%
          Seurat::FindMarkers(ident.1 = time2,
                              ident.2 = time1,
                              logfc.threshold = logfc.threshold,
                              test.use = test.use)
        if (nrow(p2pdeg) > 0) {
          p2pdeg$gene <- rownames(p2pdeg)
          p2pdeg$cell.type <- celltype
          p2pdeg$class <- paste0(time2, ".vs.", time1)
          p2pdeg %<>% subset(p_val_adj < 0.05)
          p2pdeg.com <- rbind(p2pdeg.com, p2pdeg)
        } else {
          next
        }
      }
    }
  } else {
    p2pdeg.com <- NULL
  }
  return(p2pdeg.com)
}
