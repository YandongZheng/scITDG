
#' Title
#'
#' @param object An seurat's object.
#' @param celltype.use A metadata column name to store the cell type. Currently only supported for class-level (i.e. non-quantitative) attributes.
#' @param celltype Cell type used for analysis.
#' @param time.points.use A metadata column name to store the time point. Currently only supported for class-level (i.e. non-quantitative) attributes.
#' @param time.points Variables of time point. Must be in the same order as the time point in the studied.
#' @param random.value Use random sampling. Default is true.
#' @param seed The number of seeds used for random sampling.
#' @param sample.ncell The number of random samples.
#'
#' @return
#' @export
#'
#' @examples

CreateITDGObject <- function(object,
                             celltype.use = NULL, celltype = NULL,
                             time.points.use = NULL, time.points = NULL,
                             random.value = TRUE,
                             seed = 176, sample.ncell = 200) {
  if (is.null(x = celltype.use)) {
    stop("Please specify celltype.use!")
  }
  if (is.null(x = time.points.use)) {
    stop("Please specify time.points.use!")
  }
  if (is.null(x = celltype)) {
    stop("Please specify celltype!")
  }
  if (length(x = celltype) > 1) {
    stop("Only 1 ident must be specified in `celltype`!")
  }
  if (is.null(x = time.points)) {
    stop("Please specify time.points!")
  }
  if (length(x = time.points) < 2) {
    stop("At least 2 ident must be specified in `time.points`!")
  }

  if (length(x = celltype) > 1) {
    celltype <- celltype[1]
  }

  Seurat::Idents(object = object) <- celltype.use
  cells <- which(object@active.ident == celltype)
  object.sub <- object[ , cells]
  Seurat::Idents(object = object.sub) <- time.points.use

  mat.com <- NULL
  meta.com <- NULL
  cellnum.com <- NULL
  for (t in time.points) {
    cells <- which(object.sub@active.ident == t)
    cellnum <- length(cells)
    cellnum.com %<>% c(cellnum)

    object.mat <- Seurat::GetAssayData(object.sub, slot = "counts")
    new.cellid <- paste0(celltype, "_", t, "_", c(1:sample.ncell))

    if (length(cells) > 0) {
      if (length(cells) > sample.ncell) {
        set.seed(seed)
        cellid <- sample(x = cells, size = sample.ncell, replace = F)
      } else {
        set.seed(seed)
        cellid <- sample(x = cells, size = sample.ncell, replace = T)
      }
      timepoint.mat <- object.mat[,cellid]
      colnames(timepoint.mat) <- new.cellid
    } else {
      if (random.value) {
        nrow = nrow(object.mat)
        ncol = sample.ncell
        max_num_fill <- ceiling(0.05 * nrow)
        timepoint.mat <- replicate(ncol, {
          num_fill <- sample(0:max_num_fill, 1)
          sample(c(rep(1, num_fill), rep(0, nrow - num_fill)))
        })
        colnames(timepoint.mat) <- new.cellid
        rownames(timepoint.mat) <- rownames(object.mat)
      } else {
        timepoint.mat <- matrix(0, nrow = nrow(object.mat), ncol = sample.ncell)
        colnames(timepoint.mat) <- new.cellid
        rownames(timepoint.mat) <- rownames(object.mat)
      }
    }
    mat.com %<>% cbind(timepoint.mat)

    timepoint.meta <- matrix(nrow = sample.ncell, ncol = 2) %>% data.frame
    timepoint.meta[,1] <- celltype
    timepoint.meta[,2] <-  t
    colnames(timepoint.meta) <-  c(celltype.use, time.points.use)
    rownames(timepoint.meta) <- new.cellid
    meta.com %<>% rbind(timepoint.meta)
  }

  if (0 %in% cellnum.com & !random.value) {
    mat.com <- mat.com + 1
  }

  fd <- as.data.frame(rownames(mat.com))
  colnames(fd) <- "gene_short_name"
  rownames(fd) <- fd$gene_short_name
  
  ITDGds <- methods::new("scITDGDataSet",
                         exprs = as.matrix(mat.com, "sparseMatrix"),
                         phenoData = Biobase::AnnotatedDataFrame(data = meta.com),
                         featureData = Biobase::AnnotatedDataFrame(data = fd),
                         lowerDetectionLimit = 0,
                         dispFitInfo = new.env(hash = TRUE),
                         expressionFamily = VGAM::negbinomial.size())
  return(ITDGds)
}
