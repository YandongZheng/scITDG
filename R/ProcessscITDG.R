
#' Title
#'
#' @param object
#' @param celltype.use
#' @param p2p.deg.list
#' @param cal.celltype
#' @param time.points.use
#' @param time.points
#' @param sample.ncell.use
#' @param seed.use
#' @param GeneExpCur.wd
#'
#' @return
#' @export
#'
#' @examples
ProcessscITDG <- function(object,
                          celltype.use,
                          p2p.deg.list,
                          cal.celltype,
                          time.points.use,
                          time.points,
                          sample.ncell.use,
                          seed.use,
                          GeneExpCur.wd) {
  
  if (is.null(seed.use)) {
    seed.use = 176
  }
  
  cell_degs <- p2p.deg.list %>%
    dplyr::filter(cell.type == cal.celltype) %>%
    dplyr::pull(gene) %>%
    unique()
  
  if (length(cell_degs) < 2) {
    message(paste0("Skipping ", cal.celltype,
                   ": Insufficient DEGs (n < 2)"))
    return(NULL)
  }
  

  # Create ITDG object
  ITDGds.raw <- tryCatch({
    scITDG::CreateITDGObject(
      object = object,
      celltype.use = celltype.use,
      celltype = cal.celltype,
      time.points.use = time.points.use,
      time.points = time.points,
      random.value = TRUE,
      seed = seed.use,
      sample.ncell = sample.ncell.use
    )
  }, error = function(e) {
    message(paste0("Error in ITDG object creation for ",
                   cal.celltype, ": ", e$message))
    return(NULL)
  })
  
  if (is.null(ITDGds.raw)) {
    message(paste0("Skipping ", cal.celltype, ": ITDG object creation failed"))
    return(NULL)
  }
  
  # Preprocessing
  ITDGds.raw <- scITDG::estimateSizeFactorsForMatrix(ITDGds.raw)
  ITDGds.raw <- tryCatch({
    scITDG::estimateDispersionsForMatrix(ITDGds.raw)
  }, error = function(e) {
    message(paste0(
      "Error in dispersion estimation for ",
      cal.celltype, ": ", e$message))
    return(NULL)
  })
  
  if (is.null(ITDGds.raw) || length(cell_degs) < 2) {
    return(NULL)
  }
  
  # Define scITDGDataSet class--------------------------------------------------
  methods::setClass(
    "scITDGDataSet",
    contains = "ExpressionSet",
    slots = c(
      expressionFamily = "vglmff",
      lowerDetectionLimit = "numeric"
    ),
    prototype = prototype(
      methods::new("VersionedBiobase",
                   versions = c(
                     classVersion("ExpressionSet"),
                     scITDGDataSet = "1.0.0"
                   )
      )
    )
  )
  
  # Subset-specific cell type data
  ITDGds <- ITDGds.raw[cell_degs, ] %>%
    scITDG::AddTimePointData()
  
  # Generate new point-in-time data
  new.data <- data.frame(
    Time_num = seq(
      min(Biobase::pData(ITDGds)$Time_num),
      max(Biobase::pData(ITDGds)$Time_num),
      length.out = ncol(ITDGds)
    )
  )
  
  # Analytical expression curve
  exp_cur <- tryCatch({
    scITDG::getGeneExpCur(
      cds = ITDGds,
      new_data = new.data,
      relative_expr = TRUE
    ) %>%
      scITDG::ScaleExpCur(ITDGds, .)
  }, error = function(e) {
    message(paste0("Error in gene expression curve calculation for ",
                   cal.celltype, ": ", e$message))
    return(NULL)
  })
  
  if (is.null(exp_cur)) {
    return(NULL)
  }
  
  rownames(exp_cur) <- paste(
    tissue, gsub(" ", ".", cal.celltype), rownames(exp_cur), sep = ":"
  )
  
  # save expression curve
  exp_cur_file <- file.path(
    GeneExpCur.wd,
    paste0(gsub(" ", ".", cal.celltype), ".exp.cur.rds")
  )
  
  saveRDS(exp_cur, exp_cur_file)
  return(exp_cur_file)
}
