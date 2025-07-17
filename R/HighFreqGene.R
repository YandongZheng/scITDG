

#' Title
#'
#' @param deg.list 
#' @param moduleID 
#' @param top.num 
#' @param y_l 
#' @param text.size 
#' @param anno.cols 
#'
#' @return
#' @export
#'
#' @examples
HighFreqGene <- function(deg.list, 
                         moduleID, 
                         top.num = 50, 
                         y_l = -2, 
                         text.size = 2, 
                         anno.cols = NULL) {
  # Check if the required parameters are provided
  if (missing(deg.list) || missing(moduleID)) {
    stop("The parameters 'deg.list' and 'moduleID' are required!")
  }
  
  # Data preprocessing
  data.table::setDT(deg.list)[, c("column1", "column2", "column3") := tstrsplit(gene, ":")]
  deg.list %<>% dplyr::rename("tissue" = "column1",
                              "celltype" = "column2",
                              "genename" = "column3") %>%
    as.data.frame()
  
  deg.df <- deg.list %>% 
    dplyr::select(cluster, celltype, genename) %>%
    mutate(type = paste0(celltype, ":", genename))
  
  deg.df.sub <- subset(deg.df, cluster == paste0("M", moduleID))
  plot.df <- deg.df.sub %>%
    dplyr::inner_join(
      deg.df.sub %>%
        dplyr::group_by(genename) %>%
        dplyr::summarise(num = n()) %>%
        dplyr::arrange(-num), 
      by = "genename"
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(num), genename, celltype) %>%
    dplyr::mutate(genename = factor(genename, levels = unique(genename)))
  
  select.gene <- plot.df$genename %>% 
    as.vector() %>% 
    unique() %>% 
    head(n = top.num)
  
  plot.df %<>% subset(genename %in% select.gene)
  
  df_lable <- plot.df %>%
    distinct(genename, .keep_all = TRUE) %>%
    mutate(
      num = num + min(num) / 6,
      rate = 1 / n(),
      cu_rate = cumsum(rate),
      angle = 360 * (cu_rate - rate / 2),
      hjust = if_else(angle <= 180, 0, 1),
      angle = if_else(angle > 180, 270 - angle, 90 - angle)
    )
  
  celltype.num <- deg.df$celltype %>% as.vector() %>% unique() %>% length()
  
  # Set colors
  if (is.null(anno.cols)) {
    rpcols <- colorRampPalette(colors = RColorBrewer::brewer.pal(9, "Paired"))(celltype.num)
  } else {
    rpcols <- anno.cols
  }
  
  # Create the ring plot
  ringplot <- ggplot(plot.df, aes(x = genename, fill = celltype)) +
    geom_bar() +
    coord_polar(theta = 'x') +
    ylim(y_l, NA) +
    scale_fill_manual(values = rpcols) +
    geom_text(
      data = df_lable,
      aes(y = num, label = genename, angle = angle, hjust = hjust), 
      size = text.size, color = "black"
    ) + 
    labs(title = paste0("Module ", moduleID), fill = "Cell Type") +
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank()) +
    theme(plot.title = element_text(size = 8, hjust = 0.5)) +
    theme(legend.title = element_text(size = 8),
          legend.text = element_text(size = 6))
  
  # Return the plot object
  return(ringplot)
}