#' Title
#'
#' @param clusterExpData
#' @param save.dir
#'
#' @return
#' @export
#'
#' @examples
getClusterData <- function(clusterExpData, save.dir = NULL) {
  gap <- clusterExpData$Tab[, "gap"]
  se <- clusterExpData$Tab[, "SE.sim"]
  decr <- diff(gap) <= 0
  method = "firstSEmax"
  SE.factor = 1
  stopifnot((K <- length(gap)) >= 1, K == length(se), se >= 0, SE.factor >= 0)
  fSE <- SE.factor * se

  k.num <- switch(method, firstmax = {
    decr <- diff(gap) <= 0
    if (any(decr)) which.max(decr) else K
  }, globalmax = {
    which.max(gap)
  }, Tibs2001SEmax = {
    g.s <- gap - fSE
    if (any(mp <- gap[-K] >= g.s[-1])) which.max(mp) else K
  }, firstSEmax = {
    decr <- diff(gap) <= 0
    nc <- if (any(decr)) which.max(decr) else K
    if (any(mp <- gap[seq_len(nc - 1)] >= gap[nc] - fSE[nc])) which(mp)[1] else nc
  }, globalSEmax = {
    nc <- which.max(gap)
    if (any(mp <- gap[seq_len(nc - 1)] >= gap[nc] - fSE[nc])) which(mp)[1] else nc
  })

  df <- as.data.frame(clusterExpData$Tab, stringsAsFactors = TRUE)
  df$clusters <- as.factor(1:nrow(df))
  df$ymin <- gap - se
  df$ymax <- gap + se
  rect.df <- data.frame(xmin = -Inf, xmax = k.num, ymin = -Inf, ymax = df$gap[k.num])
  p <- ggplot2::ggplot(data = df, mapping = aes(x = clusters, y = gap, group = 1)) +
    ggplot2::geom_rect(data = rect.df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = "#b2dfdb", alpha = 1, inherit.aes = FALSE) +
    ggplot2::geom_line(color = "black") +
    ggplot2::geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.2, color = "black") +
    ggplot2::geom_vline(xintercept = k.num, linetype = 2, color = "black") +
    ggplot2::geom_hline(yintercept = df$gap[k.num], linetype = 2, color = "black") +
    ggplot2::geom_point(color = "#009688", size = 1.5) +
    labs(y = "Gap statistic (k)", x = "Number of clusters k", title = "Optimal number of clusters") +
    scale_x_discrete(expand = c(0.02, 0.02)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(axis.text = element_text(size = 8, colour = "black")) +
    theme(axis.title = element_text(size = 8, colour = "black")) +
    theme(plot.title = element_text(size = 10, colour = "black")) +
    theme(aspect.ratio = 0.618)
  if (is.null(save.dir)) {
    ggplot2::ggsave(filename = "./Optimal.num.clusters.pdf", plot = p, width = 4, height = 3)
  } else {
    ggplot2::ggsave(filename = paste0(save.dir, "/Optimal.num.clusters.pdf"), plot = p, width = 4, height = 3)
  }
  list(cluster.df = df,k.num = k.num)
}
