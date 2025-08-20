#' Title
#'
#' @param object
#' @param k.num
#' @param cluster.order.list
#' @param integrate.cluster.list
#' @param show.trajectory
#' @param show.term
#' @param OrgDb
#' @param save.wd
#' @param use_raster
#' @param plot.height
#' @param heatmap.cols
#' @param scale.adj
#'
#' @return
#' @export
#'
#' @examples
scITDGIntegratePlot <- function(object,
                       k.num = 6,
                       cluster.order.list = NULL,
                       integrate.cluster.list = NULL,
                       show.trajectory = TRUE,
                       show.term = TRUE,
                       OrgDb = org.Hs.eg.db,
                       save.wd = NULL,
                       use_raster = TRUE,
                       plot.height = NULL,
                       heatmap.cols = NULL,
                       scale.adj = NULL) {
  if (is.null(save.wd)) {
    save.wd <- "./"
  }
  
  exp.cur = object$exp.cur
  object.gene = object$object.gene
  row.dist = object$row.dist
  time.points = object$time.points
  sample.ncell.use = object$sample.ncell.use
  
  
  mat = as.matrix(exp.cur)
  tree.row = hclust(row.dist, method = "ward.D2")
  mat = mat[tree.row$order, , drop = FALSE]
  labels.row = rownames(mat)
  gene.cluster = cutree(tree.row, k.num)[tree.row$order]
  
  order.mat.df <- data.frame(
    cluster = paste0("M", as.vector(gene.cluster)),
    gene = labels.row,
    id = NA) %>%
    group_by(cluster) %>%
    mutate(id = seq(length(id))) %>%
    ungroup %>%
    as.data.frame
  
  save.type <- ""
  
  if (!is.null(integrate.cluster.list)) {
    order.mat.df.integrate <- NULL
    for (i in 1:length(integrate.cluster.list)) {
      sel.cluster <- integrate.cluster.list[[i]]
      sel.gene <- subset(order.mat.df, cluster %in% paste0("M", sel.cluster)) %>% .$gene
      mat.sub <- mat[sel.gene, ] %>% as.matrix
      row.dist.sub <- as.dist((1 - cor(Matrix::t(mat.sub)))/2)
      row.dist.sub[is.na(row.dist.sub)] <- 1
      tree.row.sub = hclust(row.dist.sub, method = "ward.D2")
      mat.sub = mat.sub[tree.row.sub$order, , drop = FALSE]
      labels.row.sub = rownames(mat.sub)
      order.mat.df.sub <- data.frame(
        cluster = paste0("M", i),
        gene = labels.row.sub) %>%
        mutate(id = seq(length(labels.row.sub))) %>%
        ungroup %>%
        as.data.frame
      order.mat.df.integrate %<>% rbind(order.mat.df.sub)
    }
    
    k.num <- length(integrate.cluster.list)
    order.mat.df.integrate$cluster %<>% as.vector %>% factor(levels = paste0("M", seq_len(k.num)))
    order.mat.df <- order.mat.df.integrate
    save.type <- "Integrate"
  }
  
  if (is.null(plot.height)) {
    plot.height = 5 + (k.num - 4) * 0.75
  }
  
  if (is.null(scale.adj)) {
    if (k.num == 2) { scale.adj = 0.2 }
    if (k.num == 3) { scale.adj = 0.12 }
    if (k.num == 4) { scale.adj = 0.091 }
    if (k.num == 5) { scale.adj = 0.072 }
    if (k.num == 6) { scale.adj = 0.06 }
    if (k.num == 7) { scale.adj = 0.05 }
    if (k.num == 8) { scale.adj = 0.04 }
    if (k.num == 9) { scale.adj = 0.03 }
    if (k.num == 10) { scale.adj = 0.01 }
  } else {
    scale.adj = 0.06
  }
  
  if (!is.null(cluster.order.list)) {
    order.mat.df$cluster %<>% as.vector %>% factor(levels = paste0("M", cluster.order.list))
    order.mat.df %<>% arrange(cluster, id)
  }else{
    order.mat.df$cluster %<>% as.vector %>% factor(levels = paste0("M", seq_len(k.num)))
    order.mat.df %<>% arrange(cluster, id)
  }
  
  mat %<>% .[order.mat.df$gene, ]
  
  # heatmap color
  bks <- seq(-3.1, 3.1, by = 0.1)
  if (is.null(heatmap.cols)) {
    hmcols <- colorRampPalette(
      colors = RColorBrewer::brewer.pal(11, "Spectral"))(length(bks) - 1) %>%
      rev()
  } else {
    hmcols = colorRampPalette(colors = heatmap.cols)(length(bks) - 1)
  }
  
  # time point color
  time.cols = usecol(pal_seeblau, n = length(time.points))
  names(time.cols) <- time.points
  
  # module color
  cluster.cols <- colorRampPalette(brewer.pal(8,"Dark2"))(k.num)
  names(cluster.cols) <- paste0("M",seq_len(k.num))
  
  
  # row_split = annotation_row$cluster
  module.anno = rowAnnotation(foo = anno_block(
    gp = gpar(fill = cluster.cols[levels(order.mat.df$cluster)], col = NA), # 这里注意填充顺序
    labels = levels(order.mat.df$cluster),
    labels_gp = gpar(col = "black", fontsize = 8),
    width = unit(3, "mm"))) # 设置注释调宽度)
  
  # 右侧延长注释线条
  anno.data.mark <- table(order.mat.df$cluster) %>% as.data.frame %>%
    dplyr::rename(cluster = "Var1", num = "Freq")
  for (i in seq_len(nrow(anno.data.mark))) {
    num.median <- ceiling(anno.data.mark$num[i] / 2)
    if (i == 1) {
      num.median.all <- num.median
    } else {
      num.median <- num.median + cumsum(anno.data.mark$num[1:(i - 1)]) %>% max
      num.median.all %<>% c(num.median)
    }
  }
  
  mark.anno = rowAnnotation(
    link = anno_mark(
      at = num.median.all,
      link_width = unit(3.5, "mm"), # 控制长短
      #link_height = unit(1, "mm"),
      labels = levels(order.mat.df$cluster),
      which = "row",
      link_gp = gpar(col = cluster.cols[levels(order.mat.df$cluster)], lwd = 2.8), # lwd 控制粗细
      labels_gp = gpar(fontsize = 0)))
  
  # 绘制热图
  suppressMessages(
    heatmap.plot <- ComplexHeatmap::Heatmap(mat, #name = "Expression",
                                            col = hmcols, #定义热图由低值到高值的渐变颜色 viridis(length(bks))
                                            show_row_names = FALSE,  #不展示基因名称
                                            show_column_names = FALSE,  #不展示基因名称
                                            cluster_columns = FALSE,
                                            cluster_rows = FALSE,
                                            column_title = NULL,
                                            row_title = NULL,
                                            left_annotation = module.anno,
                                            #right_annotation = mark.anno,
                                            row_split = order.mat.df$cluster,
                                            column_names_gp = gpar(fontsize = 6),
                                            row_names_gp = gpar(fontsize = 6),
                                            row_gap = unit(2, "mm"),
                                            show_heatmap_legend = F)
  )
  
  lgd.top = Legend(title = "Time Point", labels = time.points,
                   legend_gp = gpar(fill = time.cols),
                   labels_gp = gpar(fontsize = 8),
                   title_gp = gpar(fontsize = 8, fontface = "bold"))
  
  lgd.row = Legend(title = "Module", labels = paste0("M",seq_len(k.num)),
                   legend_gp = gpar(fill = cluster.cols),
                   labels_gp = gpar(fontsize = 8),
                   title_gp = gpar(fontsize = 8, fontface = "bold"))
  
  col_fun_prop = circlize::colorRamp2(breaks = seq(-3.1, 3.1, by = 6.2/61), colors = hmcols)
  lgd.exp = Legend(col_fun = col_fun_prop, title = "Rel.Exp", at = c(-3, 0, 3),
                   labels = c("Low", "Median", "High"), #break_dist = 1,
                   labels_gp = gpar(fontsize = 8),
                   title_gp = gpar(fontsize = 8, fontface = "bold"))
  
  lgd = packLegend(list = list(lgd.top, lgd.row, lgd.exp))
  
  # bar图theme
  bar.theme <- function() {
    theme(panel.grid = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_text(size = 8, colour = NA),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none")
  }
  
  anno.data.col <- data.frame(num = 0.2, age.group = rep(time.points, each = sample.ncell.use),
                              cell = seq_len(length(time.points) * sample.ncell.use))
  anno.data.col$age.group %<>% factor(., levels = time.points)
  
  anno.bar.col <- ggplot() +
    geom_col(anno.data.col, mapping = aes(x = cell, y = num, fill = age.group),
             show.legend = F, width = 1) +
    scale_fill_manual(values = time.cols) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    bar.theme() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))
  
  anno.data.row <- data.frame(num = 0.2, Expression = bks)
  anno.bar.row <- ggplot() +
    geom_col(anno.data.row, mapping = aes(x = Expression, y = num, fill = Expression),
             show.legend = F, width = 0.1) + ##注意，宽度不能设置太大，否则末端不合理
    scale_fill_gradientn(colours = hmcols) +
    scale_x_continuous(expand = c(0.02, 0.02)) +
    coord_flip() +
    bar.theme() +
    theme(plot.margin = unit(c(0,0,0,0), "cm"))
  
  
  # 是否显示轨迹
  if (show.trajectory) {
    order.cluster <- order.mat.df$cluster %>% levels
    colnames(mat) <- seq_len(ncol(mat))
    
    pb <- progress_bar$new(format = "(:spin) [:bar]:percent [ET: :elapsedfull || ETA: :eta]",
                           total = length(order.cluster),
                           complete = "=",   # Completion bar character
                           incomplete = "-", # Incomplete bar character
                           current = ">",    # Current bar character
                           clear = TRUE,    # If TRUE, clears the bar when finish
                           width = 100)      # Width of the progress bar
    # 绘制平滑轨迹图
    plot.smooth.list <- list()
    for (i in 1:length(order.cluster)) {
      # 显示进度条
      pb$tick()
      
      newgene <- subset(order.mat.df, cluster == order.cluster[i]) %>% .$gene %>% as.vector
      newdata <- mat[newgene,]
      tmp_data <- t(newdata)
      
      myFUN <- function(x){
        df = data.frame(count = x, cell = rownames(tmp_data))
        return(df)
      }
      
      line_data = do.call(rbind, apply(tmp_data,2,myFUN))
      line_data$group = rep(colnames(tmp_data), each = length(rownames(tmp_data)))
      line_data$cell = factor(line_data$cell, levels = seq_len(length(time.points) * sample.ncell.use))
      line_data$smooth = rep('smooth', length(rownames(line_data)))
      
      line_data$cell = as.numeric(line_data$cell)
      line_data$age.group <- rep(rep(time.points, each = sample.ncell.use), ncol(tmp_data))
      
      gene_num <- length(newgene)
      
      g <- ggplot(data = line_data,aes(x = cell,y = count, group = group))
      title <- paste0(order.cluster[i], " (n = ", gene_num, ")") %>% gsub("M", "Module ", .)
      x.limits = c(1, length(time.points) * sample.ncell.use)
      x.breaks = seq((length(time.points) / 2), length(time.points) * sample.ncell.use, sample.ncell.use)
      p.line <- g +
        # rasterise(
        #   stat_smooth(method = lm, formula = y~poly(x,3), se = F, colour = "#DCDCDC", size = 0.5, alpha = 0.2),
        #   dpi = 300) +
        geom_smooth(aes(group = smooth), colour = cluster.cols[order.cluster[i]], formula = y ~ poly(x, 3),
                    method = lm, size = 1.5, alpha = 1, se = F) +
        labs(x = "", y = "") +
        annotate("text", x = (length(time.points) * sample.ncell.use)/2,
                 y = 2.6, label = title, size = 2.5, colour = "black") +
        coord_cartesian(ylim = c(-3, 3)) +
        scale_y_continuous(breaks = seq(-3, 3, 1), expand = c(0.02, 0.02)) +
        scale_x_continuous(limits = x.limits, breaks = x.breaks, expand = c(0.02, 0.02),
                           labels = time.points) +
        theme(panel.grid = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank(),
              plot.title = element_text(size = 8, colour = "black", hjust = 0.5),
              axis.text.x = element_text(size = 8, colour = NA),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_text(size = 8, colour = "black")) +
        theme(panel.border = element_rect(fill = NA, color = "black", size = 0.75, linetype = "solid")) +
        theme(plot.margin = unit(c(0,0,0,0), "cm")) # 左边和上边页边距为0
      
      # 创建一个等边三角形
      anno.data.triangle <- data.frame(
        x = c(0, 0, sqrt(3) / 2),
        y = c(0, 1, 0.5)
      ) %>% .[c(2, 3, 1), ] # 将三角形逆时针旋转90度
      anno.data.triangle$x %<>% -.
      p.tri <- ggplot(anno.data.triangle, aes(x = x, y = y)) +
        geom_polygon(fill = cluster.cols[order.cluster[i]]) +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_continuous(expand = c(0, 0)) +
        coord_fixed() +
        bar.theme() +
        theme(plot.margin = unit(c(0,0,0,0), "cm"))
      
      p.final <- cowplot::plot_grid(p.tri, anno.bar.row, NULL, p.line, nrow = 1,
                                    rel_widths = c(2, 1, -1.2, 12), align = "h")
      plot.smooth.list[[i]] <- p.final
    }
    plot.smooth <- cowplot::plot_grid(plotlist = plot.smooth.list, ncol = 1, align = "hv")
    
    
    # 是否显示term，绘制term图
    if (show.term) {
      get_wraper <- function(width) {
        function(x) {
          lapply(strwrap(x, width = width, simplify = FALSE), paste, collapse = "\n")
        }
      }
      nchar.cut <- 35
      
      pb <- progress_bar$new(format = "(:spin) [:bar]:percent [ET: :elapsedfull || ETA: :eta]",
                             total = length(order.cluster),
                             complete = "=",   # Completion bar character
                             incomplete = "-", # Incomplete bar character
                             current = ">",    # Current bar character
                             clear = TRUE,    # If TRUE, clears the bar when finish
                             width = 100)      # Width of the progress bar
      
      plot.term.list <- list()
      for (i in 1:length(order.cluster)) {
        pb$tick()
        
        newgene <- subset(order.mat.df, cluster == order.cluster[i]) %>% .$gene %>% as.vector
        gene.go <- as.vector(newgene)
        
        gene.go.all <- NULL
        for (g in gene.go) {
          gene.select <- strsplit(g, ":") %>% unlist %>% tail(., n = 1)
          gene.select %<>% strsplit(., "-") %>% unlist %>% head(., n = 1)
          gene.go.all %<>% c(gene.select)
        }
        
        suppressMessages(
          gene.df <- bitr(gene.go.all, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = OrgDb)
        )
        
        universe.df.all <- NULL
        for (g in object.gene) {
          gene.select <- strsplit(g, "-") %>% unlist %>% head(., n = 1)
          universe.df.all %<>% c(gene.select)
        }
        
        suppressMessages(
          universe.df <- bitr(universe.df.all, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = OrgDb)
        )
        
        ego <- enrichGO(gene = gene.df$ENTREZID,
                        universe = universe.df$ENTREZID, #背景基因集
                        OrgDb = OrgDb, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                        #keyType = 'SYMBOL', # ENSEMBL
                        ont = "BP", #也可以是 CC  BP  MF中的一种
                        pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                        pvalueCutoff = 0.05, #P值会过滤掉很多，可以全部输出
                        qvalueCutoff = 0.2,
                        readable = TRUE) #Gene ID 转成gene Symbol ，易读
        
        term.df <- fortify(ego, showCategory = 3) %>%
          dplyr::select(Description, p.adjust) %>%
          dplyr::mutate(neg_LogP = -log10(p.adjust)) %>%
          dplyr::arrange(desc(neg_LogP))
        term.df$neg_LogP[which(term.df$neg_LogP > 50)] <- 50
        term.df$Description %<>% as.vector
        
        nchar.term <- nchar(term.df$Description)
        if (!max(nchar.term) > nchar.cut) {
          term.df$Description[1] <- paste0(term.df$Description[1], paste(rep(" ", (nchar.cut - nchar.term[1])), collapse = ""))
        }
        term.df$Description %<>% capitalize %>% factor(., levels = rev(.))
        
        p.term <- ggplot(data = term.df, mapping = aes(x = neg_LogP, y = Description)) +
          geom_bar(stat = "identity", width = 0.5, fill = cluster.cols[order.cluster[i]]) +
          labs(x = "", y = "") +
          theme(panel.grid = element_blank(),
                panel.background = element_blank(),
                plot.background = element_blank(),
                axis.title = element_text(size = 8, colour = "black"),
                axis.text = element_text(size = 8, colour = "black")) +
          theme(panel.border = element_rect(fill = NA, color = "black", size = 0.75, linetype = "solid")) +
          theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) + #  上 右 下 左
          scale_x_reverse() +
          scale_y_discrete(labels = get_wraper(nchar.cut), position = "right")
        plot.term.list[[i]] <- p.term
      }
      
      plot.term <- cowplot::plot_grid(plotlist = plot.term.list, ncol = 1, align = "hv")
      plot.ggplot <- cowplot::plot_grid(plot.smooth, NULL, plot.term, nrow = 1,
                                        rel_widths = c(1.65, 0.05, 3), align = "hv")
      
      # -log10注释文字
      plot.term.text <- ggplot() +
        annotate("text", x = 0, y = 0, label = "-log10 (p.adjust)", size = 2.8, colour = "black") +
        theme(panel.grid = element_blank(),
              panel.background = element_blank(),
              plot.background = element_blank(),
              axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank()) +
        theme(plot.margin = unit(c(0,0,0,0), "cm")) # 左边和上边页边距为0
      
      
      # 绘制对应线
      all.label = rownames(mat)
      mark.label = rownames(mat)[num.median.all]
      
      gap.height <- round(2 / ((plot.height * 25.4 - 10 - 2*(k.num-1))/nrow(mat)))
      
      cluster.cumsum <- cumsum(anno.data.mark$num)
      cluster.num <- anno.data.mark$num
      all.label.new <- all.label
      for (i in 1:(k.num-1)) {
        int.gap <- paste0("gap:", i, ":", 1:gap.height)
        all.label.new <- append(x = all.label.new, values = int.gap, after = cluster.cumsum[i])
        cluster.cumsum[i+1] <- cluster.cumsum[i + 1] + gap.height*i
        cluster.num[i+1] <- cluster.num[i + 1] + gap.height
      }
      all.label.new <- c(paste0("gap:top:", 1:gap.height),
                         all.label.new,
                         paste0("gap:bottom:", 1:c(gap.height * 4)))
      all.label.new %<>% rev
      
      top.scale <- 1 / (plot.height * 25.4) * 2 #热图顶部2mm
      bottom.scale <- 1 / (plot.height * 25.4) * 8 #热图顶部8mm
      mark.scale = c(0 + bottom.scale + scale.adj, 1 - top.scale - scale.adj)
      
      idx <- match(mark.label, all.label.new)
      pos.raw <- scales::rescale(idx, to = c(0, 1), from = c(0.5, length(all.label.new) + 0.5))
      pos.pos <- scales::rescale(1:length(idx), to = mark.scale)
      pos.idx <- match(idx, sort(idx))
      pos.new <- pos.pos[pos.idx]
      
      seg.x0 = 0.75
      link.line.length = 0.035
      seg.y0 <- pos.raw
      seg.x1 <- seg.x0 + link.line.length
      seg.y1 <- pos.raw
      seg.x2 <- seg.x1 + link.line.length + 0.01
      seg.y2 <- pos.new
      
      # viewport
      link.line.gp = gpar(col = cluster.cols[levels(order.mat.df$cluster)], lwd = 2.8)
      raw.annoline.subvp <- grid::viewport(x = 0.225, y = 0.5, width = 0.35, height = 1) # x越小越靠左边，y越小越靠下
      
      link.1.grob <- segmentsGrob(x0 = seg.x0, x1 = seg.x1,
                                  y0 = seg.y0, y1 = seg.y0,
                                  vp = raw.annoline.subvp,
                                  gp = link.line.gp)
      link.2.grob <- segmentsGrob(x0 = seg.x1, x1 = seg.x1,
                                  y0 = seg.y1, y1 = seg.y2,
                                  vp = raw.annoline.subvp,
                                  gp = link.line.gp)
      link.3.grob <- segmentsGrob(x0 = seg.x1, x1 = seg.x2,
                                  y0 = seg.y2, y1 = seg.y2,
                                  vp = raw.annoline.subvp,
                                  gp = link.line.gp)
      
      pdf(paste0(save.wd, tissue, "_scITDGPlot_", save.type,".pdf"),
          width = 8, height = plot.height, bg = "transparent")
      draw(heatmap.plot, padding = unit(c(8, 2, 2, 140), "mm")) ## see right heatmap in following 底部、左侧、顶部和右
      
      grid.draw(grid::gList(link.1.grob, link.2.grob, link.3.grob))
      
      line.subvp <- grid::viewport(x = 0.632, y = 0.493, width = 0.6, height = 0.985) # x越小越靠左边，y越小越靠下
      print(plot.ggplot, vp = line.subvp)
      
      col.bar.subvp <- grid::viewport(x = 0.1675, y = 0.018, width = 0.2975, height = 0.06)
      print(anno.bar.col, vp = col.bar.subvp)
      
      col.bar.subvp <- grid::viewport(x = 0.460, y = 0.018, width = 0.170, height = 0.06)
      print(anno.bar.col, vp = col.bar.subvp)
      
      bar.text.subvp <- grid::viewport(x = 0.6, y = 0.015, width = 0.2, height = 0.06)
      print(plot.term.text, vp = bar.text.subvp)
      
      draw(lgd, x = unit(0.99, "npc"), y = unit(0.5, "npc"), just = c("right")) # 1 npc为最右边
      dev.off()
      
    } else {
      plot.ggplot <- plot.smooth
      
      # 绘制对应线
      all.label = rownames(mat)
      mark.label = rownames(mat)[num.median.all]
      
      gap.height <- round(2 / ((plot.height * 25.4 - 10 - 2*(k.num-1))/nrow(mat)))
      
      cluster.cumsum <- cumsum(anno.data.mark$num)
      cluster.num <- anno.data.mark$num
      all.label.new <- all.label
      for (i in 1:(k.num-1)) {
        int.gap <- paste0("gap:", i, ":", 1:gap.height)
        all.label.new <- append(x = all.label.new, values = int.gap, after = cluster.cumsum[i])
        cluster.cumsum[i+1] <- cluster.cumsum[i + 1] + gap.height*i
        cluster.num[i+1] <- cluster.num[i + 1] + gap.height
      }
      all.label.new <- c(paste0("gap:top:", 1:gap.height),
                         all.label.new,
                         paste0("gap:bottom:", 1:c(gap.height * 4)))
      all.label.new %<>% rev
      
      top.scale <- 1 / (plot.height * 25.4) * 2 #热图顶部2mm
      bottom.scale <- 1 / (plot.height * 25.4) * 8 #热图顶部8mm
      mark.scale = c(0 + bottom.scale + scale.adj, 1 - top.scale - scale.adj)
      
      idx <- match(mark.label, all.label.new)
      pos.raw <- scales::rescale(idx, to = c(0, 1), from = c(0.5, length(all.label.new) + 0.5))
      pos.pos <- scales::rescale(1:length(idx), to = mark.scale)
      pos.idx <- match(idx, sort(idx))
      pos.new <- pos.pos[pos.idx]
      
      seg.x0 = 0.75
      link.line.length = 0.035
      seg.y0 <- pos.raw
      seg.x1 <- seg.x0 + link.line.length
      seg.y1 <- pos.raw
      seg.x2 <- seg.x1 + link.line.length + 0.03
      seg.y2 <- pos.new
      
      # viewport
      link.line.gp = gpar(col = cluster.cols[levels(order.mat.df$cluster)], lwd = 2.8)
      raw.annoline.subvp <- grid::viewport(x = 0.365, y = 0.5, width = 0.35, height = 1) # x越小越靠左边，y越小越靠下
      
      link.1.grob <- segmentsGrob(x0 = seg.x0, x1 = seg.x1,
                                  y0 = seg.y0, y1 = seg.y0,
                                  vp = raw.annoline.subvp,
                                  gp = link.line.gp)
      link.2.grob <- segmentsGrob(x0 = seg.x1, x1 = seg.x1,
                                  y0 = seg.y1, y1 = seg.y2,
                                  vp = raw.annoline.subvp,
                                  gp = link.line.gp)
      link.3.grob <- segmentsGrob(x0 = seg.x1, x1 = seg.x2,
                                  y0 = seg.y2, y1 = seg.y2,
                                  vp = raw.annoline.subvp,
                                  gp = link.line.gp)
      
      pdf(paste0(save.wd, tissue, "_scITDGPlot_NoTerm_", save.type, ".pdf"),
          width = 5, height = plot.height, bg = "transparent")
      draw(heatmap.plot, padding = unit(c(8, 2, 2, 70), "mm")) ## see right heatmap in following 底部、左侧、顶部和右
      
      grid.draw(grid::gList(link.1.grob, link.2.grob, link.3.grob))
      
      line.subvp <- grid::viewport(x = 0.65, y = 0.493, width = 0.35, height = 0.985) # x越小越靠左边，y越小越靠下
      print(plot.ggplot, vp = line.subvp)
      
      col.bar.subvp <- grid::viewport(x = 0.244, y = 0.018, width = 0.425, height = 0.06)
      print(anno.bar.col, vp = col.bar.subvp)
      
      col.bar.subvp <- grid::viewport(x = 0.686, y = 0.018, width = 0.282, height = 0.06)
      print(anno.bar.col, vp = col.bar.subvp)
      
      draw(lgd, x = unit(0.99, "npc"), y = unit(0.5, "npc"), just = c("right")) # 1 npc为最右边
      dev.off()
      
    }
    
  } else {
    if (!is.null(cluster.order.list)) {
      only_heatmap_name <- paste0(save.wd, tissue, "_scITDGPlot_NoTrajectory_NoTerm", save.type, ".pdf")
    } else {
      only_heatmap_name <- paste0(save.wd, tissue, "_scITDGPlot_NoTrajectory_NoTerm_NoOrder", save.type, ".pdf")
    }
    
    pdf(only_heatmap_name, width = 4, height = 5, bg = "transparent")
    draw(heatmap.plot, padding = unit(c(8, 2, 2, 20), "mm")) ## see right heatmap in following 底部、左侧、顶部和右
    col.bar.subvp <- grid::viewport(x = 0.425, y = 0.018, width = 0.785, height = 0.06)
    print(anno.bar.col, vp = col.bar.subvp)
    draw(lgd, x = unit(0.99, "npc"), y = unit(0.5, "npc"), just = c("right")) # 1 npc为最右边
    dev.off()
  }
  saveRDS(order.mat.df, file = paste0(save.wd, tissue, "_scITDG_Cluster_GeneName", save.type, ".rds"))
}

