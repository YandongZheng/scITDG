#### title: "Identify time-dependent genes in single-cell transcriptome sequencing data"
#### author: "Yandong Zheng"
#### date: "12/10/2020"


suppressPackageStartupMessages({
  library(Seurat)       #(version 3.2.0)
  library(monocle)      #(version 2.14.0)
  library(htmlwidgets)  #(version 1.3)
  library(ggpubr)       #(version 0.4.0)
  library(pheatmap)     #(version 1.0.12)
  library(colorRamps)   #(version 2.3)
  library(dplyr)        #(version 1.0.2)
  library(RColorBrewer) #(version 1.1-2)
  library(magrittr)     #(version 2.0.1)
})

options(stringsAsFactors=F)
print("--------------All packages are loaded ------------------")

#### object is an S4 object processed by seurat software
object <- readRDS("Preprocessed_object_by_seurat.rds")
Idents(object) <- "Time"
print("--------------object are loaded ------------------")


#### input DEGs list
deg_list <- read.csv("CC_point_to_point_DEGs_list.csv")
time_select_gene <- unique(as.vector(deg_list$gene))

object_matrix <- GetAssayData(object = object)[time_select_gene, ]


#### Make all time points have the same number of cells by downsampling
time <- c("1M","2M","5M","12M","15M") ### time point
CellID <- c()
Sample_ID <- c()
Sample_Size <- 300 ### Set the sample size
for (m in time) {
  sub_sample <- subset(object@meta.data,Time%in%m) %>% rownames(.) %>% as.vector()
  set.seed(176)
  sample_id <- sample(x=sub_sample,size=Sample_Size,replace=T) ### Set the sample size for each sample 
  Sample_ID <- c(Sample_ID,sample_id)
  cellid <- paste0(m,"_",c(1:Sample_Size))
  CellID <- c(CellID,cellid)
}

anno_cell_sample <- data.frame(Sample_ID,CellID)

all_sample_matrix <- c()
all_cell_num <- length(time)*Sample_Size
for (i in 1:all_cell_num) {
  if (i==1) {
    all_sample_matrix <- object_matrix[,which(colnames(object_matrix)%in%Sample_ID[i])] %>% as.matrix()
  }else{
    sample_matrix <- object_matrix[,which(colnames(object_matrix)%in%Sample_ID[i])] %>% as.matrix()
    all_sample_matrix <- cbind(all_sample_matrix,sample_matrix)
  }
}

colnames(all_sample_matrix) <- CellID

meta <- data.frame(Sample_ID=rownames(object@meta.data), Time=object@meta.data$Time)
for (i in 1:all_cell_num) {
  if (i==1) {
    all_sample_meta <- meta[which(meta$Sample_ID%in%Sample_ID[i]),] %>% as.data.frame()
  }else{
    sample_meta <- meta[which(meta$Sample_ID%in%Sample_ID[i]),] %>% as.data.frame()
    all_sample_meta <- rbind(all_sample_meta,sample_meta)
  }
}
all_sample_meta <- merge(all_sample_meta, anno_cell_sample, by="Sample_ID") %>% .[!duplicated(.),]
rownames(all_sample_meta) <- as.vector(all_sample_meta$CellID)
all_sample_meta$CellID <- factor(all_sample_meta$CellID, levels = CellID)
all_sample_meta <- all_sample_meta[order(all_sample_meta$CellID,decreasing = F),]

count_matrix <- all_sample_matrix


### Smoothing single-cell sequencing data

##create pd
pd <- all_sample_meta
##create fd
fd <- as.data.frame(rownames(all_sample_matrix))
colnames(fd) <- "gene_short_name"
rownames(fd) <- fd$gene_short_name
my_cds <- newCellDataSet(count_matrix,
                         phenoData = new("AnnotatedDataFrame", data = pd),
                         featureData = new("AnnotatedDataFrame", data = fd))
DelayedArray:::set_verbose_block_processing(TRUE)

#Accessors and generic functions used in the context of count datasets
my_cds <- estimateSizeFactors(my_cds)
my_cds <- estimateDispersions(my_cds)

time_select_gene <- as.data.frame(time_select_gene)
colnames(time_select_gene) <- "gene"

length.out <- nrow(pData(my_cds))
times.df <- seq(0.1,100,length.out = length.out)
pData(my_cds)$age_time <- times.df


hclust_method = "ward.D2"
num_clusters = 8 #### Set the number of clusters
scale_max = 3
scale_min = -3
cores = 20

num_clusters <- min(num_clusters, nrow(my_cds))
pseudocount <- 1
newdata <- data.frame(age_time = seq(min(pData(my_cds)$age_time),
                                     max(pData(my_cds)$age_time), length.out = all_cell_num))

### Fit smooth spline curves and return the response matrix
smooth_mat <- genSmoothCurves(my_cds, cores = cores, trend_formula = "~sm.ns(age_time, df=3)",
                     relative_expr = T, new_data = newdata)
smooth_mat = smooth_mat[!apply(smooth_mat, 1, sum) == 0, ]
smooth_mat = vstExprs(my_cds, expr_matrix = smooth_mat)
smooth_mat = smooth_mat[!apply(smooth_mat, 1, sd) == 0, ]
smooth_mat = Matrix::t(scale(Matrix::t(smooth_mat), center = TRUE))
smooth_mat = smooth_mat[is.na(row.names(smooth_mat)) == FALSE, ]
smooth_mat[is.nan(smooth_mat)] = 0
smooth_mat[smooth_mat > scale_max] = scale_max
smooth_mat[smooth_mat < scale_min] = scale_min

heatmap_matrix <- smooth_mat
row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
row_dist[is.na(row_dist)] <- 1

bks <- seq(-3.1, 3.1, by = 0.1)
hmcols <- blue2green2red(length(bks) - 1)
ph <- pheatmap(heatmap_matrix, useRaster = TRUE, cluster_cols = FALSE,
               cluster_rows = TRUE, show_rownames = FALSE, show_colnames = FALSE,
               clustering_distance_rows = row_dist, clustering_method = hclust_method,
               cutree_rows = num_clusters, silent = TRUE, filename = NA,
               breaks = bks, border_color = NA, color = hmcols)

annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, num_clusters)))

feature_label <- as.character(fData(my_cds)[row.names(heatmap_matrix),"gene_short_name"])
feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
row.names(heatmap_matrix) <- feature_label
colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))

annotation_col <- data.frame(Time = rep(time,each=Sample_Size))
annotation_col$Time <- factor(annotation_col$Time,levels = time)
row.names(annotation_col) <- 1:all_cell_num

cluster_cols <- colorRampPalette(brewer.pal(8,"Dark2"))(num_clusters)
names(cluster_cols) <- 1:num_clusters

ann_colors = list(Time = c("1M"="red",
                           "2M"="#FF8C00",
                           "5M"="#8A2BE2",
                           "12M"="#00A587",
                           "15M"="#374F8B"),
                  Cluster = cluster_cols)

ph_res <- pheatmap(heatmap_matrix[, ], useRaster = TRUE, cluster_cols = FALSE,
                   cluster_rows = TRUE, show_rownames = FALSE,
                   show_colnames = FALSE, clustering_distance_rows = row_dist,
                   clustering_method = hclust_method, cutree_rows = num_clusters,
                   annotation_row = annotation_row, annotation_col = annotation_col,
                   annotation_colors = ann_colors,annotation_names_row = FALSE,annotation_names_col = FALSE,
                   treeheight_row = 0, breaks = bks, fontsize = 8, color = hmcols,
                   border_color = NA, silent = TRUE, filename = NA)
pdf("Time_dependent_DEGs_heatmap_v1.pdf",height = 7,width = 5)
ph_res
dev.off()

### Get clustering for each gene
clusters <- cutree(ph_res$tree_row, k = num_clusters)
clustering <- data.frame(clusters)
clustering[,1] <- as.numeric(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
clustering$Gene <- rownames(clustering)
clustering <- clustering[order(clustering$Gene_Clusters,decreasing = F),]

write.csv(clustering,"Time_dependent_DEGs_clustering.csv",row.names = F)