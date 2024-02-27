# scITDG
`scITDG` is a R package for identifying time-dependent genes at the single-cell resolution
<img src="https://github.com/YandongZheng/scITDG/raw/main/logo.png" alt="Screenshot" style="zoom: 25%;" />

## Requirements
    install.packages(c("devtools", "data.table", "dplyr", "VGAM", "Seurat", "magrittr", "methods", "Biobase", "BiocGenerics", "ggplot2", "Matrix", "ComplexHeatmap", "RColorBrewer", "circlize", "progress", "Hmisc", "unikn", "clusterProfiler", "cowplot", "ggrastr"))


## Install
```R
devtools::install_github("YandongZheng/scITDG")
```

## Quick Start
`scITDG` generally supports the identification and visualization of time-dependent gene expression patterns at the single-cell resolution. 


### 1. Load packages and demo data
The demo data is sourced from the Tabula Muris Senis open-access dataset, which includes mouse limb muscle datasets at 1 month, 3 months, 18 months, 21 months, 24 months, and 30 months. To conserve computational resources, we have sampled only 10,000 single cells for this demonstration. The scanpy object can be download [here](https://figshare.com/ndownloader/files/23873036).


```R
load(Limb_Muscle_10000)

library(Seurat)
library(magrittr)
library(VGAM)
library(stringr)
library(tibble)
library(dplyr)
library(BiocGenerics)
library(Biobase)
library(Matrix)
```

```R
order.list <- c("1m", "3m", "18m", "21m", "24m", "30m")
Limb_Muscle_10000@meta.data$Group <- Limb_Muscle_10000@meta.data$age %>% as.vector
Limb_Muscle_10000@meta.data$Group %<>% factor(levels = order.list)

time.points <- order.list
celltype.list <- seurat.object@meta.data$cell_ontology_class %>% as.vector %>% unique
```



### 2. Seurat was used for  Pairwise differential expression analysis with Seurat

```R
p2pdeg.com <- do.call(rbind, lapply(X = celltype.list, FUN = function(celltype) {
  getP2PDEGs(object = Limb_Muscle_10000, 
             celltype.use = "cell_ontology_class", # celltype anno
             celltype = celltype, 
             time.points.use = "Group",
             time.points = order.list)
})) 
```

```
          p_val avg_log2FC pct.1 pct.2     p_val_adj   gene        cell.type
1 2.508238e-114  2.8103633 0.788 0.002 5.051090e-110   Xist endothelial cell
2  1.754946e-38 -0.9293894 0.159 0.667  3.534111e-34 Col1a2 endothelial cell
3  4.454621e-30 -0.7681975 0.071 0.512  8.970716e-26 Col1a1 endothelial cell
4  7.974473e-27 -0.7233127 0.150 0.562  1.605899e-22 Col3a1 endothelial cell
5  9.754965e-26 -4.1202840 0.049 0.440  1.964455e-21    Dbp endothelial cell
6  2.091209e-24 -2.2712751 0.000 0.356  4.211277e-20  Ddx3y endothelial cell
     class
1 3m.vs.1m
2 3m.vs.1m
3 3m.vs.1m
4 3m.vs.1m
5 3m.vs.1m
6 3m.vs.1m
```

```R
p2p.deg.list <- p2pdeg.com %>% subset(abs(avg_log2FC) >= 0.5) %>% dplyr::select(gene, cell.type, class)
```



### 3. Create the scITDG object and preprocess

```R
sample.ncell.use = 200 # The number of cells sampled

methods::setClass("scITDGDataSet",
                  contains = "ExpressionSet",
                  slots = c(expressionFamily = "vglmff", 
                            lowerDetectionLimit = "numeric",
                            dispFitInfo = "environment"),
                  prototype = prototype(
                    methods::new("VersionedBiobase",
                                 versions = c(classVersion("ExpressionSet"),
                                              scITDGDataSet = "1.2.0" )))
)

new.celltype.list <- NULL
GeneExpCur.wd <- paste0(save.wd, "GeneExpCur/")
dir.create(GeneExpCur.wd)

for (cal.celltype in celltype.list) {
  p2pdeg <- subset(p2p.deg.list, cell.type == cal.celltype) %>%
    .$gene %>% as.vector %>% unique
  ITDGds.raw <- CreateITDGObject(object = Limb_Muscle_10000,
                                 celltype.use = "cell_ontology_class", celltype = cal.celltype,
                                 time.points.use = "Group", time.points = order.list,
                                 random.value = TRUE, seed = 176,
                                 sample.ncell = sample.ncell.use)
  ITDGds.raw <- estimateSizeFactorsForMatrix(ITDGds.raw)
  
  x.inv <- try(
    ITDGds.raw <- estimateDispersionsForMatrix(ITDGds.raw), 
    silent = TRUE)
  if ('try-error' %in% class(x.inv) | length(p2pdeg) < 2) {
    next
  } else {
    ITDGds <- ITDGds.raw[p2pdeg, ]
    ITDGds <- AddTimePointData(ITDGds)
    new.data <- data.frame(Time_num = seq(min(pData(ITDGds)$Time_num),
                                          max(pData(ITDGds)$Time_num),
                                          length.out = ncol(ITDGds)))
    system.time({
      exp.cur <- getGeneExpCur(cds = ITDGds, new_data = new.data, relative_expr = T)
    })
    exp.cur <- ScaleExpCur(ITDGds, exp.cur)
    cal.celltype %<>% gsub(" ", ".", .)
    save(exp.cur, file = paste0(GeneExpCur.wd, cal.celltype, ".exp.cur.rdata"))
    new.celltype.list %<>% c(cal.celltype)
  }
}
```

```R
select.celltype <- new.celltype.list

exp.cur.all <- NULL
for (cal.celltype in select.celltype) {
  load(paste0(GeneExpCur.wd, cal.celltype, ".exp.cur.rdata"))
  rownames(exp.cur) %<>% paste0(tissue, ":", cal.celltype, ":", .)
  exp.cur.all %<>% rbind(exp.cur)
}

mat = as.matrix(exp.cur.all)
row.dist <- as.dist((1 - cor(Matrix::t(mat)))/2)
row.dist[is.na(row.dist)] <- 1
```



### 4. Visualize 

```R
k.num = 6
OrgDb = org.Mm.eg.db
scITDGPlot(exp.cur = exp.cur.all, row.dist = row.dist, k.num = 6, 
           time.points = time.points, 
           sample.ncell.use = sample.ncell.use, 
           show.trajectory = FALSE, show.term = FALSE, OrgDb = org.Mm.eg.db,
           save.wd = save.wd)
```

<img src="https://github.com/YandongZheng/scITDG/raw/main/NoTrajectory_NoTerm_NoOrder.png" alt="Screenshot" style="zoom: 25%;" />



```R
cluster.order.list = c(2, 6, 1, 5, 3, 4)

scITDGPlot(exp.cur = exp.cur.all, row.dist = row.dist, k.num = k.num, 
           time.points = time.points, 
           sample.ncell.use = sample.ncell.use, 
           cluster.order.list = cluster.order.list,
           show.trajectory = FALSE, show.term = FALSE, OrgDb = OrgDb,
           save.wd = save.wd)
```

<img src="https://github.com/YandongZheng/scITDG/raw/main/NoTrajectory_NoTerm.png" alt="Screenshot" style="zoom: 25%;" />

```R
scITDGPlot(exp.cur = exp.cur.all, row.dist = row.dist, k.num = k.num, 
           time.points = time.points, 
           sample.ncell.use = sample.ncell.use, 
           cluster.order.list = cluster.order.list,
           show.trajectory = TRUE, show.term = TRUE, OrgDb = OrgDb,
           save.wd = save.wd)
(-) [====================>------------------------------------------] 33% [ET: 00:00:04]
(\) [===============================>-------------------------------] 50% [ET: 00:00:07]
(|) [=========================================>---------------------] 67% [ET: 00:00:08]
(/) [===================================================>-----------] 83% [ET: 00:00:10]
(-) [===============================================================]100% [ET: 00:00:12]
```

<img src="https://github.com/YandongZheng/scITDG/raw/main/scITDG_Plot.png" alt="Screenshot" style="zoom:25%;" />