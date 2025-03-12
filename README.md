# scITDG

## A tool for identifying time-dependent genes in single-cell transcriptome sequencing data 
scITDG, a tool designed for the analysis of time-dependent gene expression in single-cell transcriptomic sequencing data. scITDG stands out by its unique ability to identify dynamic gene expression patterns across various time points at single-cell resolution, which is pivotal for deciphering complex biological processes such as aging and tissue regeneration. The tool is fully compatible with widely used single-cell analysis platforms like Seurat and Scanpy. By integrating natural cubic splines regression with bootstrapping resampling, scITDG not only enhances the functionality of these platforms but also broadens their applicability.

<p align="center">
  <img src="https://github.com/YandongZheng/scITDG/raw/main/images/scITDG_framework_overview.png" alt="Description" width=20%;" />
</p>



Furthermore, the versatility of scITDG makes it applicable to a wide range of biological contexts, such as:
<ul>
    <li>Aging
    <li>Regeneration
    <li>Development
    <li>Disease Progression
    <li>Therapeutic Responses
</ul>  

In addition, the [Seurat](https://satijalab.org/seurat/), or [Scanpy](https://scanpy.readthedocs.io/) packages must be installed for scITDG to take Seurat / scanpy objects as input, respectively.

scITDG has been tested with R version 4.0.0 and higher.

## Installation

To install scITDG, first install the devtools and Seurat package, if it is not already installed: 

```R
install.packages("devtools")
install.packages("Seurat")
```

Finally, install scITDG from GitHub: 
```R
devtools::install_github("YandongZheng/scITDG")
```

## Quick Start
`scITDG` generally supports the identification and visualization of time-dependent gene expression patterns at the single-cell resolution. 


### 1. Load Packages and Demo Data
The demo data (Limb_Muscle_10000) is sourced from the Tabula Muris Senis open-access dataset, which includes mouse limb muscle datasets at 1 month, 3 months, 18 months, 21 months, 24 months, and 30 months. To conserve computational resources, we have sampled only 10,000 single cells for this demonstration. The scanpy object can be download [here](https://figshare.com/ndownloader/files/23873036).

#### Load Packages
```R
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
#### Load Demo Data

##### `Limb_Muscle_10000` has already been converted into a Seurat object.
```R
load(Limb_Muscle_10000)
```



#### Preparation of Time Points and Cell Type List

`time.points` is a vector of time points sorted in ascending order. It must correspond to the time sequence in the `meta.data` of the Seurat object. For example, in the `Limb_Muscle_10000` object, the time sequence is sorted as "1m", "3m", "18m", "21m", "24m", and "30m".
The time information is stored in the `age` column of the `meta.data`. To facilitate subsequent analyses, we redirect the `age` column to a new column named `Group` and convert it into a factor with levels corresponding to the `time.points`.
`celltype.list` is a list of all unique cell types present in the Seurat object, extracted from the `cell_ontology_class` column in the `meta.data`.

```R
# Define the time points in the correct order
time.points <- c("1m", "3m", "18m", "21m", "24m", "30m")

# Redirect the 'age' column to 'Group' and convert it to a factor with specified levels
Limb_Muscle_10000@meta.data$Group <- Limb_Muscle_10000@meta.data$age %>% as.vector
Limb_Muscle_10000@meta.data$Group %<>% factor(levels = time.points)

# Extract the unique list of cell types from the Seurat object's meta.data
celltype.list <- Limb_Muscle_10000@meta.data$cell_ontology_class %>% as.vector %>% unique
```



### 2. Seurat was used for  Pairwise differential expression analysis with Seurat

#### Pairwise Differential Expression Analysis with Seurat

We perform pairwise differential expression analysis using the Seurat package. The function `getP2PDEGs` is applied to each cell type in `celltype.list`. This function identifies differentially expressed genes between time points specified in `time.points`. The results are combined into a single data frame (`p2pdeg.com`) for further analysis.

Explanation of Parameters:
`celltype.use`: Specifies the column name in the Seurat object's meta.data that contains the cell type annotations.
`celltype`: The specific cell type being processed in the current iteration. This is used to subset the analysis to the desired cell type.
`time.points.use`: Specifies the column name in the Seurat object's meta.data that contains the time information.
`time.points`: A vector of time points in ascending order, used to define the sequence of pairwise comparisons.

```R
p2pdeg.com <- do.call(rbind, lapply(X = celltype.list, FUN = function(celltype) {
  getP2PDEGs(object = Limb_Muscle_10000, 
             celltype.use = "cell_ontology_class",  # Column name in Seurat object's meta.data that annotates cell types
             celltype = celltype,                   # Specific cell type to be processed in this iteration
             time.points.use = "Group",             # Column name in Seurat object's meta.data that annotates time information
             time.points = time.points)             # Vector of time points in ascending order
}))
```

#### Example Output of Differential Expression Analysis

The resulting data frame (`p2pdeg.com`) contains differentially expressed genes along with their associated statistics and cell types. The columns include p-values, log2 fold changes, percentage of cells expressing the gene in each group, adjusted p-values, gene names, cell types, and the comparison class.

```
          p_val avg_log2FC pct.1 pct.2     p_val_adj   gene        cell.type     class
1 2.508238e-114  2.8103633 0.788 0.002 5.051090e-110   Xist endothelial cell  3m.vs.1m
2  1.754946e-38 -0.9293894 0.159 0.667  3.534111e-34 Col1a2 endothelial cell  3m.vs.1m
3  4.454621e-30 -0.7681975 0.071 0.512  8.970716e-26 Col1a1 endothelial cell  3m.vs.1m
4  7.974473e-27 -0.7233127 0.150 0.562  1.605899e-22 Col3a1 endothelial cell  3m.vs.1m
5  9.754965e-26 -4.1202840 0.049 0.440  1.964455e-21    Dbp endothelial cell  3m.vs.1m
6  2.091209e-24 -2.2712751 0.000 0.356  4.211277e-20  Ddx3y endothelial cell  3m.vs.1m
```

#### Saving the Differential Expression Results

The resulting list of differentially expressed genes can be saved as an RDS file for future reference or further analysis.
```R
saveRDS(object = p2pdeg.com, file = paste0(save.wd, "p2pdeg.com.rds"))
```



### 3. Create the scITDG Object and Perform the Main Calculation

Explanation of Parameters:
`object`: The Seurat object containing the single-cell data (Limb_Muscle_10000 in this case).
`celltype.list`: A list of unique cell types to be processed, extracted from the Seurat object's meta.data.
`logfc.threshold`: The log2 fold change threshold for identifying differentially expressed genes (Default is 0.5).
`pvalue.threshold`: The p-value threshold for identifying differentially expressed genes (Default is 0.05).
`p2p.deg.use`: The precomputed pairwise differentially expressed genes from the previous analysis (p2pdeg.com).
`celltype.use`: The column name in the Seurat object's meta.data that contains cell type annotations (cell_ontology_class).
`time.points.use`: The column name in the Seurat object's meta.data that contains time point annotations (Group).
`time.points`: A vector of time points in ascending order, used to define the sequence of analyses.
`sample.ncell.use`: The number of cells to sample for each cell type (Default is 200).
`save.wd`: The working directory where results will be saved.
`n.cores`: The number of CPU cores to use for parallel processing (Default is 1).


```R
scitdg <- CalculatescITDG(
  object = Limb_Muscle_10000, 
  p2p.deg.use = p2pdeg.com, 
  celltype.list = celltype.list, 
  logfc.threshold = 0.5, 
  pvalue.threshold = 0.05, 
  celltype.use = "cell_ontology_class",  # Column name in Seurat object's meta.data for cell type annotation
  time.points.use = "Group",              # Column name in Seurat object's meta.data for time point annotation
  time.points = time.points,              # Vector of time points in ascending order
  sample.ncell.use = 200,                 # Number of cells to sample for each cell type
  save.wd = save.wd,                      # Working directory for saving results
  n.cores = 8                             # Number of CPU cores to use for parallel processing
)
```

### 4. Visualize 

#### Display Only the Heatmap
To display only the heatmap of gene expression patterns over time without additional trajectory or enrichment analysis annotations, set `show.trajectory = FALSE` and `show.term = FALSE`. The `k.num` parameter determines the number of distinct expression patterns (clusters) to identify among the genes.


Explanation of Parameters:
`k.num`: The number of clusters (patterns) to identify among the genes. This parameter controls how many distinct expression patterns will be shown in the heatmap (Default is 6).
`show.trajectory`: Set to FALSE to hide trajectory information (Default is FALSE).
`show.term`: Set to FALSE to disable enrichment analysis using clusterProfiler. When show.term = TRUE, the function will perform enrichment analysis on the identified gene clusters (Default is FALSE).


```R
scITDGPlot(object = scitdg, 
           k.num = 5, 
           show.trajectory = FALSE, 
           show.term = FALSE, 
           save.wd = save.wd)
```

<p align="center">
  <img src="https://github.com/YandongZheng/scITDG/raw/main/images/Limb_Muscle_scITDGPlot_NoTrajectory_NoTerm_NoOrder.png" alt="Description" width=20%;" />
</p>


#### Display Heatmap with Gene Expression Trajectories

Setting `show.trajectory = TRUE` will fit the heatmap with expression trajectories, showing how gene expression changes over time.

```R
scITDGPlot(object = scitdg, 
           k.num = 5, 
           show.trajectory = TRUE, 
           show.term = FALSE, 
           save.wd = save.wd)
(-) [====================>------------------------------------------] 33% [ET: 00:00:04]
(\) [===============================>-------------------------------] 50% [ET: 00:00:07]
(|) [=========================================>---------------------] 67% [ET: 00:00:08]
(/) [===================================================>-----------] 83% [ET: 00:00:10]
(-) [===============================================================]100% [ET: 00:00:12]
```

<p align="center">
  <img src="https://github.com/YandongZheng/scITDG/raw/main/images/Limb_Muscle_scITDGPlot_NoTerm.png" alt="Description" width=20%;" />
</p>

#### Display Heatmap with Gene Expression Trajectories and Functional Enrichment Analysis

Additionally, scITDG can perform gene enrichment analysis for each pattern by setting show.term = TRUE, using [clusterProfiler](https://yulab-smu.top/biomedical-knowledge-mining-book/). Before this, you need to load species-specific gene annotation packages from Bioconductor (e.g., org.Hs.eg.db for human, org.Mm.eg.db for mouse, org.Dm.eg.db for Drosophila). Packages for other species can be loaded accordingly from [Bioconductor](https://bioconductor.org/packages/release/BiocViews.html#___OrgDb).


```R
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db")
```

```R
library(clusterProfiler)
library(DOSE)
library(org.Mm.eg.db)
library(ggExtra)
library(ggrastr)
library(Hmisc) 
library(cowplot)
library(data.table)
```


```R
scITDGPlot(object = scitdg, 
           k.num = 5, 
           show.trajectory = TRUE, 
           show.term = TRUE, 
           OrgDb = org.Mm.eg.db,
           save.wd = save.wd)
```
