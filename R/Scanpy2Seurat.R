#' Title
#'
#' @param file.name
#' @param data.type
#'
#' @return
#' @export
#'
#' @examples
Scanpy2Seurat <- function(file.name, data.type) {
  if (is.null(x = file.name)) {
    stop("Please input file!")
  }
  if (is.null(x = data.type)) {
    stop("Please input data.type!")
  }

  adata = scanpy$read_h5ad(file.name) ###载入scanpy输出的h5ad文件

  #######导出基因名和样本信息################
  # .obs存储meta信息，.var存储基因信息
  meta = adata$obs
  gene <- adata$var

  #############导出矩阵并转置，scanpy和Seurat的行列是反的#############
  # .X存储矩阵信息，.T存储的是转置信息，.T$X为转置矩阵

  # Normalized data
  mtx.data <- tryCatch({
    mtx <- pandas$DataFrame(data = adata$T$X$todense(),
                            index = adata$var_names, #行名
                            columns = adata$obs_names) # 列名
    return(mtx)
  }, error = function(e) {
    mtx <- pandas$DataFrame(data = adata$T$X %>% as.matrix,
                            index = adata$var_names, #行名
                            columns = adata$obs_names) # 列名
    return(mtx)
  })


  # raw count
  mtx.count <- tryCatch({
    mtx <- pandas$DataFrame(data = adata$raw$X$todense() %>% t,
                            index = adata$raw$var_names, #行名
                            columns = adata$raw$obs_names) # 列名
    return(mtx)
  }, error = function(e) {
    mtx <- pandas$DataFrame(data = adata$raw$X %>% as.matrix %>% t,
                            index = adata$raw$var_names, #行名
                            columns = adata$raw$obs_names) # 列名
    return(mtx)
  })


  ### 如果直接在Seurat4中使用的话会报错，先转为稠密矩阵，再转dgCMatrix############
  # Error in as(object = x, Class = "dgCMatrix") :
  #   no method or default for coercing “dgRMatrix” to “dgCMatrix”

  mtx.data %<>% as.matrix %>% as("dgCMatrix")
  mtx.count %<>% as.matrix %>% as("dgCMatrix")

  seurat.object <- CreateSeuratObject(counts = mtx.count)
  seurat.object <- AddMetaData(seurat.object, meta)
  seurat.object[["RNA"]][["n_cells"]] <- adata$var["n_cells"]

  normal.data <- CreateAssayObject(counts = mtx.data, min.cells = 0, min.features = 0)
  seurat.object@assays$RNA@data <- normal.data@data

  # Add embedding
  embedding <- adata$obsm["X_umap"]
  rownames(embedding) <- adata$obs_names$to_list()
  colnames(embedding) <- c("umap_1", "umap_2")
  seurat.object[["umap"]] <- CreateDimReducObject(embedding, key = "umap_",assay = "RNA")

  seurat.object$data.type <- data.type
  return(seurat.object)
}
