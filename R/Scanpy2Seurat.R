#' Title
#'
#' @param file.name
#' @param data.type
#'
#' @return
#' @export
#'
#' @examples
Scanpy2Seurat <- function(file.name, data.type = "droplet") {
  if (is.null(x = file.name)) {
    stop("Please input file!")
  }
  adata = scanpy$read_h5ad(file.name) 

  meta = adata$obs
  gene <- adata$var

  # output Normalized data
  mtx.data <- tryCatch({
    mtx <- pandas$DataFrame(data = adata$T$X$todense(),
                            index = adata$var_names, 
                            columns = adata$obs_names) 
    return(mtx)
  }, error = function(e) {
    mtx <- pandas$DataFrame(data = adata$T$X %>% as.matrix,
                            index = adata$var_names, 
                            columns = adata$obs_names) 
    return(mtx)
  })


  # output raw count
  mtx.count <- tryCatch({
    mtx <- pandas$DataFrame(data = adata$raw$X$todense() %>% t,
                            index = adata$raw$var_names, 
                            columns = adata$raw$obs_names) 
    return(mtx)
  }, error = function(e) {
    mtx <- pandas$DataFrame(data = adata$raw$X %>% as.matrix %>% t,
                            index = adata$raw$var_names, #行名
                            columns = adata$raw$obs_names) # 列名
    return(mtx)
  })

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
