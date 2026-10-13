library(Seurat)

subset_seurat_by_umap <- function(
    seurat_obj, 
    celltype_col = "Celltype",  # customizable cell type column name (e.g. "celltype")
    target_celltype = "MC",     # target cell type (the part to be removed)
    umap1_condition = c("<", -5), 
    umap2_condition = c(">", -5),
    operator = "&"              # supports "&" (AND) or "|" (OR)
) {
  # Validate input
  if (!inherits(seurat_obj, "Seurat")) stop("Input must be a Seurat object.")
  if (!celltype_col %in% names(seurat_obj@meta.data)) {
    stop(paste0("Missing column in metadata: '", celltype_col, "'"))
  }
  if (!"umap" %in% names(seurat_obj@reductions)) stop("UMAP reduction not found in the Seurat object.")
  
  # Get all cell names
  all_cells <- colnames(seurat_obj)
  
  # Extract target cell type (cells to be checked against the removal criteria)
  cells_target <- all_cells[seurat_obj@meta.data[[celltype_col]] == target_celltype]
  
  # Non-target cells (kept as is)
  cells_non_target <- setdiff(all_cells, cells_target)
  
  # If there are no target cells, return the original object
  if (length(cells_target) == 0) {
    message("Target cell type not found: ", target_celltype)
    return(seurat_obj)
  }
  
  # Get UMAP coordinates of the target cells
  umap_data <- Embeddings(seurat_obj, "umap")[cells_target, ]
  umap1 <- umap_data[, "umap_1"]
  umap2 <- umap_data[, "umap_2"]
  
  # Parse the conditions
  op1 <- umap1_condition[1]
  val1 <- as.numeric(umap1_condition[2])
  op2 <- umap2_condition[1]
  val2 <- as.numeric(umap2_condition[2])
  
  # Build the logical conditions
  condition1 <- switch(
    op1,
    "<" = umap1 < val1,
    ">" = umap1 > val1,
    "==" = umap1 == val1,
    "!=" = umap1 != val1,
    stop("Unsupported operator: ", op1)
  )
  condition2 <- switch(
    op2,
    "<" = umap2 < val2,
    ">" = umap2 > val2,
    "==" = umap2 == val2,
    "!=" = umap2 != val2,
    stop("Unsupported operator: ", op2)
  )
  
  # Combine the conditions
  if (operator == "&") {
    cells_to_remove <- cells_target[condition1 & condition2]
  } else if (operator == "|") {
    cells_to_remove <- cells_target[condition1 | condition2]
  } else {
    stop("operator must be '&' or '|'")
  }
  
  # Cells to keep = non-target cells + target cells not removed
  cells_to_keep <- setdiff(all_cells, cells_to_remove)
  
  # Return the filtered Seurat object
  subset(seurat_obj, cells = cells_to_keep)
}