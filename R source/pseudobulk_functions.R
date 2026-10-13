cell_proportion_by_sample <- function(
    seu,
    group_col   = "group",
    type_col    = "Type",
    cell_col    = "Cell_Anno"
) {
  # 1) Extract meta.data
  meta <- seu@meta.data
  
  # 2) Check required columns
  required_cols <- c(group_col, type_col, cell_col)
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in meta.data: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # 3) Number of cells per cell type in each sample
  cell_counts <- aggregate(
    x   = meta[[cell_col]],
    by  = list(
      group     = meta[[group_col]],
      Type      = meta[[type_col]],
      Cell_Anno = meta[[cell_col]]
    ),
    FUN = length
  )
  colnames(cell_counts)[4] <- "n"
  
  # 4) Total number of cells per sample
  total_counts <- aggregate(
    x   = cell_counts$n,
    by  = list(group = cell_counts$group),
    FUN = sum
  )
  colnames(total_counts)[2] <- "total"
  
  # 5) Merge and compute proportions
  cell_counts <- merge(cell_counts, total_counts, by = "group")
  cell_counts$prop <- cell_counts$n / cell_counts$total
  
  # 6) Reshape to wide format
  prop_df <- reshape(
    cell_counts[, c("group", "Type", "Cell_Anno", "prop")],
    idvar = c("group", "Type"),
    timevar = "Cell_Anno",
    direction = "wide"
  )
  
  colnames(prop_df) <- sub("^prop\\.", "", colnames(prop_df))
  prop_df[is.na(prop_df)] <- 0
  
  # 7) Sort
  prop_df <- prop_df[order(prop_df$group), ]
  
  return(prop_df)
}
################################################
################################################
################################################
plot_cell_proportion <- function(
    prop_df,
    cell_type,
    colors = c("#1B9E77", "#7570B3"),
    my_comparisons
) {
  library(ggplot2)
  library(ggpubr)
  
  ggplot(prop_df, aes(x = Type, y = .data[[cell_type]], fill = Type)) +
    theme_bw() +
    RotatedAxis() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    ) +
    labs(x = NULL, y = NULL, title = cell_type) +
    
    geom_jitter_rast(
      col = "#33333366",
      pch = 19,
      cex = 2,
      position = position_jitter(0.2)
    ) +
    
    geom_boxplot(position = position_dodge(0)) +
    scale_fill_manual(values = colors) +
    
    geom_boxplot(
      position = position_dodge(0),
      color = "black",
      outlier.colour = NA,
      outlier.fill = NA,
      outlier.shape = NA
    ) +
    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
    
    stat_compare_means(
      method = "t.test",
      hide.ns = FALSE,
      comparisons = my_comparisons,
      label = "p.signif",
      bracket.size = 0.8,
      size = 6
    )
}


############################################################
############################################################
############################################################



save_all_cell_proportions <- function(
    prop_df,
    colors = c("#1B9E77", "#7570B3"),
    my_comparisons,
    output_dir = "CellProportion_PNGs",
    width = 4,
    height = 5,
    dpi = 300
) {
  library(ggplot2)
  library(ggpubr)
  
  # All cell types (excluding non-numeric columns)
  cell_types <- setdiff(colnames(prop_df), c("group", "Type"))
  
  # Create the output directory
  dir.create(output_dir, showWarnings = FALSE)
  
  for (ct in cell_types) {
    
    p <- plot_cell_proportion(
      prop_df = prop_df,
      cell_type = ct,
      colors = colors,
      my_comparisons = my_comparisons
    )
    
    safe_name <- gsub("[[:space:]/]", "_", ct)
    
    ggsave(
      filename = file.path(output_dir, paste0(safe_name, "_proportion.png")),
      plot = p,
      width = width,
      height = height,
      dpi = dpi
    )
  }
}
############################################################
############################################################
############################################################
create_pseudobulk <- function(
    seu,
    group.by = "group",
    cell.type.col = NULL,
    cell.type = NULL,
    assay = "RNA",
    layer = "counts",
    min.total = 500,
    min.cells = 3,
    verbose = TRUE
) {
  # ---- checks ----
  if (!group.by %in% colnames(seu@meta.data)) {
    stop("group.by column '", group.by, "' not found in meta.data")
  }
  
  if (!is.null(cell.type.col) && !cell.type.col %in% colnames(seu@meta.data)) {
    stop("cell.type.col '", cell.type.col, "' not found in meta.data")
  }
  
  if (!assay %in% names(seu@assays)) {
    stop("Assay '", assay, "' not found in Seurat object")
  }
  
  # ---- Bypass the low-level safety extraction of GetAssayData entirely ----
  assay_obj <- seu[[assay]]
  mat <- NULL
  
  # 1. Try the standard Seurat v5 layers list or [[ extraction
  if (.hasSlot(assay_obj, "layers") && layer %in% names(assay_obj@layers)) {
    mat <- assay_obj@layers[[layer]]
    # Seurat v5 layers sometimes lack row/column names, so assign them manually
    if (is.null(rownames(mat))) rownames(mat) <- rownames(assay_obj)
    if (is.null(colnames(mat))) colnames(mat) <- colnames(assay_obj)
  } else if (!is.null(tryCatch(assay_obj[[layer]], error = function(e) NULL))) {
    mat <- assay_obj[[layer]]
  } else if (.hasSlot(assay_obj, layer)) {
    # 2. Fall back to the Seurat v4 slot (e.g. @counts)
    mat <- slot(assay_obj, layer)
  } else {
    stop("Cannot find layer/slot '", layer, "' in assay '", assay, "'")
  }
  
  # Ensure the column names match the cell barcodes exactly
  if (ncol(mat) != ncol(seu)) {
    colnames(mat) <- colnames(seu)
  }
  
  meta <- seu@meta.data
  
  # ---- filter by cell type (if requested) ----
  if (!is.null(cell.type.col) && !is.null(cell.type)) {
    keep.cells <- meta[[cell.type.col]] %in% cell.type
    mat <- mat[, keep.cells, drop = FALSE]
    meta <- meta[keep.cells, , drop = FALSE]
    
    if (verbose) {
      message("Filtered to cell type: ", paste(cell.type, collapse = ", "))
      message("Remaining cells: ", ncol(mat))
    }
  }
  
  # ---- build group key ----
  groups <- meta[[group.by]]
  
  if (verbose) {
    message("Aggregating pseudobulk by: ", group.by)
  }
  
  # ---- aggregate ----
  pb.list <- split(colnames(mat), groups)
  
  pb.mat <- sapply(pb.list, function(cell.set) {
    Matrix::rowSums(mat[, cell.set, drop = FALSE])
  })
  
  # ---- filter low-quality pseudobulk samples ----
  n.cells <- sapply(pb.list, length)
  tot.umi <- colSums(pb.mat)
  
  keep <- (tot.umi >= min.total) & (n.cells >= min.cells)
  
  pb.mat <- pb.mat[, keep, drop = FALSE]
  n.cells <- n.cells[keep]
  tot.umi <- tot.umi[keep]
  
  if (verbose) {
    message("Pseudobulk summary:")
    message("  Groups retained: ", ncol(pb.mat))
    message("  Min cells/group: ", if (length(n.cells) > 0) min(n.cells) else 0)
    message("  Max cells/group: ", if (length(n.cells) > 0) max(n.cells) else 0)
  }
  
  # ---- metadata ----
  pb.meta <- data.frame(
    group = colnames(pb.mat),
    n_cells = n.cells,
    total_umi = tot.umi,
    stringsAsFactors = FALSE
  )
  
  return(list(
    counts = pb.mat,
    meta   = pb.meta
  ))
}
#######################################
#####################################
##################################
plot_pseudobulk_feature <- function(seurat_obj, 
                                    feature, 
                                    assay = "Spatial",        
                                    sample_col = "group",     
                                    type_col = "Type",        
                                    comparisons_list = NULL,  
                                    min_cells_per_sample = 5, 
                                    colors = c('#F3AE63','#73558B','#93CBAE')) {
  
  require(dplyr)
  require(ggplot2)
  require(ggpubr)
  require(Seurat)
  
  # 1. Safety checks
  is_meta <- feature %in% colnames(seurat_obj@meta.data)
  
  is_gene <- FALSE
  if (assay %in% names(seurat_obj@assays)) {
    is_gene <- feature %in% rownames(seurat_obj[[assay]])
  }
  
  if (!is_meta && !is_gene) {
    stop(paste("Error: not found in metadata or the specified assay (", assay, "):", feature, "\nPlease check that the feature name is spelled correctly, or pass the correct assay parameter."))
  }
  
  pb_df <- data.frame()
  
  # Build an exact "sample -> group" lookup dictionary
  meta_sub <- unique(seurat_obj@meta.data[, c(sample_col, type_col)])
  mapping_dict <- setNames(as.character(meta_sub[[type_col]]), as.character(meta_sub[[sample_col]]))
  
  # ---- Case A: metadata column ----
  if (is_meta) {
    message(paste("--> Metadata column detected, aggregating the mean per sample for", feature, "..."))
    
    meta_data <- data.frame(
      Sample = seurat_obj@meta.data[[sample_col]],
      Type   = seurat_obj@meta.data[[type_col]],
      Value  = as.numeric(as.character(seurat_obj@meta.data[[feature]]))
    )
    
    # 💡[Fix] Force dplyr:: to prevent hijacking by other packages such as plyr
    pb_df <- meta_data %>%
      dplyr::group_by(Sample, Type) %>%
      dplyr::summarise(Expression = mean(Value, na.rm = TRUE), .groups = 'drop')
    
    y_label <- paste("Mean", feature, "(Sample level)")
    
    # ---- Case B: gene name ----
  } else {
    message(paste("--> Gene detected, computing pseudo-bulk (CPM) based on assay:", assay, "..."))
    
    cell_counts <- table(seurat_obj@meta.data[[sample_col]])
    valid_samples <- names(cell_counts)[cell_counts >= min_cells_per_sample]
    
    sub_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj@meta.data[[sample_col]] %in% valid_samples])
    
    cell_totals <- colSums(sub_obj[[assay]]$counts) 
    true_sample_totals <- tapply(cell_totals, sub_obj@meta.data[[sample_col]], sum) 
    
    pb <- AggregateExpression(
      sub_obj, 
      assays = assay, 
      features = feature, 
      return.seurat = FALSE, 
      group.by = sample_col,
      slot = "counts"
    )
    
    pb_matrix <- pb[[assay]]
    
    pb_df <- data.frame(
      Sample = colnames(pb_matrix),
      Raw_Count = as.numeric(pb_matrix[feature, ]),
      stringsAsFactors = FALSE
    ) %>%
      dplyr::mutate(
        Type = mapping_dict[Sample] 
      )
    
    pb_df$Expression <- log2((pb_df$Raw_Count / true_sample_totals[pb_df$Sample]) * 10^4 + 1)
    
    y_label <- paste(feature, "Expression (log2 CPM)")
  }
  
  # ---- 2. Dynamic axis space calculation ----
  y_max_val <- max(pb_df$Expression, na.rm = TRUE)
  y_min_val <- min(pb_df$Expression, na.rm = TRUE)
  y_range   <- y_max_val - y_min_val
  if (is.na(y_range) || y_range == 0) y_range <- 1 
  
  start_label_y <- y_max_val + (y_range * 0.08)
  num_comps <- ifelse(is.null(comparisons_list), 0, length(comparisons_list))
  extended_y_max <- y_max_val + (y_range * (0.12 + 0.12 * num_comps))
  
  # 💡[Key fix] Drop the conflict-prone filter() and use the most native, safest Base R subsetting
  pb_df <- pb_df[!is.na(pb_df$Type), ]
  pb_df$Type <- factor(trimws(pb_df$Type))
  
  # ---- 3. ggplot2 plotting ----
  p <- ggplot(pb_df, aes(x = Type, y = Expression, fill = Type)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5, linewidth = 0.8, color = "gray20") + 
    geom_jitter(width = 0.15, height = 0, size = 3, shape = 21, color = "black", stroke = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = colors) +
    ylim(y_min_val - (y_range * 0.05), extended_y_max) +
    labs(title = paste("Pseudo-bulk:", feature), x = NULL, y = y_label) +
    theme_bw() +
    theme(
      plot.margin = margin(t = 20, r = 15, b = 15, l = 15, unit = "pt"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title.y = element_text(size = 11, face = "bold"),
      axis.text = element_text(size = 10, colour = "black"),
      panel.grid = element_blank(),
      legend.position = "none"
    )
  
  # ---- 4. Add statistical test brackets ----
  if (!is.null(comparisons_list) && nrow(pb_df) > 0) {
    step_increase <- y_range * 0.13
    label_positions <- seq(from = start_label_y, by = step_increase, length.out = length(comparisons_list))
    
    p <- p + stat_compare_means(
      comparisons = comparisons_list, 
      method = "t.test", 
      label = "p.signif", 
      label.y = label_positions,
      bracket.size = 0.6,
      size = 5,
      symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns"))
    )
  }
  
  return(p)
}
#####################
####################
######################
library(edgeR)
library(limma)
library(dplyr)
library(edgeR)
library(DESeq2)

run_pseudobulk_dea <- function(pb_output, 
                               sample_mapping, 
                               contrast_str,       # pass a single contrast at a time, e.g. "Y - M"
                               method = c("limma", "DESeq2", "edgeR"), # method selection
                               filter_genes = TRUE) {
  
  # Match the selected method, default limma
  method <- match.arg(method)
  
  exprSet <- as.matrix(pb_output$counts)
  pb_meta <- pb_output$meta
  
  # 1. Check that the metadata and expression matrix are aligned
  if (!all(colnames(exprSet) == pb_meta$group)) {
    stop("Error: column names of the expression matrix do not match the order of pb_meta$group!")
  }
  
  # 2. Map sample names to the real experimental groups
  if (is.data.frame(sample_mapping)) {
    mapping_vec <- setNames(as.character(sample_mapping[[2]]), as.character(sample_mapping[[1]]))
  } else {
    mapping_vec <- sample_mapping
  }
  
  mapped_groups <- mapping_vec[pb_meta$group]
  if (any(is.na(mapped_groups))) {
    stop("Error: sample_mapping is missing group information for some samples!")
  }
  
  group_list <- factor(mapped_groups)
  
  # Parse the contrast string (extract the Up and Down groups for DESeq2 and console printing)
  clean_contrast <- gsub(" ", "", contrast_str) 
  if (grepl("-", clean_contrast)) {
    groups_split <- unlist(strsplit(clean_contrast, "-"))
    up_group <- groups_split[1]
    down_group <- groups_split[2]
  } else {
    up_group <- clean_contrast
    down_group <- "control"
  }
  
  # 3. Gene filtering (common step)
  dge <- DGEList(counts = exprSet, group = group_list)
  if (filter_genes) {
    # 💡[Fine-tuning] Relax the filtering threshold for single-cell/spatial pseudo-bulk
    keep <- filterByExpr(dge, group = group_list, min.count = 2, min.total.count = 10)
    dge <- dge[keep, , keep.lib.sizes = FALSE]
    exprSet <- exprSet[keep, , drop = FALSE] # keep the raw matrix in sync for DESeq2
  }
  
  # ==========================================
  # 4. Dispatch the differential analysis according to the selected method
  # ==========================================
  
  if (method == "limma") {
    
    dge <- calcNormFactors(dge) # TMM normalization
    design <- model.matrix(~ 0 + group_list)
    colnames(design) <- levels(group_list)
    
    v <- voom(dge, design, plot = FALSE)
    
    # 💡 Improved the way arguments are passed to makeContrasts to avoid the pitfalls of eval
    contrast.matrix <- makeContrasts(contrasts = clean_contrast, levels = design)
    
    fit <- lmFit(v, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    
    target_coef <- colnames(contrast.matrix)[1] 
    tempOutput  <- topTable(fit2, coef = target_coef, n = Inf)
    nrDEG       <- na.omit(tempOutput)
    
  } else if (method == "edgeR") {
    
    dge <- calcNormFactors(dge) # TMM normalization
    design <- model.matrix(~ 0 + group_list)
    colnames(design) <- levels(group_list)
    
    # Estimate dispersion and fit a generalized linear model (QLF)
    dge <- estimateDisp(dge, design)
    fit <- glmQLFit(dge, design)
    
    contrast.matrix <- makeContrasts(contrasts = clean_contrast, levels = design)
    target_coef <- colnames(contrast.matrix)[1]
    
    qlf <- glmQLFTest(fit, contrast = contrast.matrix[, 1])
    tempOutput <- topTags(qlf, n = Inf)$table
    nrDEG      <- na.omit(tempOutput)
    
  } else if (method == "DESeq2") {
    
    # Check DESeq2-specific format support
    if (!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("Please install the DESeq2 package first! (BiocManager::install('DESeq2'))")
    }
    
    # DESeq2 needs a non-normalized integer count matrix and estimates size factors internally
    colData <- data.frame(group = group_list, row.names = colnames(exprSet))
    
    # Coerce floating-point counts (which may contain decimals when extracted from Seurat) to integers
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(exprSet), 
                                          colData = colData, 
                                          design = ~ group)
    
    dds <- DESeq2::DESeq(dds, quiet = TRUE)
    
    # Extract results; the contrast order is: c("Condition", "numerator/treatment group", "denominator/control group")
    res <- DESeq2::results(dds, contrast = c("group", up_group, down_group))
    tempOutput <- as.data.frame(res)
    nrDEG <- na.omit(tempOutput)
  }
  
  # 5. Clean console output
  message("==========================================================")
  message(sprintf("📊 Differential analysis contrast [%s] computed successfully!", clean_contrast))
  message(sprintf("⚙️  Statistical model used: %s", method))
  message(sprintf("💡 [Sign interpretation]: a positive logFC 🟢 means higher expression in the %s group", up_group))
  message(sprintf("                a negative logFC 🔴 means higher expression in the %s group", down_group))
  message(sprintf("🔍 Significant genes produced (P value not filtered): %d", nrow(nrDEG)))
  message("==========================================================")
  
  print(head(nrDEG))
  
  return(nrDEG)
}

##########################
##########################
#########################
#' Volcano plot for differential analysis results
#'
#' @param deg_result Data frame output by limma (must contain logFC and adj.P.Val columns)
#' @param logFC_cutoff Absolute value threshold for logFC, default 1
#' @param pval_cutoff Significance threshold for adj.P.Val, default 0.05
#' @param title_prefix Prefix text for the title (optional)
#'
#' @return Returns a ggplot object and prints the plot automatically
plot_volcano <- function(deg_result, logFC_cutoff = 1, pval_cutoff = 0.05, title_prefix = "") {
  
  suppressMessages(library(ggplot2))
  
  # Copy the data to avoid modifying the original, and ensure it is a data.frame
  DEG <- as.data.frame(deg_result)
  
  # ==========================================
  # 0. Automatically adapt to the column names of the three tools (limma, edgeR, DESeq2)
  # ==========================================
  # Adapt the logFC column
  if ("log2FoldChange" %in% colnames(DEG)) {
    DEG$logFC <- DEG$log2FoldChange
  } else if (!"logFC" %in% colnames(DEG)) {
    stop("Error: neither logFC nor log2FoldChange column found, please check the input data!")
  }
  
  # Adapt the p-value / FDR column (prefer the adjusted p value)
  if ("padj" %in% colnames(DEG)) {
    DEG$adj.P.Val <- DEG$padj
  } else if ("FDR" %in% colnames(DEG)) {
    DEG$adj.P.Val <- DEG$FDR
  } else if ("P.Value" %in% colnames(DEG) && !"adj.P.Val" %in% colnames(DEG)) { 
    # Fallback: if only the unadjusted P.Value is present (limma)
    DEG$adj.P.Val <- DEG$P.Value
  } else if ("PValue" %in% colnames(DEG) && !"adj.P.Val" %in% colnames(DEG)) { 
    # Fallback: if only the unadjusted PValue is present (edgeR)
    DEG$adj.P.Val <- DEG$PValue
  } else if ("pvalue" %in% colnames(DEG) && !"adj.P.Val" %in% colnames(DEG)) { 
    # Fallback: if only the unadjusted pvalue is present (DESeq2)
    DEG$adj.P.Val <- DEG$pvalue
  } else if (!"adj.P.Val" %in% colnames(DEG)) {
    stop("Error: none of the p value/FDR/padj/adj.P.Val columns found, please check the input data!")
  }
  
  # 💡[Key fix] DESeq2 results often have NA padj (lowly expressed genes are filtered out automatically)
  # Rows containing NA must be removed, otherwise ifelse or ggplot plotting will fail
  DEG <- subset(DEG, !is.na(adj.P.Val) & !is.na(logFC))
  
  # ==========================================
  
  # 1. Classify genes by threshold (UP, DOWN, NOT)
  DEG$change <- factor(
    ifelse(DEG$adj.P.Val < pval_cutoff & abs(DEG$logFC) > logFC_cutoff,
           ifelse(DEG$logFC > logFC_cutoff, 'UP', 'DOWN'), 
           'NOT'),
    levels = c('DOWN', 'NOT', 'UP') # Force a fixed order so the color mapping stays correct
  )
  
  # 2. Count the up- and down-regulated genes
  num_up <- sum(DEG$change == 'UP', na.rm = TRUE)
  num_down <- sum(DEG$change == 'DOWN', na.rm = TRUE)
  
  # 3. Build the dynamic title
  this_title <- paste0(title_prefix,
                       'Cutoff for logFC is ', round(logFC_cutoff, 3),
                       '\nThe number of up genes is ', num_up,
                       '\nThe number of down genes is ', num_down)
  
  # 4. Define the color dictionary (DOWN=blue, NOT=grey, UP=red)
  color_mapping <- c("DOWN" = "#42AFFC", "NOT" = "grey", "UP" = "tomato")
  
  # 5. Plot with ggplot2
  g <- ggplot(data = DEG, aes(x = logFC, y = -log10(adj.P.Val), color = change)) +
    geom_point(alpha = 0.6, size = 2) +  # slightly reduced alpha and size to avoid over-dense points
    theme_bw(base_size = 12) +
    xlab("log2 fold change") + 
    ylab("-log10(adj.P.Val / FDR)") +  # dynamically adjust the Y axis label to cover FDR
    ggtitle(this_title) + 
    theme(
      plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) +
    scale_colour_manual(values = color_mapping) +
    # Add dashed guide lines for the thresholds
    geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), linetype = "dashed", color = "black", alpha = 0.4) +
    geom_hline(yintercept = -log10(pval_cutoff), linetype = "dashed", color = "black", alpha = 0.4)
  
  # Print the plot to the device
  print(g)
  
  # Return the object so it can be saved with ggsave as PDF/PNG
  invisible(g)
}
##############################
############################
##########################
plot_sc_violin_contrast <- function(seurat_obj, 
                                    feature, 
                                    group_col = "Type", 
                                    comparisons_list = NULL, 
                                    colors = c('#F3AE63','#73558B','#93CBAE')) {
  
  # 1. Automatically check whether feature exists in genes or metadata
  # FetchData is powerful and recognizes gene names or meta.data column names automatically
  all_vars <- c(group_col, feature)
  
  df1 <- tryCatch({
    FetchData(seurat_obj, vars = all_vars)
  }, error = function(e) {
    stop(paste("Error: the specified feature or group column was not found in the object. Please check whether", feature, "or", group_col, "is correct."))
  })
  
  # Rename the columns so ggplot recognizes them consistently
  colnames(df1) <- c("Group", "Value")
  
  # 2. Base violin plot
  p <- ggplot(df1, aes(x = Group, y = Value, fill = Group)) +
    geom_violin(color = 'black', size = 1) + 
    theme_classic() +   
    theme(text = element_text(size = 10, colour = "black")) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 15, face = "bold"), 
      axis.text.x = element_text(colour = "black", size = 12),
      axis.text.y = element_text(colour = "black", size = 10),
      axis.title.y = element_text(color = 'black', size = 12),
      axis.line = element_line(size = 1),
      legend.position = "none"
    ) +  
    scale_fill_manual(values = colors) + 
    labs(title = paste("Single-cell Contrast:", feature), x = group_col, y = feature)
  
  # 3. Overlay a mini boxplot
  p <- p + geom_boxplot(
    width = 0.05, 
    fill = "#F5F5F5", 
    color = "gray30", 
    linewidth = 0.8, 
    alpha = 0.8,
    outlier.shape = NA
  )
  
  # 4. Dynamically add significance comparisons
  if (!is.null(comparisons_list)) {
    # Compute the data range for smart step increases (avoids overlapping brackets with multiple comparisons)
    y_max_val <- max(df1$Value, na.rm = TRUE)
    y_min_val <- min(df1$Value, na.rm = TRUE)
    y_range   <- y_max_val - y_min_val
    if (y_range == 0) y_range <- 1
    
    start_label_y <- y_max_val + (y_range * 0.05)
    step_increase <- y_range * 0.12
    
    # Dynamically compute the y positions for each comparison
    label_positions <- seq(from = start_label_y, 
                           by = step_increase, 
                           length.out = length(comparisons_list))
    
    p <- p + stat_compare_means(
      method = "t.test",
      hide.ns = FALSE,  
      comparisons = comparisons_list, 
      label = "p.signif",                
      bracket.size = 0.8,   
      size = 6,
      label.y = label_positions # key fix: prevents brackets from overlapping in multi-group comparisons
    )
  }
  
  # 5. Expand the flexible space above the y axis
  p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  
  return(p)
}


##########################
# Correlation between genes and metadata (pseudo-bulk)
ks_pseudobulk_cor <- function(seurat_obj, 
                              var1, 
                              var2, 
                              group_by = "group", 
                              type_by = "Type",
                              cor_method = "pearson",
                              point_palette = NULL,
                              title = NULL) {
  library(Seurat)
  library(ggplot2)
  library(ggpubr)
  
  # 1. Check and extract the first variable (var1)
  if (var1 %in% colnames(seurat_obj@meta.data)) {
    val1 <- seurat_obj@meta.data[[var1]]
  } else if (var1 %in% rownames(seurat_obj)) {
    val1 <- as.numeric(FetchData(seurat_obj, vars = var1)[, 1])
  } else {
    stop(paste("var1 not found in metadata or genes:", var1))
  }
  
  # 2. Check and extract the second variable (var2)
  if (var2 %in% colnames(seurat_obj@meta.data)) {
    val2 <- seurat_obj@meta.data[[var2]]
  } else if (var2 %in% rownames(seurat_obj)) {
    val2 <- as.numeric(FetchData(seurat_obj, vars = var2)[, 1])
  } else {
    stop(paste("var2 not found in metadata or genes:", var2))
  }
  
  # 3. Extract the sample (Group) and experimental group (Type) information
  if (!group_by %in% colnames(seurat_obj@meta.data)) stop(paste("Group column not found:", group_by))
  if (!type_by %in% colnames(seurat_obj@meta.data)) stop(paste("Type column not found:", type_by))
  
  groups <- seurat_obj@meta.data[[group_by]]
  types <- seurat_obj@meta.data[[type_by]]
  
  # 4. Build the single-cell level data frame
  sc_df <- data.frame(
    V1 = val1,
    V2 = val2,
    Group = groups,
    Type = types,
    stringsAsFactors = FALSE
  )
  
  # 5. Aggregate to pseudo-bulk with the base split mechanism (immune to all package conflicts!)
  bulk_list <- split(sc_df, list(sc_df$Group, sc_df$Type), drop = TRUE)
  
  bulk_df <- do.call(rbind, lapply(bulk_list, function(sub_df) {
    data.frame(
      Group = sub_df$Group[1],
      Type = sub_df$Type[1],
      Mean_V1 = mean(sub_df$V1, na.rm = TRUE),
      Mean_V2 = mean(sub_df$V2, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  
  # 6. Handle the title and axis labels
  plot_title <- if(is.null(title)) paste("Pseudobulk Correlation:", var1, "vs", var2) else title
  
  # 7. Draw the point plot in theme_classic style
  p <- ggplot(bulk_df, aes(x = Mean_V1, y = Mean_V2)) +
    # Trend fitting line (fit on all samples)
    geom_smooth(method = "lm", color = "#333333", fill = "#E6F5D0", alpha = 0.2, size = 1) +
    # Sample scatter points (no labels)
    geom_point(aes(color = Type), size = 4.5, alpha = 0.9) +
    # Correlation statistics annotation
    stat_cor(method = cor_method, label.x.npc = "left", label.y.npc = "top", size = 4.5) +
    # Classic academic theme
    theme_classic() +
    labs(
      x = paste("Pseudobulk Mean of", var1),
      y = paste("Pseudobulk Mean of", var2),
      title = plot_title,
      color = type_by
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black", size = 10),
      legend.title = element_text(face = "bold")
    )
  
  # 8. Custom color mapping
  if (!is.null(point_palette)) {
    p <- p + scale_color_manual(values = point_palette)
  }
  
  return(p)
}
###########################
######################
## Compute the convex hull area of regional cells to assess cell density
#' Compute and compare spatial macro-region density (supports merging multiple regions)
CompareRegionalDensity <- function(seu_obj, 
                                   celltype_col, target_celltype, 
                                   region_col = NULL, target_region = NULL,
                                   sample_col = "group", condition_col = "Type",
                                   reduction = "spatial",
                                   comparisons = NULL,                 
                                   test_method = "t.test",             
                                   plot_colors = c('#FF9AA2', '#B5EAD7', '#A0E7E5', '#FFB7B2')) {
  
  if (!(reduction %in% names(seu_obj@reductions))) stop(paste("Reduction named", reduction, "not found."))
  
  meta <- seu_obj@meta.data
  coords <- Embeddings(seu_obj, reduction = reduction)
  sp_x <- colnames(coords)[1]
  sp_y <- colnames(coords)[2]
  df <- cbind(meta, coords)
  
  # ✨ Key improvement 1: use %in% to support passing multiple regions (e.g. c("niche_2", "niche_3"))
  if (!is.null(region_col) && !is.null(target_region)) {
    df <- df %>% filter(.data[[region_col]] %in% target_region)
    if (nrow(df) == 0) stop("No cells found in the specified anatomical region!")
  }
  
  calc_polygon_area <- function(x, y) {
    if (length(x) < 3) return(NA)
    hull_idx <- chull(x, y) 
    hx <- x[hull_idx]
    hy <- y[hull_idx]
    n <- length(hx)
    return(0.5 * abs(sum(hx[-n] * hy[-1]) + hx[n] * hy[1] - sum(hy[-n] * hx[-1]) - hy[n] * hx[1]))
  }
  
  samples <- unique(df[[sample_col]])
  result_list <- list()
  
  for (samp in samples) {
    samp_df <- df %>% filter(.data[[sample_col]] == samp)
    if(nrow(samp_df) == 0) next
    
    cond <- samp_df[[condition_col]][1]
    area <- calc_polygon_area(samp_df[[sp_x]], samp_df[[sp_y]])
    count <- sum(samp_df[[celltype_col]] == target_celltype, na.rm = TRUE)
    density <- count / area
    
    result_list[[samp]] <- data.frame(Sample = samp, Condition = cond, Area = area, Count = count, Density = density)
  }
  
  result_df <- do.call(rbind, result_list) %>% filter(!is.na(Area))
  
  # ✨ Key improvement 2: elegantly show the merged region names in the plot title
  region_title <- ifelse(!is.null(target_region), 
                         paste("Region:", paste(target_region, collapse = " + ")), 
                         "Whole Sample")
  
  p <- ggplot(result_df, aes(x = Condition, y = Density, fill = Condition)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.2, size = 2, alpha = 0.8, color = "black") +
    scale_fill_manual(values = plot_colors) +
    labs(
      title = paste0(target_celltype, " Density Difference"),
      subtitle = region_title, # use the new merged title
      y = "Density",
      x = ""
    ) +
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) 
  
  if (var(result_df$Density) > 0 && length(unique(result_df$Condition)) >= 2) {
    if (!is.null(comparisons)) {
      p <- p + stat_compare_means(comparisons = comparisons, method = test_method, 
                                  label = "p.signif", tip.length = 0.02)
    } else {
      p <- p + stat_compare_means(method = test_method, label = "p.signif", label.x.npc = "center")
    }
  } else {
    message("⚠️ Warning: no statistical test performed because the density variance is 0 or only one group is present.")
  }
  
  return(list(Data = result_df, Plot = p))
}
#########################
ks_plot_density_summary <- function(obj, 
                                    feature = "factor7", 
                                    group_by = "Type", 
                                    custom_colors = NULL,
                                    y_limits = NULL,     # optional: limit the Y axis range, e.g. c(0, 5)
                                    x_limits = NULL) {   # optional: limit the X axis range, e.g. c(0, 10)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  
  # 1. Auto-compatible: accepts a Seurat object directly or a meta data frame
  if (inherits(obj, "Seurat")) {
    meta_df <- obj@meta.data
  } else if (is.data.frame(obj)) {
    meta_df <- obj
  } else {
    stop("Invalid input object: must be a Seurat object or a metadata data.frame")
  }
  
  # Check whether the column names exist
  if (!feature %in% colnames(meta_df)) stop(paste("Feature column not found in the data:", feature))
  if (!group_by %in% colnames(meta_df)) stop(paste("Group column not found in the data:", group_by))
  
  # 2. Summary statistics computation (fully aligned to the specified feature)
  summary_table <- meta_df %>%
    dplyr::group_by(across(all_of(group_by))) %>%
    dplyr::summarise(
      mean_score   = mean(.data[[feature]], na.rm = TRUE),
      median_score = median(.data[[feature]], na.rm = TRUE),
      sd_score     = sd(.data[[feature]], na.rm = TRUE),
      cell_count   = dplyr::n(),
      .groups = "drop"
    )
  
  # Print the tidy summary statistics to the console
  cat("\n📊 --- Summary statistics of", feature, "across", group_by, "groups ---\n")
  print(as.data.frame(summary_table))
  cat("----------------------------------------------------\n\n")
  
  # 3. Draw an academic-style density plot
  p <- ggplot(meta_df, aes(x = .data[[feature]], fill = .data[[group_by]])) +
    geom_density(alpha = 0.4, color = "#333333", linewidth = 0.5) + # note: in recent ggplot versions size should be replaced with linewidth
    labs(
      x = paste(feature, "Score"),
      y = "Density",
      title = paste("Density Distribution of", feature),
      fill = group_by
    ) +
    theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
      axis.title = element_text(face = "bold", size = 11),
      axis.text = element_text(color = "black", size = 10),
      legend.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  
  # 4. Custom color control
  if (!is.null(custom_colors)) {
    p <- p + scale_fill_manual(values = custom_colors)
  } else {
    p <- p + scale_fill_brewer(palette = "Set2")
  }
  
  # 5. Axis limiting (Zoom in - key modification)
  # Use coord_cartesian to zoom in safely without dropping data points
  if (!is.null(y_limits) || !is.null(x_limits)) {
    p <- p + coord_cartesian(ylim = y_limits, xlim = x_limits)
  }
  
  return(p)
}
