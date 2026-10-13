#' Plot a UMAP with unchull
#'
#' @param seurat_obj Seurat object
#' @param reduction Dimensional reduction method, default "umap"
#' @param cluster_col Name of the column holding the cluster labels (e.g. "Type")
#' @param colors Color vector (must match the cluster labels)
#' @param point_size Point size, default 0.01
#' @param unchull_alpha unchull transparency, default 0.00
#' @param unchull_size unchull border size, default 0.25
#' @param unchull_delta unchull range parameter, default 0.3
#' @param axis_trunc Axis truncation length (unit: cm), default 2
#' @param aspect_ratio Aspect ratio, default 1
#' @param show_legend Whether to draw the legend, default TRUE
#' @param arrow_length Axis arrow length (unit: inches), default 0.2
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' # Example usage
#' p <- plot_umap_with_unchull(
#'   seurat_obj = afd,
#'   cluster_col = "Type",
#'   colors = TNBC_color,
#'   point_size = 0.3
#' )
#' print(p)
plot_umap_with_unchull <- function(
    seurat_obj,
    reduction = "umap",
    cluster_col = "Type",
    colors = NULL,
    point_size = 0.01,
    unchull_alpha = 0.00,
    unchull_size = 0.25,
    unchull_delta = 0.3,
    axis_trunc = 2,
    aspect_ratio = 1,
    show_legend = TRUE,
    arrow_length = 0.2
) {
  # Check required packages
  if (!requireNamespace("ggunchull", quietly = TRUE)) {
    stop("Please install the ggunchull package first: devtools::install_github('const-ae/ggunchull')")
  }
  if (!requireNamespace("ggh4x", quietly = TRUE)) {
    stop("Please install the ggh4x package first: install.packages('ggh4x')")
  }
  
  # Check input
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop("Specified reduction not found in the Seurat object: ", reduction)
  }
  if (!cluster_col %in% names(seurat_obj@meta.data)) {
    stop("Specified cluster column not found in metadata: ", cluster_col)
  }
  axis <- ggh4x::guide_axis_truncated(
    trunc_lower = unit(0, "npc"),
    trunc_upper = unit(3, "cm")
  )
  library(ggunchull)
  library(ggh4x)
  # Extract UMAP coordinates and cluster labels
  plotData <- as.data.frame(seurat_obj[[reduction]]@cell.embeddings)
  Cluster <- seurat_obj@meta.data[[cluster_col]]
  colnames(plotData) <- c("umap_1", "umap_2")
  # Check the color vector
  if (is.null(colors)) {
    n_clusters <- length(unique(Cluster))
    colors <- scales::hue_pal()(n_clusters)
    warning("No color vector provided, using the default ggplot2 palette")
  }
  
  # Create the base plot
  p <- ggplot(plotData, aes(x = .data[[paste0(reduction, "_1")]], 
                            y = .data[[paste0(reduction, "_2")]], 
                            fill = Cluster, 
                            color = Cluster)) +
    stat_unchull(
      alpha = unchull_alpha, 
      size = unchull_size, 
      lty = 0, 
      delta = unchull_delta
    ) +
    geom_point(
      size = point_size, 
      show.legend = show_legend
    ) +
    theme(
      aspect.ratio = aspect_ratio,
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(
        arrow = arrow(type = "closed", length = unit(arrow_length, "inches"))
      ),
      axis.title = element_text(hjust = 0.05, face = "italic")
    ) +
    guides(
      color = guide_legend(override.aes = list(alpha = 1)),
      fill = guide_legend(override.aes = list(alpha = 1))
    ) +
    scale_x_continuous(
      breaks = NULL,
      guide = guide_axis_truncated(
        trunc_lower = unit(0, "npc"), 
        trunc_upper = unit(axis_trunc, "cm")
      )
    ) +
    scale_y_continuous(
      breaks = NULL,
      guide = guide_axis_truncated(
        trunc_lower = unit(0, "npc"), 
        trunc_upper = unit(axis_trunc, "cm")
      )
    )
  
  # Add color scales
  if (!is.null(colors)) {
    p <- p + 
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors)
  }
  
  # Set axis labels
  p <- p + 
    xlab(paste0(toupper(reduction), "_1")) +
    ylab(paste0(toupper(reduction), "_2"))
  
  return(p)
}


plotUMAP <- function(seurat_obj, reduction = "umap", cluster_col, color_palette, 
                     axis_line_size = 0.8, arrow_length = 0.1, point_size = 0.03,
                     legend_point_size = 3, legend_text_size = 5,
                     legend_ncol = 1, legend_alpha = 0.8) {
  library(ggplot2)
  library(cowplot)
  
  # Extract data
  df <- as.data.frame(seurat_obj[[reduction]]@cell.embeddings)
  colnames(df) <- c("UMAP_1", "UMAP_2")
  df$Cluster <- as.factor(seurat_obj@meta.data[[cluster_col]])
  
  # Create the main plot
  p_main <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = point_size, shape = 16, stroke = 0) +
    scale_color_manual(values = color_palette) +
    theme_classic() +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_line(colour = "black", size = axis_line_size, 
                               arrow = arrow(length = unit(arrow_length, "cm"))),
      legend.position = "none",
      aspect.ratio = 1
    ) +
    scale_x_continuous(limits = c(min(df$UMAP_1), max(df$UMAP_1))) +
    scale_y_continuous(limits = c(min(df$UMAP_2), max(df$UMAP_2)))
  
  # Create a separate legend
  p_legend <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = Cluster)) +
    geom_point(size = 0, shape = 16, stroke = 0, alpha = 0) +
    theme_void() +
    scale_color_manual(values = color_palette, name = '') +
    theme(
      legend.position = 'right',
      plot.margin = margin(0, 0, 0, 0)
    ) +
    guides(color = guide_legend(ncol = legend_ncol, 
                                override.aes = list(size = legend_point_size, 
                                                    shape = 16, 
                                                    alpha = legend_alpha))) +
    theme(
      legend.text = element_text(size = legend_text_size),
      legend.spacing.y = unit(0.1, 'cm'),
      legend.spacing.x = unit(0.1, 'cm'),
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.3, "cm"),
      legend.box.spacing = unit(0, 'cm'),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0)
    )
  
  # Extract the legend
  legend <- get_legend(p_legend)
  
  # Combine the plots, legend on the right
  combined_plot <- plot_grid(
    p_main, 
    legend, 
    ncol = 2, 
    rel_widths = c(0.85, 0.15)
  )
  
  return(combined_plot)
}
