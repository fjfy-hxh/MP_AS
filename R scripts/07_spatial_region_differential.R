# =============================================================================
# 07_spatial_region_differential.R
# -----------------------------------------------------------------------------
# Spatial regional cell proportion and regional cell density differential analysis
#
# Functions
#   - Compute specified niche / region cell-type proportions per sample and plot group-wise differential boxplots
#   - Batch-loop output of proportion comparison plots for all cell types
#   - Pseudobulk single-gene inter-group comparison (plot_pseudobulk_feature)
#   - Pseudobulk correlation between gene expression and metadata scores (e.g. pathway scores)
#   - Single-cell level multi-group pairwise comparison violin plots
#   - Compute target-cell convex hull area / count / density by region and compare across groups
#
# Main inputs
#   - Seurat spatial object (with niche, Cell_Anno / Celltype, group, Type metadata, etc.)
#
# Main outputs
#   - Cell proportion boxplots, niche score violin plots, fibroblast_density_diff.png
#
# Dependencies
#   In-repo custom functions: pseudobulk_functions.R
#   Main R packages: Seurat, ggpubr, ggplot2
#
# Notes
#   This script was organized from analysis drafts, keeping the original analysis logic intact with no changes to algorithms or workflow.
#   Before running, adjust the input / output paths as needed (the original drafts mostly used local absolute paths),
#   and replace the example sample names, cell type names, and comparison groups with your own data.
#   It is recommended to first set the repository root as the working directory (open
#   sc-spatial-analysis-scripts.Rproj in RStudio) so that source("./R/...") resolves correctly.
# =============================================================================


source(file = "./R/pseudobulk_functions.R")
library(ggpubr)
#########################
table(af.ST$group)
af.cellpro<- subset(af.ST, niche %in%  c("niche_1",'niche_2','niche_3','niche_5'))
#af.epi$subcluster
prop_df <- cell_proportion_by_sample(
  seu      = af.cellpro,
  group_col = "group",#sample
  type_col  = "Type",#grouping basis
  cell_col  = "Cell_Anno"#####Y axis
)

#set comparison groups
my_comparisons1 <- list(c("A", "M")) 
my_comparisons2 <- list(c("A", "Y"))
my_comparisons3 <- list(c("M", "Y"))

colors<-c('#F3AE63','#73558B','#93CBAE')

table(af.ST$Cell_Anno)
#make boxplot
p <- plot_cell_proportion(
  prop_df = prop_df,
  cell_type = "Epi",
  colors=c('#F3AE63','#73558B','#93CBAE'),
  my_comparisons = c(my_comparisons1,my_comparisons2,my_comparisons3)
)
p
ggsave("Fibroblasts niche1 3 5 ration.pdf",width = 5,height = 4.5)

##########################
###loop output for all cells
##########################
#all cells except these two columns
cell_types <- setdiff(colnames(prop_df), c("group", "Type"))
save_all_cell_proportions(
  prop_df = prop_df,
  colors=c('#F3AE63','#73558B','#93CBAE'),
  my_comparisons = c(my_comparisons1,my_comparisons2,my_comparisons3)
)

##################################################################
af.fib$subcluster
##############################################################################
##############################################################################
##############################################################################
af.pseubluk<- subset(af.fib, subcluster %in% c('Fib_1'))
table(af.pseubluk$Type)
# define your comparison group list
my_comparisons1 <- list(c("A", "M")) 
my_comparisons2 <- list(c("A", "Y"))
my_comparisons3 <- list(c("M", "Y"))
my_colors <- c('#F3AE63', '#73558B', '#93CBAE')

## pseudobulk comparison of a feature
af.niche2$Type
p_niche <- plot_pseudobulk_feature(
  seurat_obj = af.pseubluk,
  feature = "CXCL12",
  assay = "Spatial",                # specify the SCT assay
  sample_col = "group",
  type_col = "Type",
  comparisons_list = c(my_comparisons1,my_comparisons2,my_comparisons3),
  colors = colors_Type
)
p_niche

pdf('niche1 Epi EDNRB  .pdf',height = 4,width = 3)
p_niche
dev.off()

####################
#correlation between genes and metadata, pseudobulk
#define your preferred SCI colors; names must match af.cellpro$Type exactly
#my_colors <- c("Young" = "#4DBBD5FF", "Middle" = "#00A087FF", "Elderly" = "#E64B35FF")

all_types <- unique(af.cellpro$Type)
uni_color <- rep("#4682B4", length(all_types))
names(uni_color) <- all_types

ks_pseudobulk_cor(
  seurat_obj = af.cellpro,
  var1 = "FOXE1",                                  # gene name
  var2 = "REACTOME_CELLULAR_SENESCENCE_score",      # metadata column name
  group_by = "group",                              # sample ID column
  type_by = "Type",                                # experimental group column
  cor_method = "pearson",                          # options: "pearson" or "spearman"
  point_palette = uni_color                        # pass in your custom colors
)


af.ST$REACTOME_CELLULAR_SENESCENCE_score
#####compute differences directly at single-cell level
# define a large list of pairwise comparisons

af.niche<- subset(af.ST, niche %in% c('niche_4','niche_6'))


#all_comparisons <- list(c("Epi_1", "Epi_4"))

p_multi <- plot_sc_violin_contrast(
  seurat_obj = af.niche,
  feature = "factor1",
  group_col = "Type",
  comparisons_list = c(my_comparisons1,my_comparisons2,my_comparisons3), # pass in multiple group combinations
  colors = colors_Type
)

print(p_multi)
ggsave("niche_factor1_score.pdf",width = 6,height = 5)


##########################
##view the distribution of values for a feature
ks_plot_density_summary(af.ST_Fib1, feature = "Plasma_cells", group_by = "Type",
                        custom_colors=colors_Type,x_limits = c(0,0.2)
)

af.ST_Fib1$Plasma_cells





##############################
###compute regional cell convex hull area

# run the function
my_comparisons <- list(
  c("A", "M"), 
  c("A", "Y"), 
  c("M", "Y")
)
# 2. set your preferred colors (count must be >= number of groups)
my_colors <- c("A" = "#F3AE63", "M" = "#73558B", "Y" = "#93CBAE") 
# 3. run the function
density_res <- CompareRegionalDensity(
  seu_obj = af.ST,
  celltype_col = "Cell_Anno",      
  target_celltype = "Fibroblasts",  
  region_col = "niche",            
  target_region = c("niche_1","niche_3",'niche_5'),       
  sample_col = "group",            
  condition_col = "Type",
  comparisons = my_comparisons,    # 👈 pass in the pairwise comparison list
  plot_colors = my_colors,         # 👈 pass in custom colors
  test_method = "t.test"           # 👈 default t-test; for skewed data change to "wilcox.test"
)

# view the computed statistics table (area, count, density per sample)
head(density_res$Data)

# show the generated differential comparison plot
density_res$Plot
ggsave(filename = "fibroblast_density_diff.png",  width  = 3.5,  height   = 5,  units = "in",dpi= 300)
