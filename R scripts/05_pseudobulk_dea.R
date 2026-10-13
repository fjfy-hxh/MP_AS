# =============================================================================
# 05_pseudobulk_dea.R
# -----------------------------------------------------------------------------
# Pseudobulk construction and differential expression analysis (DESeq2 / limma / edgeR)
#
# Features
#   - Aggregate by sample and cell type to build pseudobulk expression profiles
#   - run_pseudobulk_dea for between-group differential analysis (DESeq2, limma, edgeR)
#   - Volcano plot visualization of differential results
#
# Main inputs
#   - Seurat object (af.ST, with metadata such as niche, Celltype / Type, sample groups)
#   - Sample-to-group mapping (specified manually via mapping_df in this script)
#
# Main outputs
#   - pseudobulk_xxx_MP_vs_CON.pdf, differential gene result tables
#
# Dependencies
#   Repo-internal custom functions: pseudobulk_functions.R
#   Main R packages: DESeq2 / limma / edgeR (called inside the functions), Seurat, ggplot2
#
# Notes
#   Organized from analysis drafts; the original logic, algorithm and workflow are unchanged.
#   Before running, adjust the input / output paths as needed (the original drafts were
#   mostly local absolute paths), and replace the example sample names, cell type names
#   and comparison groups with your own data. Run from the repository root (open
#   sc-spatial-analysis-scripts.Rproj in RStudio) so source("./R/...") resolves correctly.
# =============================================================================


source(file = "./R/pseudobulk_functions.R")

#af.epi1<- subset(af.ST, niche %in% c("niche_2",'niche_3','niche_1'))
table(af.peusobulk$Celltype)
af.peusobulk<- subset(af.ST, niche %in% c("niche2"))
af.peusobulk$oldType
#table(af$Cell_Anno)
pb_res <- create_pseudobulk(
  seu = af.peusobulk,
  group.by = "oldType",
  cell.type.col = "Celltype",
  cell.type = "SC",
  assay = "RNA",
  layer = "counts",
  min.total = 500,
  min.cells = 3
)
table(af.peusobulk$Type)
sample_names <- pb_res$meta$group  # e.g. here it is c("A11", "A2", "M1", "Y10" ...)
sample_names
#group_names <- substr(sample_names, 1, 1) # extract the first letter: "A", "A", "M", "Y"
group_names <- c(rep("CON", 4), rep("MP", 5))  # establish the grouping
group_names

mapping_df <- data.frame(
  sample = sample_names,
  condition = group_names
)

res_data <- run_pseudobulk_dea(
  pb_output = pb_res, 
  sample_mapping = mapping_df, 
  method="DESeq2",  #"limma", "DESeq2", "edgeR"
  contrast_str = "MP - CON"##########
)
p2 <- plot_volcano(deg_result = res_data,
                   logFC_cutoff = 0.5, pval_cutoff = 0.05,
                   title_prefix = "MP - CON: ")
res_data
table(af.ST$group )
######### save figure #####

pdf('pseudobulk_FBs_C1_MP_vs_CON.pdf',height = 3,width = 4)
p2
dev.off()
