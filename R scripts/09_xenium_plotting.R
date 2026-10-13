# =============================================================================
# 09_xenium_plotting.R
# -----------------------------------------------------------------------------
# Xenium spatial visualization: density plots / spatial scatter / cross-sample gene expression / CSV export
#
# Functions
#   - Export a cell-annotation CSV readable by Xenium Explorer
#   - RNA and protein (e.g. CD38) spatial density plots (Nebulosa plot_density)
#   - spatial DimPlot highlighting target cell populations
#   - Cross-sample spatial scatter of genes / continuous variables (e.g. vertical distance)
#   - Joint density plot of multiple markers
#
# Main inputs
#   - Xenium / CODEX merged Seurat object (with Protein assay)
#   - annotation CSV such as xenium_cellcharter_CN.csv
#
# Main outputs
#   - Celltype_Fibsub.csv, S2_A4_*_density.pdf, vertical distance.png, etc.
#
# Dependencies
#   In-repo custom functions: spatial_distance_analysis.R
#   Main R packages: Seurat, Nebulosa, ggplot2
#
# Notes
#   This script was organized from analysis drafts, keeping the original analysis logic intact with no changes to algorithms or workflow.
#   Before running, adjust the input / output paths as needed (the original drafts mostly used local absolute paths),
#   and replace the example sample names, cell type names, and comparison groups with your own data.
#   It is recommended to first set the repository root as the working directory (open
#   sc-spatial-analysis-scripts.Rproj in RStudio) so that source("./R/...") resolves correctly.
# =============================================================================


af.output<-  subset(af.stroma, CN == "CN_5")

# output a CSV file readable by Xenium Explorer
# create a data frame containing cell names and domain information
export_data <- data.frame(
  cell_id = colnames(af.ST),    # first column: cell name
  group = af.ST$Celltype_Epi_Fib  # second column: domain information
)

# write the data frame to a CSV file
write.csv(export_data, 
          file = "Celltype_Fibsub.csv", # specify file name
          row.names = FALSE)         # do not write row names (usually numeric indices) to the file

###########################################################################
source("./R/spatial_distance_analysis.R")
#######################


af.LEC <- add_metadata_from_csv(  ##the first column read just needs to be the cell name
  csv_file = "./xenium_cellcharter_CN.csv",  # CSV file containing cell information
  seurat_obj = af.LEC,         # input Seurat object
  metadata_cols = c("group")  # column names to add
)
table(af.ST$group)
af.ST_radius

AA<- qread("./slide1 data/output-XETG00295__0066285__XP-1__20250924__051644/data/Xenium_codex_merge.qs")
###########################################################################
#RNA plotting
###########################################################################
DefaultAssay(af.ST) <- "Protein"
rownames(af.ST)
table(af.ST$group)
af.sub<-  subset(af.ST, group == "slide2_A4")

pdf(file = "S2_A4_CD38_RNA_density.pdf",width =6,height = 5)

plot_density(afd, features = c("POSTN"),reduction = 'umap',slot='data',
             pal = 'magma',method='ks',raster = T, size = 0.8) &
  theme_blank()&#this theme comes from ggnetwork and removes all ggplot backgrounds
  theme(legend.frame = element_rect(colour = "black"),#modify legend border
        legend.ticks = element_line(colour = "black", linewidth  = 0),#modify legend ticks
        legend.key.width = unit(0.3, "cm"),#modify legend width
        legend.key.height = unit(0.8, "cm"),#modify legend height
        legend.title = element_text(color = 'black', face = "bold", size=8))#modify legend title
dev.off()

###########################################################################
#Protein plotting
############################################################################
af.ST2<-af.ST

DefaultAssay(af.ST2) <- "Protein"

af.ST2 <- NormalizeData(af.ST2, normalization.method = "LogNormalize", scale.factor = 10000)
table(af.ST2$group)
af.sub2<-  subset(af.ST2, group %in% c('slide2_A4' ) )
rownames(af.ST2)


pdf(file = "S2_A4_CD38-protein_density.pdf",width =6,height = 5)
plot_density(af.sub2, features = c("CD38-protein"),reduction = 'spatial',
             slot='data',method='wkde',
             pal = 'magma', raster = T, size = 0.8) &
  theme_blank()&#this theme comes from ggnetwork and removes all ggplot backgrounds
  theme(legend.frame = element_rect(colour = "black"),#modify legend border
        legend.ticks = element_line(colour = "black", linewidth  = 0),#modify legend ticks
        legend.key.width = unit(0.3, "cm"),#modify legend width
        legend.key.height = unit(0.8, "cm"),#modify legend height
        legend.title = element_text(color = 'black', face = "bold", size=8))#modify legend title
dev.off()


#######################


af.sub2$Cell_Anno

colors_prcc <- c('#5d2b84',"#89558D",
                 '#d5a339','#ddb58b', "#F19294", '#c65f2d',"#D51F26",
                 "#edc951",'#7db281',
                 "#435B95","#8A9FD1",  "#4991C1",
                 "#C6307C","#7A989A","#4f372d"    )

colors_dimplot<-c("#5d2b84","grey",
                  "grey","grey","grey","grey","grey",
                  "grey","grey",
                  "grey","grey","grey",
                  "grey","grey","grey"
)

pdf(file = "S2_A4_CD38-protein_density.pdf",width =6,height = 5)
DimPlot(af.sub2 ,reduction= 'spatial', group.by = "Cell_Anno",
        label =F,pt.size = 1,
        cols= colors_dimplot)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")

dev.off()


####multiple markers###########
p <- plot_density(af_combined, features = c("CD8A",'CD3E'), reduction = 'spatial', slot='count',
                  joint = TRUE, pal = 'magma', raster = TRUE, size = 0.8)&
  theme_blank()&#this theme comes from ggnetwork and removes all ggplot backgrounds
  theme(legend.frame = element_rect(colour = "black"),#modify legend border
        legend.ticks = element_line(colour = "black", linewidth  = 0),#modify legend ticks
        legend.key.width = unit(0.3, "cm"),#modify legend width
        legend.key.height = unit(0.8, "cm"),#modify legend height
        legend.title = element_text(color = 'black', face = "bold", size=8))
ggsave("spatial_density_plot.png", p, width = 15, height = 8, dpi = 500)
  



#####dimplot##########
table(af.ST$group)


af.niche <- subset(af.ST, niche %in% c('niche_5','niche_3'))
af.ST$vertical_distance
# 2. call the function in one step to add the centroid coordinates to the reduction layer!
ks_plot_gene_across_samples(
  seurat_obj = af.ST,
  feature    = "vertical_distance",
  group_col = "group",
  group_ids  = c("A6"),
  ncol       = 1,   # arrange in one row
  pt_size    = 0.1,
  colors=c('#313536', '#61BECD','#89BF69','#E1D432','#C3400C','#821407')

)

####################
 ggsave(   filename = "vertical distance .png",   width = 7, height   = 5, units    = "in", dpi  = 300 )
 

p<-DimPlot(af.treat_filter, reduction = "spatial", label = F, label.size = 2.5,
        group.by = "Signal_Status",
        pt.size = 0.01,cols= colors3)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
p <- p + coord_fixed(ratio = 1)
p
# save: PNG (commonly used for sharing / preview)
ggsave(
  filename = "A6 celltype.png",
  plot     = p,
  width    = 7,      # inches, change as needed
  height   = 7,      # spatial plots usually have width≈height, especially with coord_fixed
  units    = "in",
  dpi      = 300
)

plot_density(af.plot, features = c("FOXE1"),reduction = 'spatial',slot='data',
             pal = 'magma',method='ks',raster = T, size = 0.8) &
  theme_blank()&#this theme comes from ggnetwork and removes all ggplot backgrounds
  theme(legend.frame = element_rect(colour = "black"),#modify legend border
        legend.ticks = element_line(colour = "black", linewidth  = 0),#modify legend ticks
        legend.key.width = unit(0.3, "cm"),#modify legend width
        legend.key.height = unit(0.8, "cm"),#modify legend height
        legend.title = element_text(color = 'black', face = "bold", size=8))#modify legend title
