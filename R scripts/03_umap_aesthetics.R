# =============================================================================
# 03_umap_aesthetics.R
# -----------------------------------------------------------------------------
# UMAP aesthetics: contour arrow plots / misty UMAP / target cell population masking / Circlize chord plot
#
# Features
#   - Aesthetic UMAP plots with contour lines (unchull) and small arrows
#   - Misty (smooth) UMAP, colored by niche
#   - generateMask to create a mask that highlights only the target cell population region
#   - Circlize chord plot from plot1cell, showing clusters and metadata tracks
#
# Main inputs
#   - Seurat objects (af / af.ST, with UMAP reduction and Celltype, niche, Type metadata)
#
# Main outputs
#   - LEC C1_C3 UMAP.pdf, ST niche UMAP.png, circlize_plot.png
#
# Dependencies
#   In-repo custom functions: umap_with_unchull.R
#   Main R packages: Seurat, plot1cell, circlize, mascarade, ggplot2
#
# Notes
#   This script was compiled from analysis drafts, retaining the original analysis logic without changing algorithms or workflow.
#   Before running, modify the input / output paths appearing in the script as needed (the original drafts mostly used local absolute paths),
#   and replace the sample names, cell type names, comparison groups, etc. in the examples with your own data.
#   It is recommended to first set the repository root as the working directory (open
#   sc-spatial-analysis-scripts.Rproj in RStudio) so that source("./R/...") resolves correctly.
# =============================================================================

#http://www.360doc.com/content/24/0430/20/82995563_1122020662.shtml


#############################
source("./R/umap_with_unchull.R")
table(af.ST$Cell_Anno,useNA = "always")
af.ST$Cell_Anno<-af.ST$Celltype
af.ST$Cell_Anno <- factor(af.ST$Cell_Anno, 
                      levels = c('EC','T cell','LAMs','Res Mac','Inf Mac',
                                 'DCs' , 'Mono' ,'Cycling Myeloid', 'SC' , 'vSMCs',
                                 'FBs C1' ,'FBs C2' ,'FBs C3','FBs C4','FBs C5',
                                 'FBs C6', 'FBs C7', 'FBs C8'
                                 ))
colors_prccLEC <- c(  "#CB95BB","#83C5BD",'#B5972D',"#3F6F76",
                   '#5d2b84',
                 '#d5a339','#ddb58b', "#F19294", '#c65f2d',"#D51F26",
                 "#edc951",'#7db281',
                 "#435B95","#8A9FD1",  "#4991C1",
                 "#C6307C","#7A989A","#4f372d"    )
#'#8e75d7', '#d5a339',

###Aesthetic plot with small arrows####
#plotData$UMAP_1
p <- plot_umap_with_unchull(
  seurat_obj = af,
  cluster_col = "Celltype",
  colors = colors3,
  point_size = 0.03#0.2-0.3
)

print(p)

#pdf(file = "spatial_major_cluster_UMAP.pdf",width =9,height = 6)
pdf(file = "LEC C1_C3 UMAP.pdf",width =6,height = 4)
p
dev.off()

###Misty UMAP
p1 <- plotUMAP(
  seurat_obj = af.ST,
  reduction = "umap",
  cluster_col = "niche",
  color_palette = colors_niche,
  point_size = 0.2,
  axis_line_size = 0.4,      # axis line thickness
  legend_point_size = 3,     # legend point size
  legend_text_size = 5,      # legend text size
  legend_ncol = 1 ,          # number of legend columns
  legend_alpha = 0.7 
)
p1
ggsave(filename ="ST niche UMAP.png",height = 3,width = 4,plot = p1,dpi = 300,bg="white")





#Highlight only the target cell population


emb <- Embeddings(af,reduction = 'umap') %>% #extract coordinate data    
  as.data.frame() %>% mutate(celltype=af@meta.data$Celltype)
maskTable <- generateMask(dims=emb[,1:2],
                          clusters=emb$celltype,#######do not change this
                          smoothSigma = 0.015,  
                          minDensity = 0.1,      
                          gridSize = 200,      
                          type = 'partition')

colors<-c("grey","#FB8072")
plot6<-DimPlot(af ,reduction= 'umap', group.by = "Type",
               label =TRUE,pt.size = 1)#+ scale_color_manual(values=colors)+NoLegend()
  +geom_path(data=maskTable[cluster=='plaque'],######select a certain cell
                                                  aes(x=umap_1, y=umap_2,group=group), 
                                                  linetype=6,size=1)
 
  
plot6


af$Type

library(plot1cell)

human_data<-af
# ============1. Circlize plot to visualize cell clustering and meta data======
#####this may require a clean environment; loading other packages may cause errors
####Prepare data for ploting

circ_data <- prepare_circlize_data(human_data, scale = 0.8)

length(names(table(Idents(human_data))))

cluster_colors<-colors7[1:8]

(length(names(table(human_data$Type))))
rep_colors<-two_group_col3
#(length(names(table(human_data$oldType))))
#Typecolor<-colors2[1:27]
###plot and save figures
png(filename =  'circlize_plot.png', width = 8, height = 6,units = 'in', res = 500)
plot_circlize(circ_data, do.label = T, pt.size = 0.2,           
              col.use = cluster_colors ,bg.color = 'white',        
              kde2d.n = 500, repel = T, label.cex = 1.2)

add_track(circ_data, group = "Type",colors = rep_colors, track_num = 2) #add legend
#add_track(circ_data, group = "oldType", colors = Typecolor, track_num = 3)
legend("topright",   
       legend = unique(human_data$Type),   
       col = rep_colors,       pch = 15,   
       cex=0.5,       pt.cex=2,    
       box.lwd = 0,  
       bg=NULL,       y.intersp = 1,    
       box.lty=0)
 legend("topleft",   
        legend = unique(human_data$Celltype),   
        col = cluster_colors,       pch = 15,   
        cex=0.5,       pt.cex=2,    
        box.lwd = 0,  
        bg=NULL,       y.intersp = 1,    
        box.lty=0)


dev.off()

