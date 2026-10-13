# =============================================================================
# 01_sc_main_analysis.R
# -----------------------------------------------------------------------------
# Single-cell main pipeline: QC / normalization / clustering / manual annotation / metabolism and drug sensitivity
#
# Features
#   - QC visualization (violin plot, metric correlation scatter), filter mitochondrial and ribosomal genes
#   - LogNormalize normalization, variable features, ScaleData, PCA, harmony batch correction
#   - FindNeighbors / FindClusters clustering, FindAllMarkers to identify marker genes per cluster
#   - UMAP visualization, manual cell type annotation by cluster
#   - Cell proportion statistics between groups and differential boxplot (based on custom ratio function)
#   - scMetabolism single-cell metabolic pathway analysis, scRank drug sensitivity analysis
#   - cNMF result enrichment, gene correlation scatter, marker gene dotplot / average expression heatmap
#
# Main inputs
#   - af_MP_scRNA.qs / af_MP_FB_subclusters.qs / af_MP_trajectory.qs / MP_spatial_overall.qs (Seurat objects, prepare them yourself)
#   - gene_spectra_score file output by cNMF (optional, used for the cNMF enrichment part)
#
# Main outputs
#   - 01.vlnplot.pdf, 01.corqc.pdf, 05-cluster.UMAP.pdf, anno.UMAP.pdf
#   - FMC cluster_markers.csv, metabolic pathway / drug sensitivity / gene correlation related plots
#
# Dependencies
#   In-repo custom functions: ks_scAverExp.R, ks_scDotplot.R, Singlecellratio_plotstat.R, ks_cellration_plot.R, Mouse_sc_metabolism.R, scRank_functions.R
#   Main R packages: Seurat, harmony, tidyverse, patchwork, Matrix, SCpubr, scMetabolism, scRank, CellChat, Nebulosa, scRNAtoolVis, scCustomize, qs, etc.
#
# Notes
#   This script was compiled from analysis drafts, retaining the original analysis logic without changing algorithms or workflow.
#   Before running, modify the input / output paths appearing in the script as needed (the original drafts mostly used local absolute paths),
#   and replace the sample names, cell type names, comparison groups, etc. in the examples with your own data.
#   It is recommended to first set the repository root as the working directory (open
#   sc-spatial-analysis-scripts.Rproj in RStudio) so that source("./R/...") resolves correctly.
# =============================================================================

rm(list = ls())
options(stringsAsFactors = F)
gc()
#  The custom function library is included in the ./R directory of this repository (see README)
af
af[["RNA"]] <- as(af[["RNA"]], "Assay")#5 to 4

af[["RNA"]] <- as(object = af[["RNA"]], Class = "Assay5")#4 to 5


colors_type<-c("royalblue",  "indianred")

nmf_color<-c('#9C0642','#5FB6A9','#614DA1','#CD3B4E','#EA6547',
             '#F79C5A','#FACF7E','#FCF2A8','#F1FBAE','#D0ED9D','#97D5A5'
             )

colors_all<-c('#B7B663','#9d9ac6',"#BBE379",'#960937',"#0F7240","#edc951")

colors_niche<-c('#FD8352', '#282568','#73c9e4',  '#156434')


colors_Celltype<-c('#d51f26',
                   '#0F7240',#T
              '#FCCC3A','#d5a339','#f3758a','#C67052','#F49E47','#C0937E',#Mac
              "#5d2b84", #SC 
              "#9d9ac6", #SMC
              "#b1ddcc",  "#00a7fa","#92B2D8",'#C9E1AF',
              "#bbe379", '#b9e6fb', '#7AA15E',"#31aaaa" #FMC
           )
colors_FB<-c( "#b1ddcc",  "#00a7fa","#92B2D8",'#C9E1AF',
              "#bbe379", '#b9e6fb', '#7AA15E',"#31aaaa" #FMC
)
colors_FBsub<-c("#00AFBB", "#E7B800", "#FC4E07")

library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(patchwork)
library(ggplot2)
library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(ggpubr)
library(Biobase)
library(Matrix)
library(RcppML) 
library(magick)
library(GeneNMF)#R>=4.3
library(scRNAtoolVis)
library(Nebulosa)
library(reshape2)
library(ggnetwork)
library(ggimage)
library(mascarade)
library(UCell)
library(ggrastr)
library(qs)


#qsave(af.ST, file = "MP_spatial_overall.qs")

#qsave(afd, file = "af_MP_trajectory.qs")

setwd("./")
getwd()



af.sc<- qread("af_MP_scRNA.qs")
afd<- qread("af_MP_trajectory.qs")# af_MP_scRNA_FBs.qs  af_MP_FB_subclusters.qs
af.FB<- qread("af_MP_FB_subclusters.qs")# af_MP_scRNA_FBs.qs  af_MP_FB_subclusters.qs

af.ST<-qread("MP_spatial_overall.qs") 

table(af$Type,af$Celltype)
table(af.ST$niche,af.ST$Cell_Anno,af.ST$Type )
##############################################
###################################Data normalization

# Visualize QC metrics with violin plots, ncol is the number of columns in the plot layout
pdf(file = "01.vlnplot.pdf",width = 8,height = 5)
VlnPlot(af.sc, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
dev.off()

# Correlation between metrics
plot1 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.mt")+ RotatedAxis()
plot2 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "percent.rb")+ RotatedAxis()
plot3 <- FeatureScatter(af, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ RotatedAxis()
#Combine plots
pdf(file = "01.corqc.pdf",width =12,height = 5)
plot1+plot2+plot3+plot_layout(ncol = 3)      #plot_layout, a patchwork function, specifies how many plots per row
dev.off()

plot3
# Remove mitochondrial and ribosomal genes
exp_m <- GetAssayData(af, layer = 'counts')
gene.all <- rownames(exp_m)
mt.gene  = grep('^MT-', gene.all, value = TRUE)#mitochondrial genes
ribosome.gene = grep('^RPL|^RPS|^MRPS|^MRPL', gene.all, value = TRUE)#ribosomal genes
features=setdiff(gene.all,c(mt.gene,ribosome.gene))#remove unwanted genes
filtered_counts <- GetAssayData(af, layer = 'counts')[features, ]# get the filtered gene matrix
af <- SetAssayData(af, layer = 'counts', new.data = filtered_counts)# update the filtered gene matrix into the af object


table(af$Celltype)
######################
#harmony works better after SCT then PCT  
############
##Normalization using the LogNormalize method
af <- NormalizeData(af, normalization.method = "LogNormalize", scale.factor = 10000)
options(future.globals.maxSize = 8000 * 1024^2)
#af <- SCTransform(af, assay = "RNA", verbose = FALSE)
DefaultAssay(af) <- "RNA"
af[["RNA"]] <- as(af[["RNA"]], "Assay")#5 to 4

af <- FindVariableFeatures(af, selection.method = "vst", nfeatures = 2000)
af <- ScaleData(af)

af <- Seurat::RunPCA(af)

pct <- af[["pca"]]@stdev / sum(af[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
ElbowPlot(af)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")

#nmf non-negative matrix factorization dimensionality reduction
#af <- runNMF(af, k = 10, assay="SCT")
#af@reductions$harmony <- af@reductions$NMF
af <- Seurat::RunPCA(af,npcs = 14)
#####################################(optional harmony batch correction and dimensionality reduction)################################################################################################
af <- RunHarmony(af, group.by.vars = "Type")

rownames(af.ST)
#####Merge samples!!!!#######################################################################
namecon="Con"     #define the control group name
nametreat="Mp"   #define the experimental group name
#First take a quick look at the current sample grouping
table(af$Type)

##Start merging
#Determine the control group; samples not in it are automatically set as the experimental group
conaf=c("Con_1","Con_2",'Con_3','Con_4')
#Modify cell grouping
af$oldType=af$Type  #save the old grouping
af$Type=factor(ifelse(af$Type %in% conaf,namecon,nametreat))


Idents(af)=af$Type

table(af$oldType)


######################################################
#Create a stage vector with the same length as the metadata, initialized to NA
table(afe$leiden)
stage_vector <- rep(NA, length(af.ST$niche))

# Assign values by condition (using the Seurat object metadata index)
af.ST
stage_vector[af.ST$niche %in% c("niche1","niche2")] <- 'lesion'
stage_vector[af.ST$niche %in% c("niche4","niche3")] <- 'other'

# Add the stage vector to the Seurat object metadata
af.ST$lesion <- stage_vector

# Verify the result
table(af.ST$lesion)


#########################################################################################################################################
# # Visualize feature genes of the first two PCs

# First build a KNN graph based on Euclidean distance in PCA space
af <- FindNeighbors(af, dims = 1:14,reduction = "pca")

#Set different resolutions to observe clustering; dim is the number of PCs selected
#af=FindNeighbors(af, reduction = "harmony")
#for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1,1.2,1.5,2,2.5,3)) {
#  af=FindClusters(af, graph.name = "RNA_snn", resolution = res, algorithm = 1)}
#apply(af@meta.data[,grep("RNA_snn_res",colnames(af@meta.data))],2,table)

# Cluster and optimize
# resolution parameter: the larger the value, the more clusters; choose based on the above
# 0.4-1.2 typically returns good results for single-cell datasets of around 3K cells
# Optimal resolution often increases for larger datasets. 


#Choose a resolution for dimensionality reduction

af <- FindClusters(af, resolution = 0.4)

# Check cluster IDs
head(Idents(af), 5)

# Check how many cells per cluster
head(af@meta.data)
table(af@meta.data$seurat_clusters)
# Identify marker genes for each cluster; only.pos: keep only upregulated DEGs

af.markers <- FindAllMarkers(af.ST, only.pos = T, min.pct = 0.25, 
                             logfc.threshold = 0.1)
write.csv(af.markers,file = "FMC cluster_markers.csv")
af.markers_cluster <- af.markers[af.markers$cluster == 8, ]

myMarkers=c("COL4A1",'TIMP2','APP','CCN1','PECAM1','SYNM','SPTBN1')
pdf(file = "05_marker_genes_per_cluster.pdf",width =4.5,height = 4)
jjVolcano(diffData = af.markers, aesCol = c('#42AFFC','tomato'),topGeneN=2,
          #myMarkers=myMarkers,
          flip = F,legend.position = c(0.9,0.2),#polar=T,
          tile.col = c('#9d9ac6','#00a7fa','#354c85','#31aaaa'))##number of clusters
dev.off()

## Visualize cells in low-dimensional space UMAP/tSNE
af <- RunUMAP(af, dims = 1:14, reduction = "pca",min.dist = 0.7)

# Visualize UMAP/tSNE
pdf(file = "05-cluster.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = F, label.size = 2.5,pt.size = 1,
       group.by = "seurat_clusters",
        cols= colors2)+#ggsci::scale_color_igv()+
  theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()

af$seurat_clusters<-Idents(af)

####################

pdf(file = "08.ann_cluster_marker.pdf",width =20,height = 7)
do_DotPlot(sample = af,features = genes,dot.scale = 12,
           legend.length = 50,
           legend.framewidth = 2, font.size =12)
dev.off()

#manual annotation

af <- FindClusters(af, resolution = 0.4)
table(af@active.ident)
afe$`MP1`
ann.ids <- c("vSMCs",  #0Injured Proximal tubular cell  not a cluster, it is my numbering, cluster+1
             "vSMCs",  #1Injured Proximal tubular cell
             "FBs like",  #2
             "FBs like", #3
             "Mono/Mac", #4
             "FBs like",#5
             "EC",#6
             "vSMCs",#7
             "FBs like",#8
             "T cell",#9
             "remove",#10
             "FBs like",#11
             "SC",#12
             "vSMCs"#13
)


Idents(af)

afidens=mapvalues(Idents(af), from = levels(Idents(af)), to = ann.ids)
Idents(af)=afidens

af$Celltype=Idents(af)
#########Visualize results after manual annotation
pdf(file = "anno.UMAP.pdf",width =7,height = 5.5)
DimPlot(af, reduction = "umap", label = F, label.size = 2.5,pt.size = 0.0001,
        group.by = "Celltype",
        cols= colors_Celltype)+#ggsci::scale_color_igv()+
  theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()



table(af$Celltype,af$Type)
####################

colors<-c("#FB8072","gray")

pdf(file = "controlneu.UMAP.pdf",width =7,height = 5.5)
DimPlot(af.neu ,reduction="umap", group.by = "Type",
        label =F)+
  scale_color_manual(values=colors)+
  NoLegend()
dev.off()

source(file = "./R/ks_scAverExp.R")
source(file = "./R/ks_scDotplot.R")

afgenes=c('Pecam1','Cd3d','Clec10a','Trem2','Fabp5',
          'Tnf','Il1b','Mrc1','Stab1',
          'Cd68','Cd14','Ly6c2','Mki67',
          'Acta2','Myh11','Dcn','Lum','Mmp14',
          'Vcam1','Lgals3','Icam1','Col1a1','Sox10','Mpz'
)

table(af$Celltype)
afgenes=c('Acta2','Myh11','Col1a1','Dcn','Cd8a',
          'Nkg7','Cd79a','Cd68','Cd14','Sox10','Mpz','Krt19')
featureSets <- list("EC" = c('Pecam1'),#, "Mzb1"
                    'T cell'=c('Cd3d'),
                    'LAMs'=c( 'Trem2','Fabp5'),
                    'Res Mac'=c('Stab1','Mrc1'),
                    "Inf Mac"=c('Tnf','Il1b'),#Smooth muscle cell
                    'DCs'=c('Cd1c','Clec10a'),
                    'Mono'=c('Cd68','Cd14','S100a8','Ly6c2'),
                    'Cycling Myeloid'=c('Mki67'),
                    'SC'=c('Sox10','Mpz'),
                    "vSMCs"= c('Acta2','Myh11'),#Macrophage
                    "FBs C1" = c('Dcn'),#neutrophil
                    "FBs C2" =c('Lum'),
                    "FBs C3" = c('Mmp14'),
                    "FBs C4" = c('Vcam1'),
                    "FBs C5" = c('Cxcl12'),
                    "FBs C6" = c('Lgals3'),
                    "FBs C7" = c('Icam1'),
                    "FBs C8" = c('Cx3cl1')#nature killer
)

sce<-af.ST
sce$celltype<-sce$Celltype
Exp_scRNA <- ks_scAverExp(obj = sce,#seurat or ArchR proj
                          features = featureSets,
                          CellTypes = unique(sce$celltype))

scRNA_dotplot = ks_scDotplot(obj = sce,#seurat or ArchR proj
                             features = featureSets,
                             CellTypes = unique(sce$celltype),
                             avgPctMat= Exp_scRNA,
                             pal = colors_Celltype )

#########Plot marker gene heatmap################

# only.pos: keep only upregulated DEGs
af.markers <- FindAllMarkers(af.ST, only.pos = T, min.pct = 0.25, logfc.threshold = 0.75)

# get top 10 genes
topaf.markers <- af.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)


# plot
pdf(file = "niche heatmap.pdf",width =7,height = 10)
AverageHeatmap(object=af.ST,
               markerGene = topaf.markers$gene,slot = "data",assays='Spatial',
               clusterAnnoName =T,annoCol=T,myanCol= colors_niche,# htRange=c(-0.5, 0, 2),
               htCol =c('#435B95','white','#C3400C')
               )



afgenes=c('Spp1','Lgals3','Vcam1','Icam1','Thbs1','Col1a1','Col3a1','Fn1',
          'Postn','Mmp2','Mmp14','Mmp3',
          'Cx3cl1','Cxcl12','Mif')

afgenes=c('Ccl7','Cxcl12','Ackr3','Ccl19','Ccl2','Sdf1' )
afgenes=c('Ccr2','Cxcr4','Cxcr3','Ptn','Ccn1' )

afgenes=c('Tnf','Il1b','Il6','Spp1','Ly6C','Fcgr4','Itgam','S100a9','S100a8',
          'Trem2','Apoe','Fabp5','Lyve1','Folr2','Nampt')
afgenes=c('Pecam1','Vwf','Cd3d','Clec10a','Trem2','Fabp5',
          'Tnf','Il1b','Mrc1','Stab1',
          'Cd68','Cd14','Ly6c2','Mki67',
          'Acta2','Myh11','Dcn','Lum','Mmp14',
          'Vcam1','Lgals3','Icam1','Col1a1','Sox10','Mpz'
)
afgenes=c('Ccn2','Cebpd','Tmsb4x','Igfbp4')

table(Idents(af.ST))
Idents(af.ST)<-af.ST$Celltype
Idents(afd)<-afd$Anno_type
Idents(af.sc)<-af.sc$Celltype
pdf(file = "chemokine_cross_cell_expression.pdf",width =6,height = 3.5)
DotPlot(afd , features = afgenes )+theme_bw()+#remove background, rotate plot
  
  theme(panel.grid = element_blank(),
        
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+#display text at 90 degrees
  
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#506AAA',"#728DC1",'#91BBDF',"#B6DBF1"
                                    ,"#F6D3D9","#E79FBF","#E17692" ))+#color gradient settings
  # c('#64A2CC',"white",'#E07B56',"#731929")
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+coord_flip()
dev.off()
  ###################

################################
Idents(af)<-af$Celltype


#################
#Gene UMAP distribution plot
afgenes=c('Tnf','Il1b','Il6','Spp1','Ly6C','Fcgr4','Itgam','S100a9','S100a8',
          'Trem2','Apoe','Fabp5','Lyve1','Folr2','Nampt')
afgenes=c('Icam1','Vcam1','Cdh5',
          'Klf2','Klf4','Nos3','Pecam1','Vwf','Acta2','Vim','Tagln',
          'Apoe')
afgenes=c('Dcn','Col1a1','Acta2','Myh11','Pdgfra','Il6','Cxcl12')
afgenes=c('Acta2','Myh11','Tagln','Lum','Col1a1','Cd68')




input.pathway <- c("Glycolysis / Gluconeogenesis",
                   "Citrate cycle (TCA cycle)",
                   "Pyruvate metabolism",
                   "Oxidative phosphorylation",
                   'Nicotinate and nicotinamide metabolism')



#afgenes=read.table("stroke_geneFA.txt",header = F,sep = "\t")[,1]
pdf(file = "major_cluster_Cxcl12.pdf",width =6,height = 5)
plot_density(af.sc, features = c("Vcam1"),reduction = "umap",
             pal = 'magma', raster = T, size = 0.8) &
  theme_blank()&#this is a theme from ggnetwork, it can remove all ggplot background
  theme(legend.frame = element_rect(colour = "black"),#modify the legend frame
        legend.ticks = element_line(colour = "black", linewidth  = 0),#modify legend ticks
        legend.key.width = unit(0.3, "cm"),#modify legend width
        legend.key.height = unit(0.8, "cm"),#modify legend height
        legend.title = element_text(color = 'black', face = "bold", size=8))#modify legend title
dev.off()

table(af$cell)


##### Cell proportion overview
source("./R/Singlecellratio_plotstat.R")
table(af.ST$Cell_Anno)

af<- subset(af.ST, niche %in% c("niche1"))

af<-af.ST

af$celltype<-af$Cell_Anno
af$orig.ident<-af$Type########modify the meaning of orig.ident here


Singlecellratio_plotstat(af, group_by = "Type",###grouping
                         meta.include = c("Type","orig.ident"),#the first is sample grouping, the second is each sample within the group; only orig.ident can be filled in
                         color_by = 'cell.type'###this one is fixed; change the metadata to be analyzed to celltype
)+scale_fill_manual(values=colors_niche )#colors_prcc  LEC_color

ggsave("human_carotid_annotation_ratio.pdf",width = 5,height = 4.5)

######Second method###############################
source(file = "./R/ks_cellration_plot.R")
library(jjAnno)


ks_scFraction(group = af.ST$niche,
              celltype = af.ST$Cell_Anno,pal = colors_Celltype ,#CN_color4<-c("#c79000", "#93CBAE" ,"#00b1da" ,"#cdb6ff")
              colorLab =  colors_niche, Total = T)#colors5 group  colors_prcc
ggsave("ST niche- ration.pdf",width = 5.5,height = 3.5)

afd<-af.plaque
###########Cell differential boxplot###########
#af.ST<- subset(af.ST, niche %in% c(  'niche2' ))

af.ST$orig.ident<-af.ST$group###here orig.ident should be changed to each sample
af.ST$celltype<-af.ST$Cell_Anno  ####change celltype to the cells to be compared
my_comparisons<-list(c("CON","MP"))
Singlecellratio_plotstat(af.ST, group_by = "Type",
                         custom_fill_colors = c("royalblue",  "indianred"),
                         meta.include = c("Type","orig.ident"),
                         comparisons = my_comparisons, color_by = 'Type',
                         group_by.point = "orig.ident",label.x = 1, pt.size = 3,
                         label = 'p.format', ncol =4)
ggsave("ST niche4 cellratio_diff.pdf",width = 10,height =8)
# Cellratio <- prop.table(table(af$Celltype, af$Type), margin = 2) #calculate the proportion of different cell groups in each sample group
# Cellratio <- as.data.frame(Cellratio)
# colnames(Cellratio)[1]="Celltype"
# colourCount = length(unique(Cellratio$Celltype))
# colaa=distinctColorPalette(100)
# ggplot(Cellratio) + 
#   geom_bar(aes(x =Var2, y= Freq, fill = Celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
#   theme_classic() +
#   labs(x='Type',y = 'Ratio')+
#   coord_flip()+
#   theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")+   # legend: "left" left, "right" right, "bottom" bottom, "top" top
#   scale_fill_manual(values=sc_col)


#####################Cell subpopulation proportion for a specific gene
geneselect <- "Cxcl12"; cellselect <- "FBs C1"

# 2. Extract gene expression of this cell group
expr <- FetchData(af.ST, vars = geneselect, cells = colnames(af.ST)[af.ST$Celltype == cellselect])[[1]]
# 3. Compute the positive median and assign (no expression and below median -> low, the rest -> high)
med <- median(expr[expr > 0], na.rm = TRUE)
tag <- ifelse(!is.na(med) & expr >= med, " high ", " low ")

af.ST$genetype <- as.character(af.ST$Celltype)
af.ST$genetype[af.ST$Celltype == cellselect] <- paste0(geneselect, tag, cellselect)
table(af.ST$genetype)



afd<- subset(af.ST, niche %in% c('niche4','niche3'))




####scRANK drug sensitivity##############
library(scRank)

drug<-(utile_database)
View(drug$Drug_Info)

obj <- my_CreateScRank(input = afd,
                    species = 'mouse',#human
                    cell_type = 'sub.cluster',###cell classification
                    target = 'Irf1'
                   )
  
obj <- Constr_net(obj,n.core_cp = 5,
                  n.core = 1)
obj <- rank_celltype(obj)
obj
obj@para$target

obj <- init_mod(obj)

obj@cell_type_rank

plot_drug_function(
  obj,
  celltype = 'sensitive',
  category = "H",
  show_leading_edge = T
)
scRank_GSEA(obj, celltype = 'sensitive', disease_gene = 'Brd4')

plot_net(obj, mode = "heatmap", celltype = "FBs")
plot_net(obj, mode = "network", celltype = "sensitive")
plot_dim(obj)



plot_net(obj, mode = "heatmap", celltype = "sensitive")



#Mouse single-cell metabolic pathways
library(stringr)
library(Seurat)
library(rsvd)
library(qs)
library(scMetabolism)
library(BiocParallel)
register(MulticoreParam(workers = 4, progressbar = TRUE))
library(Seurat)
library(nichenetr)
library(dplyr)
source("./R/Mouse_sc_metabolism.R")
#Requires V4 data
res <- Mouse.sc.metabolism(af, metabolism.type = 'KEGG')
#########
#Single-cell metabolic pathway analysis; the V4 version is required
################################
#########################
afd

af[["RNA"]] <- as(af[["RNA"]], "Assay")

sub_scRNA <- af

res <-sc.metabolism.Seurat(obj = af,
                           method = "AUCell", # VISION, AUCell, ssgsea and gsva
                           imputation =F, ncores = 12, 
                           metabolism.type = "KEGG") # KEGG and REACTOME
countexp.Seurat <- res

#extract the score result
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
#convert the dots in the score barcode to underscores
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))  #non-standard gene IDs will error; replace underscore with -
#confirm the cell barcodes match
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
#[1] TRUE
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change) )
colnames(countexp.Seurat@meta.data)
#the relevant Seurat functions can be used directly
p1 <- FeaturePlot(countexp.Seurat,features = "Pyruvate metabolism")
p2 <- VlnPlot(countexp.Seurat,features = "Pyruvate metabolism")
p1 + p2
p1

pathways <- res@assays$METABOLISM$score
head(rownames(pathways))
# [1] "Terpenoid backbone biosynthesis"                
# [2] "Caffeine metabolism"                            
# [3] "Neomycin, kanamycin and gentamicin biosynthesis"
# [4] "Metabolism of xenobiotics by cytochrome P450"   
# [5] "Drug metabolism - cytochrome P450"              
# [6] "Drug metabolism - other enzymes" 

# Dimplot plotting
DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = "Citrate cycle (TCA cycle)", 
                   dimention.reduction.type = "umap",  
                   dimention.reduction.run = F, size = 1)

# DotPlot


rownames(pathways)
########Extract pathways of interest########
input.pathway <- c("Glycolysis / Gluconeogenesis",
                   "Citrate cycle (TCA cycle)",
                   "Pyruvate metabolism",
                   "Oxidative phosphorylation",
                   'Nicotinate and nicotinamide metabolism')


#input.pathway<-head(rownames(pathways))

pdf(file = "sc_metabolic_pathway_KEGG_tumor_subgroup.pdf",width =6,height = 4.5)
DotPlot.metabolism(obj = res,
                   pathway = input.pathway,
                   phenotype = "Celltype", #this parameter should be modified as needed
                   norm = "y")
dev.off()
#BoxPlot
# Since the developer defaults obj to countexp.Seurat, it needs to be renamed again

pdf(file = "sc_metabolic_pathway_boxplot1.pdf",width =8,height = 6)
BoxPlot.metabolism(obj = countexp.Seurat,
                   pathway = input.pathway, 
                   phenotype = "Celltype", #this parameter should be modified as needed
                   ncol =2)+scale_fill_manual(values =colors2)
dev.off()


genes <- list(c("PTEN","PIK3CA","KRAS","ARID1A","RCA1","WNT5A"))
names(genes) <- 'gene'
afe<-afd
#afe <- AddModuleScore_UCell(uterus,features=genes,name="_score")

#extract the data
library(ggpubr)



df1 <- data.frame(
  Cell_ID = countexp.Seurat$orig.ident,  # cell name
  Gene_Expression = countexp.Seurat$`Pyruvate metabolism`, # Nos2 expression value
  Cell_Annotation = countexp.Seurat$predicted.id  # cell annotation
)


#set the comparison groups
my_comparisons1 <- list(c("FMC_2", "vSMC")) 
my_comparisons2 <- list(c("FMC_2", "FMC_4"))
my_comparisons3 <- list(c("FMC_2", "FMC_7"))
my_comparisons4 <- list(c("vSMC", "FMC_4"))
colors<-colors_meta
#make a violin plot
ggplot(df1,aes(x=Cell_Annotation,y=Gene_Expression,fill=Cell_Annotation))+ 
  geom_violin(color='black',size=1)+#violin  
  theme_classic() +   theme(text = element_text(size=10, colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15), 
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(color = 'black', size = 12),
        axis.line = element_line(size = 1))+   
  theme(legend.position="none") +  scale_fill_manual(values = colors)+ 
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+ #boxplot 
stat_compare_means(method="wilcox.test",hide.ns = F,  
                   comparisons =c(my_comparisons1,my_comparisons2,my_comparisons3,my_comparisons4), 
                   label="p.signif",                
                   bracket.size=0.8,   
                   size=6)#add significance comparison
ggsave(filename =" Glycolysis / Gluconeogenesis metabolic pathway violin.png",height = 5.2,width = 5.5,dpi = 300,bg="white")

input.pathway <- c("Glycolysis / Gluconeogenesis",
                   "Citrate cycle (TCA cycle)",
                   "Pyruvate metabolism",
                   "Oxidative phosphorylation",
                   'Nicotinate and nicotinamide metabolism')



####scRANK drug sensitivity##############
library(scRank)
source(file = "./R/scRank_functions.R")

drug<-(utile_database)
View(drug$Drug_Info)

# 2. Run ScRank normally, no more slot/layer errors
obj <- my_CreateScRank(
  input = af.FB,
  species = 'mouse',
  cell_type = 'sub.cluster',
  target = 'Irf1'
)


obj <- Constr_net(obj,n.core_cp=1,n.core=10)
obj <- rank_celltype(obj)

obj@para$target

obj <- init_mod(obj)

obj@cell_type_rank

plot_drug_function(
  obj,
  celltype = 'Inflammatory_State',
  category = "H",
  show_leading_edge = T
)
scRank_GSEA(obj, celltype = 'Inflammatory_State', disease_gene = 'Irf1')

plot_net(obj, mode = "heatmap", celltype = "Inflammatory_State")
plot_net(obj, mode = "network", celltype = "Inflammatory_State")
plot_dim(obj)



plot_net(obj, mode = "heatmap", celltype = "Inflammatory_State")





###Gene correlation function############
library(ggstatsplot)
library(ggpubr)


norm_matrix <- GetAssayData(object = afe, assay = "SCT", slot = "scale.data")

norm_matrix<-t(norm_matrix)
norm_matrix<-as.data.frame(as.matrix(norm_matrix))

p5<-ggscatterstats(norm_matrix, 
                   y ='Tsc22d3', 
                   x ='Slc2a1',
                   type = "spearman",
                   centrality.para = "mean",                              
                   margins = "both",                                         
                   xfill = "#009E73",
                   yfill = "#D55E00", 
                   marginal.type = "boxplot",    #the type can be density, boxplot, violin, densigram
                   title = "Relationship between")
print(p5)





####cNMF enrichment
cNMF_MP <- read.delim("./iwd/iwd.gene_spectra_score.k_12.dt_0_1.txt",
                      row.names = 1,
                      stringsAsFactors = FALSE)
cNMF_MP<-as.data.frame(t(cNMF_MP))
# iterate over each column with lapply
cNMF_enrich_list <- lapply(colnames(cNMF_MP), function(col_name) {
  # get the current column
  col_data <- cNMF_MP[, col_name, drop = FALSE]
  
  # sort by value, take the top 200
  top_200 <- head(order(col_data[, 1], decreasing = TRUE), 200)
  
  # extract gene names
  top_genes <- rownames(cNMF_MP)[top_200]
  
  # create a data frame
  data.frame(
    gene = top_genes,
    cluster = col_name,
    stringsAsFactors = FALSE
  )
})

# merge
cNMF_enrich <- do.call(rbind, cNMF_enrich_list)
Symbol <- mapIds(get("org.Mm.eg.db"), keys = cNMF_enrich$gene, keytype = "SYMBOL", column="ENTREZID")
ids <- bitr(cNMF_enrich$gene,"SYMBOL","ENTREZID", "org.Mm.eg.db")
data <- merge(cNMF_enrich, ids, by.x="gene", by.y="SYMBOL")
gcSample <- split(data$ENTREZID, data$cluster)


###################
#nuclear-to-cytoplasmic ratio

genes <- list(c("nc_ratio"))
names(genes) <- 'gene'

#afe <- AddModuleScore_UCell(uterus,features=genes,name="_score")
markers=c('Creb3l1','Sox9','Klf4','Fosl1')
#extract the data
library(ggpubr)
afe<- subset(af.ST, Cell_Anno %in% c('vSMC','FMC_2','FMC_4','FMC_7'))
af.ST$domain
#afe<-af.ST
cell_anno <- afd$Celltype
#gene_expr <-  afe$nc_ratio
gene_expr <- GetAssayData(afd, assay = "SCT", slot = "data")["FOSL1", ]
df1 <- data.frame(
  Cell_ID = names(gene_expr),  # cell name
  Gene_Expression = gene_expr, # Nos2 expression value
  Cell_Annotation = cell_anno  # cell annotation
)

#set the comparison groups
my_comparisons1 <- list(c("FMC_2", "vSMC")) 
my_comparisons2 <- list(c("FMC_2", "FMC_4")) 
my_comparisons3 <- list(c("FMC_2", "FMC_7")) 
my_comparisons4 <- list(c("FMC_4", "vSMC")) 
#violin plot

ggplot(df1, aes(x = Cell_Annotation, y = Gene_Expression, fill = Cell_Annotation)) + 
  # change geom_violin to geom_boxplot
  geom_boxplot(
    color = 'black',      # border color
    size = 0.8,           # border thickness
    width = 0.6,          # boxplot width
    outlier.shape = 16,   # outlier shape
    outlier.size = 1.5,   # outlier size
    outlier.color = "black",  # outlier color
    alpha = 0.8           # fill transparency
  ) +
 # coord_cartesian(ylim = c(0, 8)) +
  # statistical comparison
  stat_compare_means(
    method = "t.test",
    hide.ns = F,
    comparisons = c(my_comparisons1,my_comparisons2,my_comparisons3,my_comparisons4),
    label = "p.signif",
    bracket.size = 0.8,
    size = 6,y.position = 0.5 
  ) +
  
  # limit the y-axis display range
  coord_cartesian(ylim = c(-0.4,0.5)) +  # replace with your own range
  
  scale_fill_manual(values = c('#00a7fa','#354c85','#31aaaa','#9d9ac6') ) +
  
  theme_classic() 



markers=c('Myh11','Tagln','Acta2','Lum','Dcn','Bgn','Fn1','Col1a1','Col1a2',
          'Col3a1' ,'Sox9','Mmp2','Mmp14','Vcam1','Icam1','Angptl4',
          'Slc43a3','Thbs1','C3','C4b','Csf1',
          'Cx3cl1','Cxcl12','Mif','Lgals3','Spp1' ,#FMC_2
          'Tcf4','Runx1'
)
markers <- c(
  'MYH11', 'TAGLN', 'ACTA2', 'LUM', 'DCN', 'FN1', 'COL1A1', 'COL1A2',
  'COL3A1', 'VCAM1', 'ICAM1', 'SPP1','LDHA',"TCF4"  , "MIF","LGALS3"
  
)

markers=c('CREB3L1','SOX9','KLF4','FOSL1')
markers=c('Creb3l1','Sox9','Klf4','Fosl1')
markers=c('Ldha','Pkm','Hif1a','Ldhb')
library(scCustomize)
pdf(file = "FMC DOT .pdf",width =10,height = 5)
Clustered_DotPlot(seurat_object = afd, 
                  colors_use_exp = c( '#64A2CC',"white",'#E07B56',"#731929"),
                  colors_use_idents=c(  "#00a7fa" ,"#354c85","#31aaaa","#9d9ac6"  ),
                  flip=T,x_lab_rotate=90,
                  features = markers)
dev.off()


