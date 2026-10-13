# =============================================================================
# 02_sc_subclustering.R
# -----------------------------------------------------------------------------
# Cell subpopulation extraction and reclustering, subpopulation annotation, UMAP gating deletion and over-clustering detection
#
# Features
#   - Extract target cell populations for reclustering and manually annotate subpopulations (mapvalues mapping)
#   - Write subpopulation annotations back to the major object, updating Celltype
#   - Filter / delete cells by UMAP coordinate conditions (UMAP gating function)
#   - Subpopulation gene expression violin plot (with significance) and DotPlot
#   - scSHC over-clustering detection
#
# Main inputs
#   - prcc_scRNA.qs and other Seurat objects (prepare them yourself)
#
# Main outputs
#   - subcluster_anno_UMAP.pdf, gene_sig_cell_expression_diff.pdf and other plots
#
# Dependencies
#   In-repo custom functions: vnplot.R, ks_VlnExp.R, subset_seurat_by_umap.R
#   Main R packages: Seurat, scSHC, ggpubr, patchwork, plyr
#
# Notes
#   This script was compiled from analysis drafts, retaining the original analysis logic without changing algorithms or workflow.
#   Before running, modify the input / output paths appearing in the script as needed (the original drafts mostly used local absolute paths),
#   and replace the sample names, cell type names, comparison groups, etc. in the examples with your own data.
#   It is recommended to first set the repository root as the working directory (open
#   sc-spatial-analysis-scripts.Rproj in RStudio) so that source("./R/...") resolves correctly.
# =============================================================================

af <- qread("prcc_scRNA.qs")
qsave(af, file = "prcc_scRNA.qs")


###############Violin plot of gene expression in cells with significance
#group diff
source(file = "./R/vnplot.R")
gene_sig <- 'Pomc'
#cell diff

comparisons <- list()
comp=combn(names(table(af$Type)),2)
af$Type
#names(table(af$cellType))
for(j in 1:ncol(comp)){comparisons[[j]]<-comp[,j]}
afvp(afe,gene_signature = gene_sig, file_name = "Slc2a1", test_sign = comparisons,pta=0,
     cols=colors6,label="p.signif",group="Type",widplot=4,heiplot=4,ak=1)

af$fat_acid
######Second method
source("./R/ks_VlnExp.R")


# ks_VlnExp(object = uterus, group="orig.ident",
#           group_order=c("HC","AEH","EEC"),
#           features="ANXA1",comparisons=c(my_comparisons1,my_comparisons2,my_comparisons3))
# table(af$Type)
afe

table(Idents(af))
my_comparisons3 <- list(c("IWD", "CWD"))

ks_VlnExp(object = afd, group="Type",group_order=c("IWD",'CWD'),
          features='Acly',cols=c(  "indianred","royalblue"),
          comparisons=c(list(c("IWD", "CWD"))
          ))

ks_VlnExp(object = af.niche2, group="Type",group_order=c("CON","MP"),
          features=gene_sig,cols=c("royalblue" ,"indianred" ),
          comparisons=c(list(c("CON", "MP"))
                        ))
afe<-afd
ggsave("gene_sig_cell_expression_diff.pdf",width = 4.5,height = 6.5)


table(afd$Hip)

afd[["RNA"]] <- as(afd[["RNA"]], "Assay")
table(af$Celltype)
afd<- subset(af.sc, Celltype %in% c(  'FBs C1' ))


#DefaultAssay(afd) <- "RNA"

afd <- NormalizeData(afd, normalization.method = "LogNormalize", scale.factor = 10000)
afd <- FindVariableFeatures(afd, selection.method = "vst", nfeatures = 2000)
scale.genes<-VariableFeatures(afd)
afd <-ScaleData(afd,features = scale.genes)
afd <- Seurat::RunPCA(afd, features = VariableFeatures(object = afd))
pct <- afd[["pca"]]@stdev / sum(afd[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
ElbowPlot(afd)$data   %>% ggplot() +
  geom_point(aes(x = dims,y = y_data)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")

afd <- Seurat::RunPCA(afd,npcs = 14)

#afd <- RunHarmony(afd, group.by.vars = "Type")####no need for PCA anymore


afd=FindNeighbors(afd, dims = 1:14, reduction = "pca")
afd <- RunUMAP(afd, reduction = "pca", dims = 1:14)

afd <- FindClusters(afd, resolution = 0.2)

table(afd@active.ident)
afd$sub.cluster<-Idents(afd)
table(afd$sub.cluster,afd$Type)


afgenes=c('Tnf','Il1b','Spp1','S100a9','S100a8','Mfap4','Timp1','Nr4a1',
          'Trem2','Fabp5','Stab1','Mrc1','Nampt','Cd1c','Clec10a')

afgenes=c('Cxcl12','Des','Lrg1','Ccl2','Pdgfra','Acta2','Timp3','Mfap4','Pi16')

DotPlot(afd, features = afgenes)+theme_bw()+#remove background, rotate plot
  
  theme(panel.grid = element_blank(),
        
        axis.text.x=element_text(angle=45,hjust = 1,vjust=1))+#display text at 90 degrees
  
  scale_color_gradientn(values = seq(0,1,0.2),
                        colours = c('#64A2CC',"white",'#E07B56',"#731929"))+#color gradient settings
  # c('#64A2CC',"white",'#E07B56',"#731929")
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))+coord_flip()

afd.markers <- FindAllMarkers(afd, only.pos = T, min.pct = 0.25,  logfc.threshold = 0.1)

ann.ids <- c("Inf Mac",  #0Injured Proximal tubular cell  not a cluster, it is my numbering, cluster+1
             "LAMs",  #1Injured Proximal tubular cell
             "Res Mac",  #2
             "DCs", #3
             "Cycling Myeloid", #4
             "Mono"#
)
ann.ids <- c("Transitional_State",  #0Injured Proximal tubular cell  not a cluster, it is my numbering, cluster+1
             "Progenitor_State",  #1Injured Proximal tubular cell
             "Inflammatory_State",
             'remove'#2
)
ann.ids <- c("vSMCs",  #0Injured Proximal tubular cell  not a cluster, it is my numbering, cluster+1
             "vSMCs",  #1Injured Proximal tubular cell
             "FBs C1",  #2
             "FBs C2", #3
             "FBs C2", #4
             "FBs C3",#5
             "vSMCs", #6
             "FBs C4",#7
             "FBs C5",#8
             "vSMCs",#9
             "vSMCs",#10
             "FBs C6",#11
             "FBs C7",#12
             "FBs C8",#13
             "remove"#14
             
)
Idents(afd)

afidens=mapvalues(Idents(afd), from = levels(Idents(afd)), to = ann.ids)
Idents(afd)=afidens
afd$sub.cluster=Idents(afd)

pdf(file = "subcluster_anno_UMAP.pdf",width =7,height = 5.5)
DimPlot(afd, reduction = "umap", label = F, label.size = 2.5,
       group.by = "sub.cluster",
        pt.size = 1)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
dev.off()
table(Idents(af.all))

afe$sub.cluster<-Idents(afe)
#################Assign the idents from the subpopulation back to the major object
Idents(af, cells = colnames(afd)) <- Idents(afd)
af$Celltype<-Idents(af)


########################################################
selected_idents <- c("PT subcluster 1","PT subcluster 2","PT subcluster 3",
                     "PT subcluster 4","PT subcluster 5")
af<- subset(af, idents = selected_idents)
#######Delete cells based on UMAP#########################
source(file = "./R/subset_seurat_by_umap.R")

af_subset <- subset_seurat_by_umap(
  seurat_obj = af,
  celltype_col = "celltype",  # custom column name
  target_celltype = "vSMC",
  umap1_condition = c(">", 5),##!=
  umap2_condition = c("!=", -4),##!=
  operator = "&"
)
removed_cells <- setdiff(colnames(af), colnames(af_subset))
message("Removed ", length(removed_cells), " cells")


DimPlot(af_subset, reduction = "umap", label = T, group.by = "celltype",
        raster=FALSE,label.size = 2.5,
        pt.size =1,cols= colors2)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
############
##Update the af object
af<-af_subset
af$celltype
table(af$celltype)
 ####
# Give the major object metadata a new name celltype_major_Fib_sub to store the major object info containing detailed FIB annotations
af$afd_subcluster <- Idents(af)
Idents(af1)<-af1$afd_subcluster
table(Idents(af1))

af1<-af

table(Idents(af))

#######Select the subpopulation object###########
selected_idents <- c("PT subcluster 1","PT subcluster 2","PT subcluster 3",
                     "PT subcluster 4","PT subcluster 5")

selected_idents <- c("PT subcluster 1","PT subcluster 2","PT subcluster 3",
                     "PT subcluster 4","PT subcluster 5","T cell","Dendritic cell",
                     "Natural killer cell","Macrophage","Neutrophil",
                     "Monocyte", "B cell")

af<- subset(af, idents = selected_idents)
afd <- subset(afd, Celltype == "vSMC")


table(af$Celltype)
#############Delete the object of a certain subpopulation
af<- subset(af, subset = Celltype != "remove")
table(afd$ec)
#Idents(af.1)<-"Proximal tubule cell"
# Create a new Seurat object
af.ST
af<- subset(af1, idents = selected_idents)

af1$afd_subcluster<-Idents(af1)
Idents(af.bei)<-'FBs like'
##########Assign subpopulation annotations back to the major object###
Idents(af, cells = colnames(af.bei)) <- Idents(af.bei)


####Set level###################################
new_levels <- c(
  "EC", 
  "T cell", 
  "LAMs", "Res Mac", "Inf Mac", 
  "DCs", "Mono", 'Cycling Myeloid','SC','vSMCs',
  "FBs C1", "FBs C2", "FBs C3", "FBs C4", "FBs C5", "FBs C6", "FBs C7",'FBs C8'
)

# Set factor levels
af$Celltype <- factor(af$Celltype, levels = new_levels)


table(af$Celltype)


##########Identify over-clustering#########################################
library(scSHC)
scSHC <- testClusters(afd@assays$RNA$counts,cluster_ids=afd$sub.cluster
                      ,cores = 1)

scSHC <- scSHC(afd@assays$RNA$counts,num_PCs=13,
                      ,cores = 1)
table(scSHC[[1]],afd$sub.cluster)

scSHC[[1]]

afd$subcluster <- scSHC[[1]][colnames(afd)]
