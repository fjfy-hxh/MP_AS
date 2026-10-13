# =============================================================================
# 04_kegg_go_enrichment.R
# -----------------------------------------------------------------------------
# GO / KEGG enrichment analysis, GSEA, and enrichment network plots
#
# Features
#   - FindAllMarkers / FindMarkers to get marker genes, bitr conversion to ENTREZID
#   - compareCluster(enrichGO / enrichKEGG) grouped enrichment, ENTREZID back to SYMBOL
#   - Custom enrichment bar plot (-log10(p) horizontal bars + pathway and gene labels)
#   - KEGG enrichment lollipop plot for a single-group contrast (ggforce::geom_link)
#   - KS_GSEA for GSEA analysis and enrichment curve plots
#   - aPEAR enrichmentNetwork for enrichment network plots
#
# Main inputs
#   - Seurat object (with clustering results and marker genes)
#   - Species annotation databases org.Hs.eg.db / org.Mm.eg.db
#
# Main outputs
#   - GO.pdf, KEGG / GO enrichment bar plots, GSEA curve plots, enrichment network plots
#
# Dependencies
#   Repo-internal custom functions: ks_plot_gsea.R
#   Main R packages: clusterProfiler, org.Hs.eg.db, org.Mm.eg.db, aPEAR, ggforce, dplyr
#
# Notes
#   Organized from analysis drafts; the original logic, algorithm and workflow are unchanged.
#   Before running, adjust the input / output paths as needed (the original drafts were
#   mostly local absolute paths), and replace the example sample names, cell type names
#   and comparison groups with your own data. Run from the repository root (open
#   sc-spatial-analysis-scripts.Rproj in RStudio) so source("./R/...") resolves correctly.
# =============================================================================





library(Seurat)
library(SeuratData)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
afc<- subset(af, idents= c("PT subcluster 1","Other TUMOR PT subcluster"))

obj.markers <- FindAllMarkers(afe, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
obj.markers
Symbol <- mapIds(get("org.Mm.eg.db"), keys = obj.markers$gene, keytype = "SYMBOL", column="ENTREZID")
ids <- bitr(obj.markers$gene,"SYMBOL","ENTREZID", "org.Mm.eg.db")

# Merge ENTREZID into obj.markers
data <- merge(obj.markers, ids, by.x="gene", by.y="SYMBOL")


## The split() function groups vectors, matrices, and data frames by a grouping factor.
## It returns a list representing the observations for each level of the grouping variable.
gcSample <- split(data$ENTREZID, data$cluster)
gcSample

xxGO <- compareCluster(gcSample,
                     fun = "enrichGO",
                     OrgDb = "org.Mm.eg.db",
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05
)
#xx <- compareCluster(gcSample, fun="enrichGO", keytype="SYMBOL", OrgDb="org.Hs.eg.db", ont="BP", pvalueCutoff=1, qvalueCutoff=1)

xxKEGG <- compareCluster(gcSample, fun="enrichKEGG", organism="hsa", pvalueCutoff=0.05)

#dotplot(xxKEGG,showCategory=5,font.size = 5)+
#  theme(axis.text.y = element_text(size = 10))+ theme(axis.text.x = element_text(size = 10))  

res<- xxGO@compareClusterResult


## Convert ENTREZID back to SYMBOL in the enrichment results
for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
  res[i,"geneID"] = gene_names
}
summary(gcSample)
res$Description[15:250]
View(res)

enrich <- res %>% 
  group_by(Cluster) %>% 
  filter(Description %in% c( 
                            "cellular response to nutrient levels", 
                            "fatty acid oxidation",
                            'lipid oxidation',
                            'cellular ketone metabolic process',
                            'response to reactive oxygen species',
                            'glutamate receptor signaling pathway',
                            'reactive oxygen species metabolic process'
                            
                            )) %>%
  filter(Cluster %in% c("CWD"))

dt <- enrich
dt <- dt[order(dt$Cluster), ]
dt$Description <- factor(dt$Description, levels = dt$Description)
colnames(dt)
View(dt)

# First define a custom theme:
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # Remove the y-axis pathway labels in the custom theme:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)


p <- ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description), fill = Cluster)) +
  scale_fill_manual(values =c("#6bb9d2","#d55640")) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # Make the bar bottom flush with the y-axis
  labs(x = "-Log10(pvalue)", y = "IWD down regulate",#"CD16+ monocyte         CD14+ monocyte"
          title = "GO Pathway enrichment") +
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0, left-aligned
  geom_text(size=3, aes(x = 0.15, label = geneID), hjust = 0, vjust = 2.5, 
            color=rep(c('white')### if there is only one, provide a single color, '#d55640'
                      ,each=7## this each value is the number of gene sets @@@@@@@@@@@@@@@@@@@@@@
                      )) + # hjust = 0, left-aligned
  theme_classic() + 
  mytheme +
  NoLegend()

p

# Save; the width and height were adjusted here to make the result look better
#ggsave(filename = "KEGG.png", width = 5.2, height = 5.1, plot = p)




#####GO enrichment analysis######################

res<- xxGO@compareClusterResult


## Convert ENTREZID back to SYMBOL in the enrichment results
for (i in 1:dim(res)[1]) {
  arr = unlist(strsplit(as.character(res[i,"geneID"]), split="/"))
  gene_names = paste(unique(names(Symbol[Symbol %in% arr])), collapse="/")
  res[i,"geneID"] = gene_names
}
summary(gcSample)
res$Description
View(res)

enrich <- res %>% 
  group_by(Cluster) %>% 
  filter(Description %in% c( 
    "cellular response to peptide hormone stimulus", 
    "antigen processing and presentation of peptide antigen via MHC class I", 
    "positive regulation of cytokine production involved in immune response", 
    "ERBB signaling pathway"
    )) %>%
  filter(Cluster %in% c("Low PCDH7 Tumor C7", "High PCDH7 Tumor C7"))

dt <- enrich
dt <- dt[order(dt$Cluster), ]
dt$Description <- factor(dt$Description, levels = dt$Description)
colnames(dt)
View(dt)

# First define a custom theme:
mytheme <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # Remove the y-axis pathway labels in the custom theme:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  legend.title = element_text(size = 13),
  legend.text = element_text(size = 11),
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)


p <- ggplot(data = dt, aes(x = -log10(pvalue), y = rev(Description), fill = Cluster)) +
  scale_fill_manual(values =c('#6bb9d2', '#d55640')) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # Make the bar bottom flush with the y-axis
  labs(x = "-Log10(pvalue)", y = "High PCDH7 Tumor C7",#"CD16+ monocyte         CD14+ monocyte"
       title = "GO Pathway enrichment") +
  geom_text(size=3.8, aes(x = 0.05, label = Description), hjust = 0) + # hjust = 0, left-aligned
  geom_text(size=3, aes(x = 0.15, label = geneID), hjust = 0, vjust = 2.5, 
            color=rep(c('#6bb9d2')### if there is only one, provide a single color, '#d55640'
                      ,each=4## this each value is the number of gene sets @@@@@@@@@@@@@@@@@@@@@@
            )) + # hjust = 0, left-aligned
  theme_classic() + 
  mytheme +
  NoLegend()

p

# Save; the width and height were adjusted here to make the result look better
#ggsave(filename = "GO.png", width = 5.2, height = 5.1, plot = p)


table(afd$caftype)
receiver<-'my_iCAF'

#markers <- FindAllMarkers(afd, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
afd$cell_cluster = ifelse(afd$caftype==receiver, receiver,"other_cells")##### select the overall receiver cells
Idents(afd)<-afd$cell_cluster
markers<-FindMarkers(afd, ident.1=receiver, ident.2="other_cells",
                     only.pos=TRUE)
Idents(afd)<-afd$caftype##### restore defaults

library(ggforce)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)

df <- data.frame(gene =markers$gene,
                 avg_log2FC = markers$avg_log2FC,
                 p_val_adj = markers$p_val_adj)
###KEGG analysis
entrezid_all = mapIds(x = org.Hs.eg.db,  #modify according to your own species
                      keys = df$gene, 
                      keytype = "SYMBOL", 
                      column = "ENTREZID")
entrezid_all  = na.omit(entrezid_all) 
entrezid_all = data.frame(entrezid_all) 


KEGG_enrich = enrichKEGG(gene = entrezid_all[,1], 
                         keyType = "kegg",
                         pAdjustMethod = 'BH',  
                         organism= "hsa", #"hsa","mmu"
                         qvalueCutoff = 0.05,
                         pvalueCutoff=0.05)
KEGG_enrich  = KEGG_enrich@result
### select some for display
KEGG_enrich$Description

KEGG_enrich <- KEGG_enrich %>%
  filter(Description %in% c( 
    "Th17 cell differentiation", 
    "Th1 and Th2 cell differentiation", 
    "NF-kappa B signaling pathway", 
    "PI3K-Akt signaling pathway",
    'NOD-like receptor signaling pathway',
    "Metabolism of xenobiotics by cytochrome P450",
    "Notch signaling pathway",
    "Toll-like receptor signaling pathway",
    "Wnt signaling pathway",
    "T cell receptor signaling pathway",
    "VEGF signaling pathway" ,
    "Chemokine signaling pathway",
    "RIG-I-like receptor signaling pathway"
  ))


#KEGG_enrich <- KEGG_enrich[c(1:15),]#### multi-group display
KEGG_enrich$`-log10(p)` <- -log10(KEGG_enrich$pvalue)


KEGG_enrich$Description <-  sub("- .*", "", KEGG_enrich$Description)


### plotting
ggplot(KEGG_enrich) +
  ggforce::geom_link(data = KEGG_enrich, 
                     aes(x = 0,
                         y = reorder(Description, -`-log10(p)`),
                         xend =`-log10(p)`,
                         yend = Description,
                         colour="grey",
                         alpha = after_stat(index),
                         size = after_stat(index)),
                     show.legend = F)+
  geom_point(data = KEGG_enrich,
             aes(x = `-log10(p)`, y = Description),
             size=6,
             shape = 21,
             fill = "white")+
  scale_colour_manual(values =c("#96C3D8"))+
  geom_text(mapping=aes(x=`-log10(p)`,
                        y=Description,
                        label=Count),
            size =2.5,color="black")+
  geom_vline(xintercept = -log10(0.05),
             color = "grey60",
             linetype=2, cex=0.5)+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(colour = 'black'))+
  labs(x="-Log10(P-value)")

#########################################################################3
### multi-group enrichment
#########################################
markers<-FindAllMarkers(afd,  only.pos=TRUE)
# marker genes and enrichment analysis; of course, this also works for DEGs of each celltype
sce_marker<-markers


# gene enrichment analysis for each celltype group
group <- data.frame(gene=sce_marker$gene,group=sce_marker$cluster)
Gene_ID <- bitr(sce_marker$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')

data_GO <- compareCluster(ENTREZID~group, 
                          data=data, 
                          fun="enrichGO", #enrichGO  enrichKEGG
                          OrgDb="org.Hs.eg.db",
                          ont = "BP",  # comment out this line for KEGG
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.01)

data_GO_sim <- simplify(data_GO, cutoff=0.7, by="p.adjust", select_fun=min)
data_GO_sim <- data_GO_sim@compareClusterResult
View(data_GO_sim)

# select top 10 terms per celltype for visualization
data_GO_sec <- data_GO_sim %>% group_by(group) %>% top_n(5, -p.adjust)
data_GO_sec$`-log10(p)` <- -log10(data_GO_sec$p.adjust)
# display the pathways to be shown together
p2 = ggplot(data_GO_sec) +
  ggforce::geom_link(data = data_GO_sec, 
                     aes(x = 0,
                         y = reorder(Description, -`-log10(p)`),
                         xend =`-log10(p)`,
                         yend = Description,
                         colour=group,
                         alpha = after_stat(index),
                         size = after_stat(index)),
                     show.legend = F)+
  geom_point(data = data_GO_sec,
             aes(x = `-log10(p)`, y = Description),
             size=6,
             shape = 21,
             fill = "white")+
  scale_colour_manual(values = colors9)+
  geom_text(mapping=aes(x=`-log10(p)`,
                        y=Description,
                        label=Count),
            size =2.5,color="black")+
  geom_vline(xintercept = -log10(0.05),
             color = "grey60",
             linetype=2, cex=0.5)+
  theme_classic()+
  theme(axis.title.y = element_blank(),
        axis.text = element_text(colour = 'black'))+
  labs(x="-Log10(P-value)")

pdf(file = "GO.pdf",width = 8,height = 6)
p2
dev.off()


################################################################
#######GSEA enrichment
################################################################
source(file = "./R/ks_plot_gsea.R")
focus_cell
markers
table(Idents(afd))
selected_idents <- c("EC","CAF")
afd<- subset(af.stroma, idents = selected_idents)
afd<-JoinLayers(afd)
df <- FindMarkers(afd, min.pct = 0, logfc.threshold = 0,
                  group.by = "xx",ident.1 ='EC',ident.2='Mono')

df$gene = rownames(df)
View(df)

#GSEA analysis-----------------------------------------------------------------------
res <- KS_GSEA(gene = df$gene,
                   LogFC = df$avg_log2FC,
                   analysis = "GO",#'GO',"KEGG"
                   package = 'clusterProfiler',
                   OrgDb = 'org.Hs.eg.db')#"org.Hs.eg.db", "org.Mm.eg.db"

View(res@result)

res@result$Description

selected_terms <- res@result %>%
  filter(Description %in% c("NF-kappa B signaling pathway - Mus musculus (house mouse)", 
                            "TNF signaling pathway - Mus musculus (house mouse)", 
                            'IL-17 signaling pathway - Mus musculus (house mouse)',
                            'Cytokine-cytokine receptor interaction - Mus musculus (house mouse)',
                            'Toll-like receptor signaling pathway - Mus musculus (house mouse)'
                            
                            ))
selected_terms <- res@result %>%
  filter(Description %in% c("cell-cell adhesion via plasma-membrane adhesion molecules",
                            "blood vessel remodeling"
  ))
ks_plot_gsea(data = res,
             term =res@result$Description,
             rankColor="#90191B",
             labels_size=T)


ks_plot_gsea(data = res,
             term = selected_terms$Description,
             rankColor="#90191B",
             labels_size=F)



# enrichment network plot
library(aPEAR)
library(grDevices)

# the enrichment results must be numerous enough
enrichmentNetwork(res@result, drawEllipses = TRUE, colorBy = 'p.adjust',
                  colorType = 'pval',fontSize = 3)+
scale_color_gradientn(colours =colorRampPalette(c('#64A2CC',"white",'#E07B56',"#731929"))(100),
                      breaks = c(0, log(0.05), log(0.01), log(0.001),log(0.0001)),   # specify tick positions
                      limits = c(log(0.0001), log(1)),
                      labels = c("p>0.05", "p<=0.05", "p<0.01", "p<0.001", "p<0.0001"),
                      name = "Log(p.adjust)")# custom labels
