# =============================================================================
# 06_xenium_basic_analysis.R
# -----------------------------------------------------------------------------
# Xenium / VisiumHD spatial basic analysis: QC, dimensionality reduction, RCTD annotation, Banksy niche, protein re-annotation
#
# Features
#   - Read VisiumHD polygons (bin.size = 'polygons') and cell segmentation data
#   - Build a spatial reduction layer from polygon centroids, write sample group info
#   - QC filtering, normalization, highly variable genes, PCA / UMAP
#   - spacexr / RCTD spatial cell type annotation (RunRCTD_Annotations), output annotation CSV
#   - Banksy spatial niche clustering (computeBanksy / clusterBanksy), map domains
#   - EC subcluster re-annotation based on the Protein assay (custom reduction + clustering)
#
# Main inputs
#   - VisiumHD data directory, GROUP.csv, visiumHD niche.csv
#   - Xenium+CODEX merged object .qs, single-cell reference object (for RCTD)
#
# Main outputs
#   - af_ST_celltype.csv, lesion_banksy.csv, af_ST_domains.csv, etc.
#
# Dependencies
#   Repo-internal custom functions: cell_segmentation_spatial_io.R
#   Main R packages: Seurat, SeuratWrappers, Banksy, SpatialExperiment, SummarizedExperiment, spacexr, sceasy, reticulate, hdf5r, qs, BiocParallel
#
# Notes
#   Organized from analysis drafts; the original logic, algorithm and workflow are unchanged.
#   Before running, adjust the input / output paths as needed (the original drafts were
#   mostly local absolute paths), and replace the example sample names, cell type names
#   and comparison groups with your own data. Run from the repository root (open
#   sc-spatial-analysis-scripts.Rproj in RStudio) so source("./R/...") resolves correctly.
# =============================================================================

df.group<-read.csv( "./visiumHD niche.csv")

af.ST$niche <- df.group$Niche[match(colnames(af.ST),
                                    df.group$Barcode)]


table(af.ST$CN,af.ST$Type)
library(Seurat)
library(argparse)
library(dplyr)
library(ggplot2)
library(Banksy)
library(SeuratWrappers)
library(harmony)
library(qs)
library(reticulate)
library(hdf5r)
library(sceasy)
library(BiocParallel)
library(dplyr)
library(spacexr)
library(qs)
source("./R/cell_segmentation_spatial_io.R")

#qsave(af.ST, file = "MP_spatial_overall.qs")



af.ST <- Load10X_Spatial(
  data.dir = "./visiumHD" ,bin.size = c("polygons")
)

df.group<-read.csv( "./visiumHD/GROUP.csv")

af.ST$group <- df.group$GROUP[match(colnames(af.ST),
                                    df.group$Barcode)]


coords_fov1 <- as.data.frame(af.ST@images$slice1.polygons$centroids@coords)
cells_fov1 <- af.ST@images$slice1.polygons$centroids@cells
rownames(coords_fov1) <- cells_fov1
coords_fov1 <- as.matrix(coords_fov1)
colnames(coords_fov1) <- paste0("spatial_", 1:ncol(coords_fov1))

# 9. Create the spatial reduction
spatial_reduction <- Seurat::CreateDimReducObject(
  embeddings = coords_fov1,
  key = "spatial_",
  assay = DefaultAssay(af.ST)
)

# 10. Add to the Seurat object
af.ST[["spatial"]] <- spatial_reduction



af.ST <- subset(af.ST, subset = !is.na(group) & group != "NA")




# remove cells with 0 counts
af.ST <- subset(af.ST, subset = nCount_Spatial.Polygons >20)

af.ST <- NormalizeData(af.ST, normalization.method = "LogNormalize", scale.factor = 10000)
af.ST <- FindVariableFeatures(af.ST, selection.method = "vst", nfeatures = 2000)
scale.genes<-VariableFeatures(af.ST)
af.ST <-ScaleData(af.ST,features = scale.genes)
af.ST <- Seurat::RunPCA(af.ST)
pct <- af.ST[["pca"]]@stdev / sum(af.ST[["pca"]]@stdev) * 100 ; cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
ElbowPlot(af.ST)$data %>% ggplot() +
  geom_point(aes(x = dims,y = stdev)) +
  geom_vline(xintercept = pc.use, color = "darkred") +
  theme_bw() + labs(title = "Elbow plot: quantitative approach")

af.ST <- Seurat::RunPCA(af.ST,npcs = 10)
af.ST=FindNeighbors(af.ST, dims = 1:10, reduction = "pca")
af.ST <- RunUMAP(af.ST, reduction = "pca", dims = 1:10)


#############################################################
#############################################################
###########RCTD####################################################
#############################################################
table(af.ST$group)
# Extract the counts matrix (must be a sparse or regular matrix)
af
table(af$Celltype)
Idents(af)<-af$Celltype
af_RCTD<- af
#
cell_counts <- table(af_RCTD$Celltype)
table(af_RCTD$Celltype)
#celltypes_to_remove <- names(cell_counts[cell_counts < 50])
#af_RCTD <- subset(af_RCTD, subset = Celltype != 'NK cell')

counts <- GetAssayData(af_RCTD, assay = "RNA", layer = "counts")

# Extract cell types (assumed to be stored in the "cell_type" column of meta.data)
cell_types <- Idents(af_RCTD)  # or Idents(af)

# Create the Reference object
reference_RCTD <- Reference(counts, cell_types)

############# run RCTD with a custom function
af.ST <- RunRCTD_Annotations(
  seurat_obj = af.ST,
  assay_name = "Spatial",
  data_type = "pcf",#visiumHD #pcf  if reduction=spatial has already been done, pcf can be used directly
  reference_RCTD = reference_RCTD,
  UMI_min = 20,  # this threshold can be adjusted
  max_cores = 25
)

af.ST$Cell_Anno[is.na(af.ST$Cell_Anno)] <- "low_Qs"

table(af.ST$Cell_Anno)
############# remove low_Qs
af.ST <- subset(af.ST, subset = Cell_Anno != "low_Qs")



########output celltype csv###########################
af.ST_both<-sce_subset
barcodes <- colnames(af.ST)###seurat object
cell_anno_values <- af.ST$Cell_Anno

# 4. Build the data frame
output_df <- data.frame(
  cell_id = barcodes,
  group = cell_anno_values,
  row.names = NULL # remove row names to make the data frame tidier
)
head(output_df)

write.csv(output_df, "af_ST_celltype.csv", row.names = FALSE)


table(af.ST$group)

###################################################
stage_vector <- rep(NA, length(af.ST$group))
stage_vector[af.ST$group %in% c("MP1", "MP11",'MP2','MP3','MP6','MP9')] <- 'MP'
stage_vector[af.ST$group %in% c("CON2", "CON3",'CON5','CON6','CON9')] <- 'CON'
# Add the stage vector to the Seurat object metadata
af.ST$Type <- stage_vector


######################################################## 
#######Banksy##################################################
###########################################
library(Banksy)
library(SummarizedExperiment)
library(SpatialExperiment)
library(scuttle)
library(scater)
library(cowplot)
library(ggplot2)


banksy_seurat<- af.ST
Idents(banksy_seurat)<-banksy_seurat$Cell_Anno
#banksy_seurat<-af_combined####### use seurat for the banksy part

banksy_seurat <- NormalizeData(banksy_seurat, normalization.method = "LogNormalize",
                               , verbose = FALSE,scale.factor = 10000)

banksy_seurat <- FindVariableFeatures(banksy_seurat, 
                                      assay = "Spatial",          # specify the assay to use
                                      selection.method = "vst", # choose the variance-stabilizing method
                                      nfeatures = 2000)  


ST_count <- GetAssayData(
  object = banksy_seurat, 
  assay = "Spatial",  # use the RNA assay by default
  layer =  "count" # get the raw counts
)

#banksy_coords <- as.data.frame(banksy_seurat@images$fov$centroids@coords)
#banksy_coords<-as.matrix(banksy_coords)

#seurat_cells <- colnames(banksy_seurat)
#rownames(banksy_coords) <- seurat_cells

banksy_coords <-Embeddings(banksy_seurat, reduction = "spatial")



se <- SpatialExperiment(assay = list(counts = ST_count), spatialCoords = banksy_coords)


#gc()
# QC based on total counts
#qcstats <- perCellQCMetrics(se)
#thres <- quantile(qcstats$total, c(0.05, 0.98))
#keep <- (qcstats$total > thres[1]) & (qcstats$total < thres[2])
#se <- se[, keep]

# Normalization to mean library size
se <- computeLibraryFactors(se)
aname <- "normcounts"
assay(se, aname) <- normalizeCounts(se, log = FALSE)
#se
lambda <- c(0.8)
k_geom <- c(25)

se <- Banksy::computeBanksy(se, assay_name = aname, compute_agf = TRUE, k_geom = k_geom)
set.seed(1000)
gc()


se <- Banksy::runBanksyPCA(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::runBanksyUMAP(se, use_agf = TRUE, lambda = lambda)
se <- Banksy::clusterBanksy(se, use_agf = TRUE, lambda = lambda, resolution = 0.3)
#se <- Banksy::connectClusters(se)


#> clust_M1_lam0.2_k50_res1.2 --> clust_M1_lam0_k50_res1.2
cnames <- colnames(colData(se))
cnames <- cnames[grep("^clust", cnames)]
colData(se) <- cbind(colData(se), spatialCoords(se))

plotColData(se,
            x = "spatial_1", y = "spatial_2",
            point_size = 0.6, colour_by = 'clust_M1_lam0.8_k50_res0.3'
)



### extract the banksy results, assign to the seurat object and save as CSV
df.domain <- data.frame(
  cell_id = colnames(se),
  group = se$clust_M1_lam0.8_k50_res0.3
)

df.domain$group <- paste0("niche", df.domain$group)

banksy_seurat$hip <- df.domain$group[match(colnames(banksy_seurat),
                                       df.domain$cell_id)]



write.csv(df.domain, "lesion_banksy.csv", row.names = FALSE)



df.domain <- df.domain %>%
  mutate(
    domain = case_when(
      group %in% c(1,6,8) ~ "Relatively normal wall",  # when the banksy_niche value is 1 or 3, the domain is tumor
      group %in% c(2,3,10,9,8,7) ~ "Peri vascular", 
      group %in% c(4) ~ "Inflamed lesion",
      group %in% c(5) ~ "Relatively diseased wall",# when the value is 2, the domain is stroma
      TRUE ~ NA_character_                  # for other values (e.g. NA, 4), assign NA (can be changed to "unknown")
    )
  )

df.domain$group<-df.domain$domain
## assign csv to the seurat object
# set the folder path containing the CSV file
df.domain<-read_csv("./lesion_banksy.csv")


af.ST$group
af.ST$niche <- df.domain$group[match(colnames(af.ST),
                                   df.domain$Cell_ID)]
af.DG<- subset(af.DG, DG %in% c("DG mo"))
Idents(af.DG)<-af.DG$DG

af.niche1<- subset(af.hip, Hip %in% c("niche_1"))
af.niche5<- subset(af.hip, Hip %in% c("niche_5"))
Idents(af.niche1)<-af.niche1$DG
Idents(af.niche5)<-af.niche5$DG

Idents(af.hip, cells = colnames(af.DG)) <- Idents(af.DG)


table(banksy_seurat$stroma)

#############################

df.domain<-read.csv( "xenium_cellcharter_CN.csv")

af.ST$niche <- df.domain$group[match(colnames(af.ST),
                                       df.domain$cell_id)]
table(af.ST$niche)
# output csv file

# create a data frame containing cell names and domain info
export_data <- data.frame(
  cell_id = colnames(af.ST),    # first column: cell names
  domain = af.ST$domain  # second column: domain info
)

# write the data frame to a CSV file
write.csv(export_data, 
          file = "af_ST_domains.csv", # specify the file name
          row.names = FALSE)         # do not write row names (usually numeric indices) to the file






###  xenium LEC cell subcluster protein re-annotation
#########################################################################
table(af.ST_both$Cell_Anno)
af.ST_both<-sce_subset
afd<- subset(af.ST_both, Cell_Anno %in% c("EC"))
afd
DefaultAssay(afd) <- "Protein"

afd <- NormalizeData(afd, normalization.method = "LogNormalize", scale.factor = 10000)
VariableFeatures(afd)<-rownames(afd)
scale.genes<-VariableFeatures(afd)
afd <-ScaleData(afd,features = scale.genes)


##custom reduction construction######################################
features <- c('CD31-protein','Podoplanin-protein')
expr_sparse <- GetAssayData(afd, assay = "Protein", slot = "data")[features, ]
# 2. Convert the sparse matrix to a regular matrix (key step) and transpose to cell x feature format
expr_matrix <- as.matrix(expr_sparse)  # convert to a regular matrix
expr_matrix <- t(expr_matrix)  # rows=cells, cols=features
# 3. Fix the embedding matrix column names (must match the reduction key)
colnames(expr_matrix) <- c("protein_1", "protein_2")  # format: key + number
# 4. Standardize the reduction key (must end with "_")
custom_reduction <- CreateDimReducObject(
  embeddings = expr_matrix,  # regular matrix format
  key = "protein_",  # correct format: alphanumeric + underscore
  assay = "Protein"
)
# 5. Add to the Seurat object
afd@reductions$custom_protein <- custom_reduction


#ElbowPlot(afd,reduction = "harmony")
afd=FindNeighbors(afd, dims = 1:2, reduction = "custom_protein")
afd <- RunUMAP(afd, reduction = "custom_protein", dims = 1:2)

afd <- FindClusters(afd, resolution = 0.3)

table(afd@active.ident)
afd$sub.cluster<-Idents(afd)

#rownames(afd)
afd.annogenes <- c('PECAM1','VWF','PROX1','LYVE1','PDPN','Podoplanin-protein'
)
DefaultAssay(afd) <- "Protein"

DotPlot(afd, features = afd.annogenes)+coord_flip()+theme_bw()+#remove background, rotate the plot
  
  theme(panel.grid = element_blank(),
        
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+#render text at 90 degrees
  
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+#color gradient settings
  
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))

ann.ids <- c("EC",  #0Injured Proximal tubular cell  not a cluster, it is my index, cluster+1
             "EC",  #1Injured Proximal tubular cell
             "EC",  #2
             "EC",#3
             "EC",#4
             "EC",#5
             "EC",#6
             "EC",#7
             "EC",#8
             "EC",#9
             "LEC",#10
             'EC',#11
             'EC'#12
)
Idents(afd)

afidens=mapvalues(Idents(afd), from = levels(Idents(afd)), to = ann.ids)
Idents(afd)=afidens
afd$mye.anno=Idents(afd)
DimPlot(afd, reduction = "umap", label = F, label.size = 2.5,
        group.by = "mye.anno",
        pt.size = 1,cols= colors5)+theme_classic()+theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),legend.position = "right")
table(Idents(afd))

DefaultAssay(afd) <- "Protein"

Idents(af.ST_both, cells = colnames(afd)) <- Idents(afd)

