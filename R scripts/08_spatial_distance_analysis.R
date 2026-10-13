# =============================================================================
# 08_spatial_distance_analysis.R
# -----------------------------------------------------------------------------
# Spatial distance analysis (Python precomputed distances + R plotting)
#
# Functions
#   - Read Python-precomputed neighborhood / vertical / centroid distance CSVs and write them into Seurat metadata
#   - Plot cell type / domain percentage stacked bar charts and smoothed area charts by distance bins (percentile binning)
#   - Plot cell proportion curves and gene expression curves as a function of distance
#   - 2D kernel density plot (distance to intima x distance to plaque core)
#   - ggscatterhist 2D scatter + marginal boxplot, with ANOVA by group
#   - Multi-group distance distribution violin plots with significance annotation; inter-group distance difference comparison
#
# Main inputs
#   - FBs C1 radius_neighborhood.csv, vertical_distance*.csv, centroid_distance.csv, etc.
#   - Seurat spatial object (af.ST / af.plaque etc.)
#
# Main outputs
#   - Distance stacked bar charts, distance curve plots, cell 2D density plots, distance difference violin plots
#
# Dependencies
#   In-repo custom functions: spatial_distance_analysis.R
#   Main R packages: Seurat, ggplot2, dplyr, viridis, ggpointdensity, ggpubr, rstatix
#
# Notes
#   This script was organized from analysis drafts, keeping the original analysis logic intact with no changes to algorithms or workflow.
#   Before running, adjust the input / output paths as needed (the original drafts mostly used local absolute paths),
#   and replace the example sample names, cell type names, and comparison groups with your own data.
#   It is recommended to first set the repository root as the working directory (open
#   sc-spatial-analysis-scripts.Rproj in RStudio) so that source("./R/...") resolves correctly.
# =============================================================================


source("./R/spatial_distance_analysis.R")
library(ggplot2)
library(dplyr)
library(viridis)

#########################################################################
##Plot in R after Python computation
#####################################################
##Differences of cells at different distances centered on a given cell
##################################################
af.ST <- add_metadata_from_csv(
  csv_file = "./FBs C1 radius_neighborhood.csv",  # CSV file containing cell information
  seurat_obj = af.ST,         # input Seurat object
  metadata_cols = c('min_distance',"radius_range" )  # column names to add; column names must not contain spaces
)
af.ST_Fib1$Celltype_Epi_Fib
seu_radiusrange<- subset(af.plaque_radiusrange, Type %in% c('CWD'))
seu_radiusrange<- subset(seu_radiusrange, subset = domain != "Peri vascular")
#bin distance by percentile
# extract metadata
seu_radiusrange <- add_distance_percentile(seu_radiusrange,
                                           dist_col = "min_distance", 
                                           group_col = "group", 
                                           probs = seq(0, 1, #you can also customize: probs = c(0, 0.1, 0.5, 0.9, 1))
                                                       0.2),##bin every 20%
                                           new_col_name = "distance_percentile")

#########plot the percentage chart changing with distance
colors_Celltype2<-c('#f3758a','#ff87bbff',"#fbaeb5",#EC
                   '#0F7240', #T
                   '#B7B663',#B
                   '#d5a339',"#f8ff00",'#fccc3a',#Mac
                   "#9d9ac6", #SMC
                   "#bbe379",  "#00a7fa",'#b9e6fb','#354c85',
                   "#b1ddcc", "#C9E1AF", "#31aaaa"#FMC
)
c("#00f3d8" , "#8c86ff" ,"#bf9800")

plot_domain_percentage_bar(
  seurat_obj = seu_radiusrange,
  radius_range_col = "distance_percentile",
  domain_col = "Cell_Anno",show_value = F,xlab = "Distance to intima",
  angle = 45, color_palette =colors_Celltype2 ,
  x_labels = c("0-20%", "20-40%", "40-60%", "60-80%",'80-100%')
)
ggsave("CWD_celltype_distance_stacked_bar.pdf",width = 5,height = 3)


####stacked area chart with smoothed curves
custom_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")

p2 <- plot_domain_percentage(
  seurat_obj = seu_radiusrange,
  radius_range_col = "distance_percentile",##x-axis; all must be categorical variables
  domain_col = "domain",##y-axis proportion
  x_labels = c("0-20%", "20-40%", "40-60%", "60-80%",'80-100%'),#adjust the x-axis categories
  color_palette = custom_colors,
  xlab = "Radius Range",
  ylab = "Percentage (%)",
  title = " ",
  legend_title = "Domain Types"
)
print(p2)
###then analyze differentially expressed genes / pathways by distance
#####################################################
##

###################################
#compute the minimum vertical distance from all points in one niche to another niche
##X-axis########################################
af.niche234 <- add_metadata_from_csv(
  csv_file = "./vertical_distance niche234.csv",  # CSV file containing cell information
  seurat_obj = af.ST,         # input Seurat object
  metadata_cols = c("vertical_distance",'dist')  # column names to add
)

af.nicheXXX<- subset(af.niche234, niche %in% c('niche4'))
af.nicheXXX<- subset(af.niche234, Celltype %in% c('SC'))

af.nicheXXX <- subset(af.nicheXXX, subset = vertical_distance >= 0 &
                                                  vertical_distance <= 800)
#af.niche124$vertical_distance
# plot the distance curve for cell proportion
p_prop <- plot_cell_composition_by_distance(
  seurat_obj   = af.nicheXXX,
  distance_col = "vertical_distance",
  cell_col     = "Type",
  bin_num      = 6,               
  plot_type    = "proportion",#proportion  Cell
  colors=colors_type
)
print(p_prop)

#table(af.ST$Cell_Anno)
ggsave("niche4_distance_to_lesion_core.pdf",width = 5,height = 3)

# plot the distance curve for a continuous variable / a given factor

af.niche6<- subset(af.ST, Cell_Anno %in% c('Plasma cells'))
colors_Type
# define your preferred colors (assuming Type)
my_colors <- c("CON" = "royalblue", "MP" = "indianred")

# plot

plot_fast_curve(af.nicheXXX, 
                x_var = "vertical_distance",
                y_var = "Ccl8", 
                group_var = "Type",#use_relative = T,
                colors = my_colors)

ggsave("Ccl7_distance_to_plaque_core.pdf",width = 4,height = 2.5)


######################################################################
#use neighborhood distance to distinguish each lesion, and compute the distance from each cell in a niche to the centroid
#Y-axis###########################################################
af.plaque_centroid <- add_metadata_from_csv(
  csv_file = "./centroid_distance.csv",  # CSV file containing cell information
  seurat_obj = af.plaque,         # input Seurat object
  metadata_cols = c("lesion_id", "distance_to_centroid")  # column names to add
)

seu_centroid<- subset(af.plaque_centroid, Type %in% c('CWD'))
seu_centroid<- subset(seu_centroid, subset = domain != "Peri vascular")
#bin distance by percentile
# extract metadata
seu_centroid <- add_distance_percentile(seu_centroid,
                                           dist_col = "min_distance", 
                                           group_col = "group", 
                                           probs = seq(0, 1, #you can also customize: probs = c(0, 0.1, 0.5, 0.9, 1))
                                                       0.2),##bin every 20%
                                           new_col_name = "distance_to_centroid")


#########plot the percentage chart changing with distance
colors_Celltype2<-c('#f3758a','#ff87bbff',"#fbaeb5",#EC
                    '#0F7240', #T
                    '#B7B663',#B
                    '#d5a339',"#f8ff00",'#fccc3a',#Mac
                    "#9d9ac6", #SMC
                    "#bbe379",  "#00a7fa",'#b9e6fb','#354c85',
                    "#b1ddcc", "#C9E1AF", "#31aaaa"#FMC
)

plot_domain_percentage_bar(
  seurat_obj = seu_centroid,
  radius_range_col = "distance_to_centroid",
  domain_col = "Cell_Anno",show_value = F,xlab = "Distance to intima",
  angle = 45, color_palette = colors_Celltype,# 
  x_labels = c("0-20%", "20-40%", "40-60%", "60-80%",'80-100%')
)
ggsave("IWD_celltype_centroid_distance_stacked_bar.pdf",width = 5,height = 3)







################################################################################
####compute 2D coordinate plot########################
#######################################################
library(ggpointdensity)
library(ggpubr)
library(rstatix)
df.distance_x<-read.csv("./vertical_neighborhood.csv",row.names = 1, stringsAsFactors = FALSE)
df.distance_y<-read.csv("./centroid_distance.csv",row.names = 1, stringsAsFactors = FALSE)
common_rows <- intersect(rownames(df.distance_y), rownames(df.distance_x))
# take the intersection of the two data frames
df.distance_y <- df.distance_y[common_rows, ]
df.distance_x <- df.distance_x[common_rows, ]

if(!identical(rownames(df.distance_y), rownames(df.distance_x))) {
  stop("Row names are inconsistent, please check the data")
}
# create a new data frame
df.distance <- data.frame(
  Celltype = df.distance_y$Cell_Anno,
  banksy_niche = df.distance_y$domain,
  vertical_distance = df.distance_x$min_distance,
  distance_to_centroid = df.distance_y$distance_to_centroid,
  Type =df.distance_y$Type,
  group =df.distance_y$group,
  row.names = rownames(df.distance_y)  # keep original row names
)

table(df.distance$Celltype)
########################################
#####plot density map for a given cell
########################################
table(df.distance$Celltype)
x_range <- range(df.distance$vertical_distance, na.rm = TRUE)
y_range <- range(df.distance$distance_to_centroid, na.rm = TRUE)

df.filtered <- df.distance[df.distance$Type == "IWD", ]
df.filtered <- df.filtered[df.filtered$Celltype == "LAMs", ]###change cell
df.filtered <- df.filtered[, c("vertical_distance", "distance_to_centroid")]
##
p <- ggplot(df.filtered, aes(vertical_distance, distance_to_centroid)) + coord_equal()
# Method 1: expand the data range so the density estimate fills the whole area
p.density <- p + 
  stat_density_2d(
    geom = "raster",
    aes(fill = after_stat(density)),
    contour = FALSE,
    n = 300  # increase grid points for a smoother fill
  ) + 
  scale_fill_viridis_c(option = 'viridis') + 
  theme_classic() +
  coord_cartesian(xlim = x_range, ylim = y_range) +
  xlab("Distance to intimal") +  ylab("Distance to plauqe core") +    theme(
    axis.title.x = element_text(size = 16),  # X-axis title font size
    axis.title.y = element_text(size = 16),  # Y-axis title font size
    axis.text.x = element_text(size = 14),   # X-axis tick label font size
    axis.text.y = element_text(size = 14)    # Y-axis tick label font size
  )+ # set Y-axis label
  expand_limits(x = x_range, y = y_range)  # expand limits
#"magma" /"inferno" "plasma" "viridis" "cividis" rocket" "mako" "turbo" 
p.density
####
pdf(file = "LAMs IWD Cell.density.pdf",width =4,height = 3)
p.density
dev.off()


##################################################
##########PCA 2D plot######################################
##################################################
df.PCA <- df.distance[df.distance$Celltype == "FB C4", ]###change cell
df.PCA <- df.PCA[, c('Type',"vertical_distance", "distance_to_centroid")]

table(df.PCA$Type)
# set levels
df.PCA$Type <- factor(df.PCA$Type, levels = c('5w HFD','10w HFD','15w HFD' ))

# plot

p <- ggscatterhist(
  df.PCA,
  x = "vertical_distance",
  y = "distance_to_centroid",
  color = "Type",
  margin.plot = "boxplot",
 # palette = c("#00AFBB", "#E7B800", "#FC4E07"),
  margin.params = list(fill = "Type"),
  ggtheme = theme_classic(),
  ellipse = TRUE  # not shown first, for later modification
)

# 1. compute ANOVA p-value for the X-axis (vertical_distance)
x_anova <- df.PCA %>%
  anova_test(vertical_distance ~ Type) %>%
  as_tibble()
x_anova$p
# 2. compute ANOVA p-value for the Y-axis (distance_to_centroid)  
y_anova <- df.PCA %>%
  anova_test(distance_to_centroid ~ Type) %>%
  as_tibble()
y_anova$p




#################### 
##old distance computation code; new code uses Python-precomputed results
####################
# compute the minimum vertical distance from all points in one niche to another niche
##############
colnames(distance_data)
##needs group as the sample ID
# compute the vertical distance from niche_3 to niche_2
vertical_result <- calculate_grouped_euclidean_distance(
  distance_data,
  source_niche = "niche_3",
  target_niches = c("niche_1", "niche_4", "niche_6")
)

vertical_result <- vertical_result[!is.na(vertical_result$euclidean_distance), ]
# initialize the euclidean_range column
vertical_result$euclidean_range <- NA
# assign bins based on the euclidean_distance value
vertical_result$euclidean_range[vertical_result$euclidean_distance <= 100] <- "0-100"
vertical_result$euclidean_range[vertical_result$euclidean_distance > 100 & vertical_result$euclidean_distance <= 200] <- "100-200"
vertical_result$euclidean_range[vertical_result$euclidean_distance > 200] <- ">200"
# verify the result
table(vertical_result$euclidean_range)

write.csv(vertical_result, "./visiumHD/vertical_distance.csv", row.names = FALSE)
colnames(vertical_result)

af.ST_vertical_distance <- add_metadata_from_csv(
  csv_file = "./visiumHD/vertical_distance.csv",  # CSV file containing cell information
  seurat_obj = af.ST_both,         # input Seurat object
  metadata_cols = c("euclidean_distance",'euclidean_range','banksy_niche')  # column names to add
)

ImageDimPlot(
  af.ST_vertical_distance,
  fov = 'af.ST_15w_2',
  group.by = 'euclidean_range',
  dark.background = F,alpha=0.6,
  border.color='grey',
  cols = colors_ST.FB,
  border.size = 0.1
)
#distance_data is a data frame holding some spatial transcriptomics data. The cell column stores cell barcode information, x_avg stores the cell x coordinate, y_avg stores the cell y coordinate, banksy_niche stores niche information, group stores sample information, and Celltype stores the cell name. I want to run a loop that computes, separately for each group, the distance from every cell with banksy_niche=niche_3 (skip the group if none) to every cell with banksy_niche=niche_1, niche_4 or niche_6, then keep only the rows with banksy_niche=niche_3, output the vertical_distance column as the minimum distance from each cell with banksy_niche=niche_3 to cells with banksy_niche=niche_1, niche_4 or niche_6, then bin it and output the vertical_range column


#########################################################
######compute distance difference####################
#########################################################
# set the folder path containing the CSV files
folder_path <- "./CN_distance"

# get the full paths of all CSV files in the folder
file_list <- list.files(path = folder_path, 
                        pattern = "\\.csv$", 
                        full.names = TRUE)
# initialize an empty list to store the processed data frame of each file
result_list <- list()
# loop over each CSV file
for (file_path in file_list) {
  # read the CSV file
  df <- read.csv(file_path)
  
  # check whether the file contains a vertical_distance column
  if ("vertical_distance" %in% colnames(df)) {
    # extract the vertical_distance column and add the file name as the group identifier
    # use basename to extract the file name without the path
    df_processed <- data.frame(vertical_distance = df$vertical_distance,
                               group = tools::file_path_sans_ext(basename(file_path)))
    
    # append the processed data frame to the result list
    result_list <- append(result_list, list(df_processed))
    
    # print progress information
    print(paste("Successfully processed file:", basename(file_path), "containing", nrow(df), "rows"))
  } else {
    warning(paste("File", basename(file_path), "has no 'vertical_distance' column, skipped"))
  }
}
# merge all data frames in the list into one final data frame
final_df <- do.call(rbind, result_list)
final_df$vertical_distance
table(final_df$group)
colnames(final_df)
my_comparisons1 <- list(c("CN2_to_tumor", "CN4_to_tumor")) 
my_comparisons2 <- list(c("CN5_to_tumor", "CN4_to_tumor"))
my_comparisons3 <- list(c("CN3_to_tumor", "CN4_to_tumor"))


colors<-c("#4991C1",'#c79000','#7AA15E','#93CBAE',
          "#00b1da","#E45D61",'#cdb6ff')

#make violin plot
ggplot(final_df,aes(x=group,y=vertical_distance,fill=group))+ 
  geom_violin(color='black',size=1)+#violin  
  theme_classic() +   theme(text = element_text(size=10, colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15), 
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(color = 'black', size = 12),
        axis.line = element_line(size = 1))+   
  theme(legend.position="none") +  scale_fill_manual(values = colors)+ 
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA)+ #boxplot 
  stat_compare_means(method="t.test",hide.ns = F,  
                   comparisons =c(my_comparisons1,my_comparisons2,my_comparisons3), 
                   label="p.signif",                
                   bracket.size=0.8,   
                   size=6)#add significance comparison


#############################################
###subtract two distances###################
######################################
# 1. extract data of the CN2_to_tumor_distance group
group_tumor <- final_df[final_df$group == "CN5_to_tumor", "vertical_distance"]
# 2. extract data of the CN2_to_CN4_distance group
group_cn4 <- final_df[final_df$group == "CN5_to_CN4", "vertical_distance"]
# 3. verify that the two groups have the same number of rows
if (length(group_tumor) != length(group_cn4)) {
  warning("The two groups have different numbers of rows; matching will be done in order, please ensure the data correspondence is correct.")
}

# 4. compute the difference (CN2_to_tumor_distance - CN2_to_CN4_distance)
vertical_distance_diff <- group_tumor - group_cn4

# 5. create a new data frame df
df <- data.frame(
  vertical_distance = vertical_distance_diff,
  group = "CN3_to_tumor_minor_CN5_to_CN4"
)

df_CN2<-df
df_CN3<-df
df_CN5<-df


############
combined_df3 <- subset(final_df, group == "CN4_to_tumor")
combined_df2 <- rbind(df_CN2, df_CN3)
combined_df1 <- rbind(combined_df2, df_CN5)
combined_df<-rbind(combined_df1, combined_df3)

table(combined_df$group)
my_comparisons1 <- list(c("CN2_to_tumor_minor_CN2_to_CN4", "CN4_to_tumor")) 
my_comparisons2 <- list(c("CN3_to_tumor_minor_CN3_to_CN4", "CN4_to_tumor"))
my_comparisons3 <- list(c("CN3_to_tumor_minor_CN5_to_CN4", "CN4_to_tumor"))

colors<-c('#c79000','#93CBAE','#cdb6ff',"#00b1da")

#make violin plot
ggplot(combined_df,aes(x=group,y=vertical_distance,fill=group))+ 
  geom_violin(color='black',size=1)+#violin  
  theme_classic() +   theme(text = element_text(size=10, colour = "black")) +
  theme(plot.title = element_text(hjust = 0.5, size = 15), 
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(color = 'black', size = 12),
        axis.line = element_line(size = 1))+   
  theme(legend.position="none") +  scale_fill_manual(values = colors)+ 
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) #boxplot 
  stat_compare_means(method="wilcox.test",hide.ns = F,  
                     comparisons =c(my_comparisons1,my_comparisons2,my_comparisons3), 
                     label="p.signif",                
                     bracket.size=0.8,   
                     size=6)#add significance comparison
