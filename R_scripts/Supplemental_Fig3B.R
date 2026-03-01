###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 3B -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "colorRamp2",
  "ComplexHeatmap"
)

for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}


#  Set working Directory
set.seed(42)

#  Load data
load("../Output/Variables/SCP_seurat_liver.R")

#  Unsupervised clustering
seurat_liver <- FindNeighbors(seurat_liver, reduction = "integrated.mnn", dims = 1:30)
seurat_liver <- FindClusters(seurat_liver, resolution = 2.3, cluster.name = "clusters.mnn")

#  Create Heatmap Data
heatmap_data <- table(seurat_liver$clusters.mnn, seurat_liver$sample_labels)
custom_row_order <- c("3","8","6","5","2","9", "7", "4", "0", "1") 
heatmap_data <- heatmap_data[custom_row_order, ]

rownames(heatmap_data) <- 1:nrow(heatmap_data)
colnames(heatmap_data) <- 1:ncol(heatmap_data)

#  Save plot
annotation_col <- data.frame(Group = factor(1:8))
rownames(annotation_col) <- colnames(heatmap_data)
ann_colors <- list(Group = c("1" = "#FDE725", "2" = "#B4DE2C", "3" = "#6DCD59", "4" = "#35B779",
                             "5" = "#31688E", "6" = "#3E4A89", "7" = "#482878", "8" = "#440154"))

heatmap_data_matrix <- as.matrix(heatmap_data)
column_split <- colnames(heatmap_data_matrix)
color_scheme <- colorRampPalette(c("white", "darkred"))(60)
col_fun <- colorRamp2(c(min(heatmap_data_matrix), max(heatmap_data_matrix)), c("white", "darkred"))
ha_col <- HeatmapAnnotation(
  Group = annotation_col$Group,
  col = list(Group = ann_colors$Group),
  show_annotation_name = F,
  annotation_name_gp = gpar(fontsize = 18, fontface = "bold"),
  annotation_name_offset = unit(3, "mm"),
  annotation_name_side = "right",
  annotation_height = unit(6, "mm"),  
  simple_anno_size_adjust = TRUE, 
  show_legend = FALSE
  
)

file_name <- "../Output/Figure/Supplemental Fig3B.png"
png(file_name, width = 8, height = 6, units = "in", res = 300)
Heatmap(
  heatmap_data_matrix, name = "Numbers", col = col_fun,  
  show_column_names = F, show_row_names = T,  
  row_names_gp = gpar(fontsize = 18),  
  cluster_rows = FALSE, cluster_columns = FALSE, 
  column_split = column_split,   
  column_title_gp = gpar(fontsize = 18),  
  top_annotation = ha_col,  
  show_heatmap_legend = TRUE ,
  heatmap_legend_param = list(
    title = "Numbers",     
    title_gp = gpar(fontsize = 16),  
    labels_gp = gpar(fontsize = 14)  
  ),
  column_gap = unit(0.1, "mm"),  
  cell_fun = function(j, i, x, y, width, height, fill) { 
    grid.text(heatmap_data_matrix[i, j], x, y, gp = gpar(fontsize = 16)) 
  }
)
dev.off()



