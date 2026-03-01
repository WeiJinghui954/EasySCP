###########################
#### SCP Figure Code ####
###########################

#### -- Figure 5A -- ####


#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "tidyverse",       
  "clusterProfiler", 
  "org.Mm.eg.db",
  "Seurat",
  "batchelor",
  "SeuratWrappers",
  "ComplexHeatmap",
  "circlize"
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
load("../Output/Variables/Conserved_zonation_markers.R")
load("../Output/Variables/SCP_sample_info.R")

#  Read SCP liver raw data
file_path <- "../Data/SCP_Liver_Report.tsv"
dat1 <- read.delim(file_path, row.names = NULL)
dat1 <- dat1[dat1$PG.Genes != "", ] 
dat1$PG.Genes <- make.unique(dat1$PG.Genes) 
prot_cols <- grep("Quantity", colnames(dat1), value = TRUE)
prot_data <- dat1[, prot_cols] 
colnames(prot_data) <- sub(".*(\\d{8}).*(plate\\d+[_\\.][A-Za-z0-9]+).*", "\\1_\\2", colnames(prot_data))
colnames(prot_data) <- gsub("\\.", "_", colnames(prot_data))
prot_data <- as.data.frame(lapply(prot_data, function(x) as.numeric(as.character(x))))
prot_data[is.na(prot_data)] = 0
rownames(prot_data) = dat1$PG.Genes

# 458 core protein
prot_data_filtered <- prot_data[rownames(prot_data) %in% conserved_zonation_markers, ]

#  Normalization
scaling_factor = 10000
df_normalized <- prot_data_filtered %>%
  mutate(across(where(is.numeric), ~ . / sum(.) * scaling_factor)) %>% 
  mutate(across(where(is.numeric), ~ log1p(.)))    

#  Create a Seurat object
seurat_liver <- CreateSeuratObject(data =as.matrix(df_normalized), counts = as.matrix(df_normalized))
seurat_liver$sample_labels = sample_info$labels
seurat_liver$dates = sample_info$dates
seurat_liver <- subset(seurat_liver, dates != "20241220")
seurat_liver[["RNA"]] <- split(seurat_liver[["RNA"]], f = seurat_liver$dates)
seurat_liver <- ScaleData(seurat_liver)
seurat_liver <- FindVariableFeatures(seurat_liver, nfeatures = 458)
seurat_liver <- RunPCA(seurat_liver, features = VariableFeatures(seurat_liver))

#  Batch Integration (FastMNN)
seurat_liver <- IntegrateLayers(
  object = seurat_liver, 
  method = FastMNNIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.mnn",
  k = 85,      
  ndist = 1.14,    
  verbose = FALSE
)

#  Unsupervised clustering
seurat_liver <- RunPCA(seurat_liver, reduction = "integrated.mnn", reduction.name = "pca.integrated")
seurat_liver <- FindNeighbors(seurat_liver, reduction = "integrated.mnn", dims = 1:30)
seurat_liver <- FindClusters(seurat_liver, resolution = 2.5, cluster.name = "clusters.mnn")

#  Create Heatmap Data
heatmap_data <- table(seurat_liver$clusters.mnn, seurat_liver$sample_labels)
custom_row_order <- c("0", "9", "2", "6", "7", "5", "3", "4", "8", "1") 
heatmap_data <- heatmap_data[custom_row_order, ]

rownames(heatmap_data) <- 1:nrow(heatmap_data)
colnames(heatmap_data) <- 1:ncol(heatmap_data)

#  Save plot
annotation_col <- data.frame(Group = factor(rep(c("1", "2", "3", "4", "5", "6", "7", "8"), length.out = ncol(heatmap_data))))
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

file_name <- "../Output/Figure/Fig5A.png"
png(file_name, width = 8, height = 6, units = "in", res = 300)
Heatmap(
  heatmap_data_matrix,  
  name = "Numbers",  
  col = col_fun,  
  show_column_names = F,  
  show_row_names = T,  
  row_names_gp = gpar(fontsize = 18),  
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
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

#  Save variables
save(seurat_liver,file = "../Output/Variables/SCP_seurat_liver_conserved_zonation_markers.R")
  
  
