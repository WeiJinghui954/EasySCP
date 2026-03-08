###########################
#### SCP Figure Code ####
###########################

#### -- Figure 4A -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "tidyverse",    
  "ComplexHeatmap",
  "circlize",
  "viridis"
)
for (p in packages) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p)
    library(p, character.only = TRUE)
  }
}

#  Set working Directory
set.seed(42)

#  Load SCP data
load("../Output/Variables/SCP_filtered_data.R")
load("../Output/Variables/SCP_sample_info.R")

#  Kruskal-Wallis
heatmap_matrix <- as.matrix(filtered_data)
region_info <- factor(sub(".*_(.)[0-9]+$", "\\1", colnames(heatmap_matrix)))
kruskal_p_values <- apply(heatmap_matrix, 1, function(x) {
  kruskal.test(x ~ region_info)$p.value
})
kruskal_p_values <- data.frame(
  protein = rownames(heatmap_matrix),
  p_value = kruskal_p_values,
  FDR = p.adjust(kruskal_p_values, method = "BH")
)

#  Subset zonated proteins (FDR < 0.05)
zonated_proteins <- subset(kruskal_p_values, FDR < 0.05)
protein_data_filtered <- heatmap_matrix[rownames(heatmap_matrix) %in% zonated_proteins$protein, ]

#  Mean value of each gate
gate_names <- unique(sub(".*_(.)([0-9]+)$", "\\1", colnames(protein_data_filtered)))
gate_names <- sort(gate_names)

merged_data <- sapply(gate_names, function(gate) {
  gate_cols <- grep(paste0(".*_", gate, "[0-9]+$"), colnames(protein_data_filtered))
  rowMeans(protein_data_filtered[, gate_cols], na.rm = TRUE)  
})

#  Normalized to max
merged_data <- as.data.frame(merged_data)
protein_data_normalized <- apply(merged_data, 1, function(protein) {
  return(protein / max(protein))
})
norm_data <- t(protein_data_normalized)

#  COM value
com_values <- apply(norm_data, 1, function(x) {
  zones <- 1:8  
  numerator <- sum(zones * x, na.rm = TRUE)  
  denominator <- sum(x, na.rm = TRUE)  
  com_value <- numerator / denominator  
  return(com_value)
})
com_values <- data.frame(protein = rownames(norm_data), COM = com_values)
ordered_com_values <- com_values[order(com_values$COM), ]  
ordered_proteins <- ordered_com_values$protein
protein_data_sorted_by_com <- norm_data[ordered_proteins, ]

#  Save plot

col_fun <- colorRamp2(
  c(0, 0.2, 0.4, 0.6, 0.8, 1),
  viridis(6, option = "cividis")  # cividis 系列
)

gene <- c("Cyp2a5", "Cyp4a10", "Hacd3", "Nat8f2", "Gulo", "Scd1","Plpp3","Slc22a1","Hsd3b5","Hsd3b7",
          "Clic5", "Commd9", "Grb14", "Gstp1", "Mtres1", "Slc13a2",  "Tmem141", "Rida", "Oxr1", "Ppme1",
          "Cyp2e1", "Gulo", "Cyp2a5","Cyp4a10","Scd1","Plpp3","Cyp2f2","Hsd17b13","Slc13a2",
          "Rabggtb","Ndufb10","Gclc","Slc25a11","Dap3","Prdx1","Gadd45gip1","Hdgfl2","Hsd17b2","Ndufb7",
          "Gba2","Romo1","Elovl5","Dhrs9","Eef2k","Alpk1")
gene <- unique(gene)

genes_in_matrix <- gene[gene %in% rownames(protein_data_sorted_by_com)]
gene_positions <- match(genes_in_matrix, rownames(protein_data_sorted_by_com))

row_anno <- rowAnnotation(
  gene = anno_mark(at = gene_positions, labels = genes_in_matrix,
                   labels_gp = gpar(fontsize = 12),
                   side = "left",
                   link_width = unit(0.5, "cm"))  
)

colnames(protein_data_sorted_by_com) <- 1:8
column_split <- colnames(protein_data_sorted_by_com)

annotation_col <- data.frame(Group = factor(rep(c("1", "2", "3", "4", "5", "6", "7", "8"),
                                              length.out = ncol(protein_data_sorted_by_com))))
rownames(annotation_col) <- colnames(protein_data_sorted_by_com)
ann_colors <- list(Group = c("1" = "#FDE725", "2" = "#B4DE2C", "3" = "#6DCD59", "4" = "#35B779",
                             "5" = "#31688E", "6" = "#3E4A89", "7" = "#482878", "8" = "#440154"))
top_anno <- HeatmapAnnotation(
  df = annotation_col,
  col = ann_colors,
  show_annotation_name = FALSE,  
  annotation_height = unit(5, "mm"),  
  show_legend = FALSE
)

file_name <-  "../Output/Figure/Fig4A.png"
png(file_name, width = 8, height = 10, units = "in", res = 300)
ht <-Heatmap(
  protein_data_sorted_by_com,
  name = "Z-score",
  col = col_fun,
  cluster_rows = FALSE, 
  cluster_columns = FALSE, 
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_split = column_split,   
  column_title_gp = gpar(fontsize = 12),  
  top_annotation = top_anno,  
  
  left_annotation = row_anno, 
  heatmap_legend_param = list(
    title = NULL,
    labels = c("0", "0.2", "0.4", "0.6", "0.8", "1"), 
    at = c(0, 0.2, 0.4, 0.6, 0.8, 1),                
    labels_gp = gpar(fontsize = 12),                 
    legend_height = unit(23, "cm")  
  )
)
draw(ht)
dev.off()

#  Save variables
save(protein_data_sorted_by_com,file = "../Output/Variables/SCP_protein_data_sorted_by_com.R")
save(ordered_com_values,file = "../Output/Variables/SCP_protein_ordered_com_values.R")

  
  










