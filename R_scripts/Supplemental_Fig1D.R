###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 1D -- ####

#  Prepare Workspace
rm(list = ls()); gc()

packages <- c(
  "dplyr",
  "tibble",
  "circlize",
  "ComplexHeatmap"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
  library(p, character.only = TRUE)
}
#  Set working Directory
set.seed(42)

#  Read HEK293T 50 cells raw data
file_path <- "../Data/HEK293T_50cells_Report.tsv"
data <- read.delim(file_path, row.names = NULL)
data_1 <- data %>%
  select(`PG.Genes`, matches("plate1")) %>%
  filter(!is.na(`PG.Genes`)) %>%
  mutate(`PG.Genes` = make.unique(as.character(`PG.Genes`))) %>%
  column_to_rownames(var = "PG.Genes")

colnames(data_1) <- sub(".*_", "", colnames(data_1))
data_1 <- data_1[, !grepl("IBAQ", colnames(data_1))]
colnames(data_1) <- sub("\\.raw\\.PG\\.Quantity$", "", colnames(data_1))

data_1_numeric <- data_1 %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))

row_sums <- rowSums(data_1_numeric != 0, na.rm = TRUE)
data_1_numeric <- data_1_numeric[row_sums > 0, ]

colnames(data_1_numeric) <- paste0(
  "HEK293T 50 cells ", seq_len(ncol(data_1_numeric))
)

cor_matrix <- cor(data_1_numeric, method = "pearson", use = "pairwise.complete.obs")

#  Save plot
col_fun <- colorRamp2(
  c(0.80,0.84, 0.88, 0.92, 0.96, 1),
  c("#ffdfd4", "#ffbfaa", "#ff9e81", "#ff7b5a", "#ff5232", "#ff0000")
)
file_name <- paste0("../Output/Figure/SupplemeDntal Fig1D.png")
png(file_name, width = 9, height = 9, units = "in", res = 300)
Heatmap(
  cor_matrix,
  name = "Pearson Correlation",
  show_column_names = T,  
  show_row_names = T,     
  col = col_fun,               
  row_title = NULL,  
  cluster_rows = F,  
  cluster_columns = F,  
  column_gap = unit(0.2, "mm"),  
  row_gap = unit(0.2, "mm"),     
  column_names_gp = gpar(fontsize = 16),  
  row_names_gp = gpar(fontsize = 16),     
  heatmap_legend_param = list(
    title = "Pearson Correlation",        
    title_gp = gpar(fontsize = 14),       
    labels_gp = gpar(fontsize = 14),      
    legend_height = unit(4, "cm"),        
    legend_width = unit(1, "cm")          
  ),  
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y, gp = gpar(fontsize = 12))
  },
  height =unit(12, "cm"),  
  width = unit(12, "cm")
)
dev.off()

#  Save variables
save(data_1_numeric,file = "../Output/Variables/HEK293T_50cells_data.R")












