###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 1K -- ####

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

#  Load data
load("../Output/Variables/Immune_non_immune_data.R")

# Immune-only correlation
non_immune_mat <- protein_cell_matrix[, grepl("^Non-immune", colnames(protein_cell_matrix))]
cor_matrix <- cor(
  non_immune_mat,
  method = "pearson",
  use = "pairwise.complete.obs"
)

#  Save plot
col_fun <- colorRamp2(
  c(0.6,0.7, 0.8, 0.9,  1),
  c("#ffdfd4",  "#ff9e81", "#ff7b5a", "#ff5232", "#ff0000")
)

file_name <- paste0("../Output/Figure/Supplemental Fig1K.png")
png(file_name, width = 8, height = 8, units = "in", res = 300)
Heatmap(
  cor_matrix,
  name = "Pearson Correlation",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_gp = gpar(fontsize = 16),
  row_names_gp = gpar(fontsize = 16),
  heatmap_legend_param = list(
    title = NULL,
    labels_gp = gpar(fontsize = 14)
  ),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid.text(sprintf("%.2f", cor_matrix[i, j]), x, y,
              gp = gpar(fontsize = 12))
  },
  width = unit(10, "cm"),
  height = unit(10, "cm")
)
dev.off()




