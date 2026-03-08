###########################
#### SCP Figure Code ####
###########################

### -- Supplemental Figure 1J -- ####

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
file_path <- "../Data/MC38_immune_non-immune_cell_Report.tsv"
data <- read.delim(file_path, row.names = NULL)

prot_cols <- grep("Quantity", colnames(data), value = TRUE)
protein_cell_matrix <- data[, prot_cols] %>%
  mutate(across(everything(), ~ as.numeric(as.character(.)))) %>%
  replace(is.na(.), 0)

gene_names <- data$PG.Genes
keep_idx <- !is.na(gene_names) & !duplicated(gene_names)

protein_cell_matrix <- protein_cell_matrix[keep_idx, ]
rownames(protein_cell_matrix) <- gene_names[keep_idx]

# Rename columns
colnames(protein_cell_matrix) <- c(
  paste0("Immune cell ", 1:5),
  paste0("Non-immune cell ", 1:5)
)

# Immune-only correlation
immune_mat <- protein_cell_matrix[, grepl("^Immune", colnames(protein_cell_matrix))]
cor_matrix <- cor(
  immune_mat,
  method = "pearson",
  use = "pairwise.complete.obs"
)

#  Save plot
col_fun <- colorRamp2(
  c(0.2, 0.4, 0.6, 0.8, 1),
  c("#ffdfd4", "#ff9e81", "#ff7b5a", "#ff5232", "#ff0000")
)

file_name <- paste0("../Output/Figure/Supplemental Fig1J.png")
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

#  Save variables
save(protein_cell_matrix,file = "../Output/Variables/Immune_non_immune_data.R")







