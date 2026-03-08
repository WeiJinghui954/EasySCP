###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 1I -- ####

#  Prepare Workspace
rm(list = ls()); gc()

packages <- c(
  "dplyr",
  "tibble",
  "circlize",
  "ComplexHeatmap",
  "ggplot2"
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

#  Normalized
scaling_factor = 10000
df_normalized <- protein_cell_matrix %>%
  mutate(across(where(is.numeric), ~ . / sum(.) * scaling_factor)) %>% 
  mutate(across(where(is.numeric), ~ log1p(.)))     

#  PCA
mat <- as.matrix(df_normalized)
mat_t <- t(mat)
mat_t <- mat_t[, apply(mat_t, 2, var, na.rm = TRUE) != 0]
pca_res <- prcomp(mat_t, scale. = TRUE)
pca_plt <- data.frame(
  x = pca_res$x[,1],  # PC1 scores for samples
  y = pca_res$x[,2],  # PC2 scores for samples
  samples = rownames(pca_res$x)  # sample names
)
pca_plt$types <- ifelse(grepl("Immune", pca_plt$samples), 
                        "Immune cell", 
                        "Non-immune cell")

#  Save plot
p <- ggplot(pca_plt, aes(x = x, y = y, color = types, fill = types, label = samples)) +
  geom_point(shape = 21, color = "black", stroke = 0.5, size = 6) +
  scale_fill_manual(values = c("Immune cell" = "#ff8000", "Non-immune cell" = "#1fadc5")) +
  theme_minimal() +
  labs(x = "PC 1", y = "PC 2") +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", size = 1.2),
    axis.text = element_text(size = 18, color = "black"),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 18),
    legend.position = "none"
  )
ggsave(file = paste0("../Output/Figure/Supplemental Fig1I.png"), plot = p, width = 6, height = 6,bg="white")




