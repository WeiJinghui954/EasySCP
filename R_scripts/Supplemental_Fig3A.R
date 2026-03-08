###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 3A -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr"
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

#  Extract UMAP data
umap_data <- DimPlot(
  seurat_liver,
  reduction = "umap.mnn",
  group.by = "sample_labels",
  label = FALSE
)$data
label_mapping <- setNames(as.character(1:8), LETTERS[1:8])
umap_data$sample_labels <- factor(label_mapping[umap_data$sample_labels], levels = as.character(1:8))

#  Save plot
custom_colors <- c("1" = "#FDE725", "2" = "#B4DE2C", "3" = "#6DCD59", "4" = "#35B779",
                   "5" = "#31688E", "6" = "#3E4A89", "7" = "#482878", "8" = "#440154")
p <- ggplot(umap_data, aes(x = sample_labels, y = umapmnn_1, fill = sample_labels)) +
  geom_jitter(shape = 21, width = 0.2, size = 3, alpha = 1, color = "black", stroke = 0.4) +
  geom_smooth(aes(group = 1), method = "loess", color = "black", se = TRUE) +
  labs(x = "Gate", y = "UMAP 1") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

ggsave(filename = "../Output/Figure/Supplemental Fig3A.png", plot = p, dpi = 300, width = 6, height = 5,bg="white")




