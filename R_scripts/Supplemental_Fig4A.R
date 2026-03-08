###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 4A -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr", 
  "patchwork",
  "ComplexHeatmap",
  "circlize",
  "ggrepel"
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

heatmap_matrix <- as.matrix(filtered_data)
region_info <- sub(".*_(.)[0-9]+$", "\\1", colnames(heatmap_matrix))  
region_info <- factor(region_info)
mean_by_region <- t(sapply(rownames(heatmap_matrix), function(protein) {
  tapply(heatmap_matrix[protein, ], region_info, mean, na.rm = TRUE)
}))
head(mean_by_region)
mean_by_region_df <- as.data.frame(mean_by_region)
colnames(mean_by_region_df) <- paste0("SCP_", 1:8)
SCP_expression <- mean_by_region_df
merged_df_clean <- SCP_expression[, c("SCP_1", "SCP_2")]
spearman_cor <- cor(merged_df_clean$SCP_1, merged_df_clean$SCP_2, method = "spearman")
cor_text <- paste0("Spearman ρ = ", round(spearman_cor, 3))
num_genes <- nrow(merged_df_clean)
merged_df_clean$Gene <- rownames(merged_df_clean)

#  Save plot
CV_marker <- c("Gulo", "Ces1","Oat","Cpt1a",  "Fads1","Cyp2e1", "Cyp1a2","Acox1",   "Acadl",  "Ugt1a1")
plots<-ggplot(merged_df_clean, aes(x = SCP_1, y = SCP_2)) +
  geom_point(size = 3, color = "grey80") +
  geom_point(data = subset(merged_df_clean, Gene %in% CV_marker),
             aes(x = SCP_1, y = SCP_2), color = "#ff6b6f", size = 4) +
  geom_text_repel(data = subset(merged_df_clean, Gene %in% CV_marker),
                  aes(label = Gene), size = 5, color = "black", max.overlaps = Inf,
                  box.padding = 0.3,       
                  point.padding = 0.3)+       
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",size = 1) +
  scale_x_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)) +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, 2)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_line(color = "grey85", size = 0.2),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1), 
    plot.background = element_rect(fill = "white", color = NA),  
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 16)  
  ) +
  labs(x = "SCP_1", y = "SCP_2") +
  annotate("text", x = 0, y = 6, label = cor_text, size = 6, hjust = 0, color = "black")

ggsave(filename = "../Output/Figure/Supplemental Figure 4A.png", plot = plots, width = 6, height = 6,limitsize = FALSE)











