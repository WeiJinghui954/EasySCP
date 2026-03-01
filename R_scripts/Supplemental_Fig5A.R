###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 5A -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "ggplot2",
  "dplyr",
  "tibble",
  "ggrepel",
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
heatmap_matrix <- as.matrix(filtered_data)
region_info <- sub(".*_(.)[0-9]+$", "\\1", colnames(heatmap_matrix))  
region_info <- factor(region_info)
mean_by_region <- t(sapply(rownames(heatmap_matrix), function(protein) {
  tapply(heatmap_matrix[protein, ], region_info, mean, na.rm = TRUE)
}))
mean_by_region_df <- as.data.frame(mean_by_region)
colnames(mean_by_region_df) <- paste0("SCP_", 1:8)
SCP_expression <- mean_by_region_df
SCP_expression_CV_PV<- matrix(ncol = 2, nrow = nrow(SCP_expression))

# CV PV area
SCP_expression_CV_PV[, 1] <- rowMeans(SCP_expression[, c("SCP_1", "SCP_2")])
SCP_expression_CV_PV[, 2] <- rowMeans(SCP_expression[, c("SCP_7", "SCP_8")])
colnames(SCP_expression_CV_PV) <- c("SCP_CV", "SCP_PV")
rownames(SCP_expression_CV_PV) <- rownames(SCP_expression)

#  Load Guilliams_2022 data (ST 10X Visium)
load("../Data/Guilliams_2022/ST_Visium_expression_group.RData")
ST_10X_expression_group <- as.data.frame(ST_10X_expression_group)
colnames(ST_10X_expression_group) <- paste0("ST_10X_", 1:8)
ST_10X_expression_CV_PV<- matrix(ncol = 2, nrow = nrow(ST_10X_expression_group))

# CV PV area
ST_10X_expression_CV_PV[, 1] <- rowMeans(ST_10X_expression_group[, c("ST_10X_1", "ST_10X_2")])
ST_10X_expression_CV_PV[, 2] <- rowMeans(ST_10X_expression_group[, c("ST_10X_7", "ST_10X_8")])

colnames(ST_10X_expression_CV_PV) <- c("ST_10X_CV", "ST_10X_PV")
rownames(ST_10X_expression_CV_PV) <- rownames(ST_10X_expression_group)

df_list <- list(SCP_expression_CV_PV, ST_10X_expression_CV_PV) %>%
  lapply(function(df) {
    df <- as.data.frame(df)
    df %>% rownames_to_column("Symbol")
  })
merged_df <- Reduce(function(x, y) {
  merge(x, y, by = "Symbol", all = FALSE)
}, df_list)
rownames(merged_df) <- merged_df$Symbol
merged_df$Symbol <- NULL

#  Bias
merged_df$SCP_bias <- log2(merged_df$SCP_CV / merged_df$SCP_PV) 
merged_df$ST_10X_bias <- log2(merged_df$ST_10X_CV / merged_df$ST_10X_PV)
merged_df_clean_All <- merged_df[complete.cases(merged_df$ST_10X_bias), ]

spearman_cor <- cor(merged_df_clean_All$SCP_bias, merged_df_clean_All$ST_10X_bias, method = "spearman")
cor_text1 <- paste0("Spearman ρ = ", round(spearman_cor, 3))

#  Load SCP zonated data
load("../Output/Variables/SCP_protein_data_sorted_by_com.R")
SCP_zonated_protein <- rownames(protein_data_sorted_by_com) 

#  Load Guilliams_2022 (ST 10X Visium) zonated data
load("../Data/Guilliams_2022/ST_Visium_Zonated.RData")
intersection_SCP_ST10X <- Reduce(intersect, list(
  SCP_zonated_protein,
  ST_Visium_union_zonated_gene
))

df_zonated_only <- merged_df_clean_All[
  rownames(merged_df_clean_All) %in% intersection_SCP_ST10X, 
]
df_zonated_only$Group <- "Same_direction"

df_zonated_only$Group[
  df_zonated_only$SCP_bias > 0 &
    df_zonated_only$ST_10X_bias < 0
] <- "SCP_up_ST_down"
df_zonated_only$Group[
  df_zonated_only$SCP_bias < 0 &
    df_zonated_only$ST_10X_bias > 0
] <- "SCP_down_ST_up"

plot_zonated_only <- ggplot(df_zonated_only, 
                            aes(x = SCP_bias, y = ST_10X_bias)) +
  
  geom_hline(yintercept = 0, color = "grey70", size = 1) +
  geom_vline(xintercept = 0, color = "grey70", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  
  geom_point(aes(color = Group), size = 3) +
  
  scale_color_manual(values = c(
    "SCP_up_ST_down" = "#e53238",
    "SCP_down_ST_up" = "#0064d2",
    "Same_direction" = "grey60"
  )) +
  
  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 2)) +
  scale_y_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 2)) +
  
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_line(color = "grey85", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"   
  ) +
  
  labs(x = "SCP bias", 
       y = "Guilliams et al. (2022) bias")
plot_zonated_only

ggsave(filename = "../Output/Figure/Supplemental Figure 5A.png",plot = plot_zonated_only, width = 6, height = 6,bg="white")

#  Save variables
save(df_zonated_only,file = "../Output/Variables/SCP_Guilliams_gene_protein.R")

