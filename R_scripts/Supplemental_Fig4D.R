###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 4D -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "ggplot2",
  "dplyr",
  "tibble",
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
mean_by_region_df <- as.data.frame(mean_by_region)
colnames(mean_by_region_df) <- paste0("SCP_", 1:8)
SCP_expression <- mean_by_region_df
SCP_expression_CV_PV<- matrix(ncol = 2, nrow = nrow(SCP_expression))

# CV PV area
SCP_expression_CV_PV[, 1] <- rowMeans(SCP_expression[, c("SCP_1", "SCP_2")])
SCP_expression_CV_PV[, 2] <- rowMeans(SCP_expression[, c("SCP_7", "SCP_8")])
colnames(SCP_expression_CV_PV) <- c("SCP_CV", "SCP_PV")
rownames(SCP_expression_CV_PV) <- rownames(SCP_expression)

#  Load Xu_2024 data (ST stereo)
load("../Data/Xu_2024/ST_stereo_expression_group.RData")
ST_stereo_expression_group
colnames(ST_stereo_expression_group) <- paste0("ST_stereo_", 1:9)
ST_stereo_expression_CV_PV<- matrix(ncol = 2, nrow = nrow(ST_stereo_expression_group))

# CV PV area
ST_stereo_expression_CV_PV[, 1] <- rowMeans(ST_stereo_expression_group[, c("ST_stereo_1", "ST_stereo_2")])
ST_stereo_expression_CV_PV[, 2] <- rowMeans(ST_stereo_expression_group[, c("ST_stereo_8", "ST_stereo_9")])

colnames(ST_stereo_expression_CV_PV) <- c("ST_stereo_CV", "ST_stereo_PV")
rownames(ST_stereo_expression_CV_PV) <- rownames(ST_stereo_expression_group)


df_list <- list(SCP_expression_CV_PV, ST_stereo_expression_CV_PV) %>%
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
merged_df$ST_stereo_bias <- log2(merged_df$ST_stereo_CV / merged_df$ST_stereo_PV)
merged_df_clean_All <- merged_df[complete.cases(merged_df$ST_stereo_bias), ]

spearman_cor <- cor(merged_df_clean_All$SCP_bias, merged_df_clean_All$ST_stereo_bias, method = "spearman")
cor_text1 <- paste0("Spearman ρ = ", round(spearman_cor, 3))


#  Load SCP zonated data
load("../Output/Variables/SCP_protein_data_sorted_by_com.R")
SCP_zonated_protein <- rownames(protein_data_sorted_by_com) 

#  Load Xu_2024 (ST stereo) zonated data
load("../Data/Xu_2024/ST_Stereoseq_Zonated.RData")
intersection_SCP_STstereo <- Reduce(intersect, list(
  SCP_zonated_protein,
  ST_Stereoseq_union_zonated_gene
))


merged_df_clean_zonated <- merged_df_clean_All[rownames(merged_df_clean_All) %in% intersection_SCP_STstereo, ]
spearman_cor <- cor(merged_df_clean_zonated$SCP_bias, merged_df_clean_zonated$ST_stereo_bias, method = "spearman")
cor_text2 <- paste0("Spearman ρ = ", round(spearman_cor, 3))

#  Save plot
merged_df_clean_All$Gene <- rownames(merged_df_clean_All)
merged_df_clean_zonated$Gene <- rownames(merged_df_clean_zonated)

highlight_PR_Down <-c("Glul", "Gulo","Cyp2e1","Rhbg", "Slc13a3", "Oat")              
highlight_PR_Up <-c("Ftcd","Cyp2f2","Cryl1","Tkfc","Rida","Hunk" )                 
highlight_P_Up_R_Down <-c("Slc25a19","Plcb1","Psmf1")                  
highlight_P_Down_R_Up <-c("Ctsz","Pdhb","Preb")  

plots <- ggplot(merged_df_clean_All, aes(x = SCP_bias, y = ST_stereo_bias)) +
  geom_hline(yintercept = 0, color = "grey60", size = 1) +
  geom_vline(xintercept = 0, color = "grey60", size = 1) +
  geom_point(size = 3, alpha = 1, color = "grey70")+
  geom_point(data = subset(merged_df_clean_zonated, Gene %in% intersection_SCP_STstereo ),
             aes(x = SCP_bias, y = ST_stereo_bias), color = "grey50", size = 3) +
  geom_point(data = subset(merged_df_clean_All, Gene %in% highlight_PR_Up),
             aes(x = SCP_bias, y = ST_stereo_bias), color = "#e53238", size = 3) +
  geom_text_repel(data = subset(merged_df_clean_All, Gene %in% highlight_PR_Up),
                  aes(label = Gene), size = 6, color = "#e53238")  +
  geom_point(data = subset(merged_df_clean_All, Gene %in% highlight_PR_Down),
             aes(x = SCP_bias, y = ST_stereo_bias), color = "#0064d2", size = 3) +
  geom_text_repel(data = subset(merged_df_clean_All, Gene %in% highlight_PR_Down),
                  aes(label = Gene), size = 6, color = "#0064d2")+
  geom_point(data = subset(merged_df_clean_All, Gene %in% highlight_P_Up_R_Down),
             aes(x = SCP_bias, y = ST_stereo_bias), color = "#f5af02", size = 3) +
  geom_text_repel(data = subset(merged_df_clean_All, Gene %in% highlight_P_Up_R_Down),
                  aes(label = Gene), size = 6, color = "#f5af02")+
  geom_point(data = subset(merged_df_clean_All, Gene %in% highlight_P_Down_R_Up),
             aes(x = SCP_bias, y = ST_stereo_bias), color = "#86b817", size = 3) +
  geom_text_repel(data = subset(merged_df_clean_All, Gene %in% highlight_P_Down_R_Up),
                  aes(label = Gene), size = 6, color = "#86b817")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 2)) +
  scale_y_continuous(limits = c(-10, 10), breaks = seq(-10, 10, 2)) +
  theme_minimal(base_size = 18) +
  theme(
    panel.grid.major = element_line(color = "grey85", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  labs(x = "SCP bias", y = "ST_Stereo bias") +
  annotate("text", x = -10, y = 10, label = cor_text1, size = 6, hjust = 0, color = "grey70")+
  annotate("text", x = -10, y = 9, label = cor_text2, size = 6, hjust = 0, color = "grey50")

ggsave(filename = "../Output/Figure/Supplemental Figure 4D.png",plot = plots, width = 6, height = 6,bg="white")

