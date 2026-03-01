###########################
#### SCP Figure Code ####
###########################

#### -- Figure 5D -- ####


#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "tidyverse",
  "Seurat",
  "GSVA",
  "ggplot2"
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
load("../Output/Variables/SCP_HSS_gene_sets.R")
load("../Data/Rosenberger_2023/d.R")
load("../Data/Rosenberger_2023/meta_distances.R")
load("../Data/Rosenberger_2023/meta_pg.R")
load("../Data/Rosenberger_2023/SA_incl_all.R")
load("../Data/Rosenberger_2023/d_full_set.R")


#  Rosenberger_2023 data process
SA_incl_heps <- d %>%
  filter(cell_ID %in% meta_distances$cell_ID) %>%
  distinct(cell_ID) %>%
  pull(cell_ID)
d_full_set %>%
  filter(cell_ID %in% SA_incl_heps) %>%
  dplyr::select(cell_ID, int, Protein) %>%
  spread(cell_ID, int) %>%
  column_to_rownames("Protein") -> d_wide_90
protein_symbol_map <- d %>%
  distinct(Protein, Symbol) %>%
  filter(Protein %in% rownames(d_wide_90))
d_wide_90_with_symbol <- d_wide_90 %>%
  rownames_to_column("Protein") %>%
  left_join(protein_symbol_map, by = "Protein")
d_wide_90_with_symbol <- d_wide_90_with_symbol %>%
  filter(!is.na(Symbol))
d_wide_90_avg <- d_wide_90_with_symbol %>%
  dplyr::select(-Protein) %>%
  group_by(Symbol) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE))) %>%
  column_to_rownames("Symbol")
expr_matrix_NM_rawdata <- d_wide_90_avg[rownames(d_wide_90_avg) %in% conserved_zonation_markers, ]

cols_to_remove <- c("m5C_56_target8", "m5C_51_target8",  "m5C_77_target4","m5C_41_target8","m5C_67_target8")
expr_matrix_filtered <- expr_matrix_NM_rawdata[, !(colnames(expr_matrix_NM_rawdata) %in% cols_to_remove)]
valid_cells <- meta_distances %>%
  filter(sum_distance <= 1650) %>%
  pull(cell_ID)
valid_cells_exist <- intersect(valid_cells, colnames(expr_matrix_filtered))
expr_matrix_filtered <- expr_matrix_filtered %>%
  dplyr::select(all_of(valid_cells_exist))

expr_matrix_filtered[is.na(expr_matrix_filtered)] = 0
scaling_factor = 10000
expr_matrix_NM_normalized <- expr_matrix_filtered %>%
  mutate(across(where(is.numeric), ~ . / sum(.) * scaling_factor)) %>% 
  mutate(across(where(is.numeric), ~ log1p(.)))

#  Rosenberger_2023 HSS score
expr_matrix <- as.matrix(expr_matrix_NM_normalized)
params <- gsvaParam(
  expr_matrix,                 
  HSS_gene_sets,       
  kcdf = "Gaussian"                 
)
scores <- gsva(params)
HSS_score <- scores["CV", ] - scores["PV", ]
sorted_scores <- sort(HSS_score, decreasing = F)
sorted_scores_matrix <- cbind(Cell = names(sorted_scores),
                              Score = as.numeric(sorted_scores))

#  Real_distance
cell_distance_subset <- meta_distances[, c("cell_ID","dist_to_central_vein(pixels)","dist_to_portal_vein(pixels)","sum_distance")]
sorted_scores_df <- as.data.frame(sorted_scores_matrix, stringsAsFactors = FALSE)
merged_df <- merge(cell_distance_subset, sorted_scores_df,by.x = "cell_ID", by.y = "Cell")
merged_df$Score <- as.numeric(merged_df$Score)
merged_df$Ratio_CV_PV <- merged_df$`dist_to_portal_vein(pixels)`/merged_df$sum_distance
merged_df$Score <- as.numeric(merged_df$Score)
merged_df$Ratio_CV_PV <- as.numeric(merged_df$Ratio_CV_PV)

#  Spearman
cor_spearman <- cor.test(merged_df$Ratio_CV_PV, merged_df$Score, method = "spearman")
spearman_r <- round(cor_spearman$estimate, 3)
spearman_p <- signif(cor_spearman$p.value, 3)

#  Save plot
fit <- lm(Score ~ Ratio_CV_PV, data = merged_df)
newdata <- data.frame(Ratio_CV_PV = seq(0,1,length.out=200))
pred_ci <- predict(fit, newdata, se.fit = TRUE)
newdata$fit <- pred_ci$fit
newdata$se <- pred_ci$se.fit
newdata$ci_lower <- newdata$fit - qt(0.975, df=fit$df.residual) * newdata$se
newdata$ci_upper <- newdata$fit + qt(0.975, df=fit$df.residual) * newdata$se
sigma <- sqrt(sum(residuals(fit)^2) / fit$df.residual)
newdata$pi_lower <- newdata$fit - qt(0.975, df=fit$df.residual) * sqrt(newdata$se^2 + sigma^2)
newdata$pi_upper <- newdata$fit + qt(0.975, df=fit$df.residual) * sqrt(newdata$se^2 + sigma^2)
merged_df$in_PI <- ifelse(merged_df$Score >= newdata$pi_lower[findInterval(merged_df$Ratio_CV_PV, newdata$Ratio_CV_PV)] &
                            merged_df$Score <= newdata$pi_upper[findInterval(merged_df$Ratio_CV_PV, newdata$Ratio_CV_PV)],
                          TRUE, FALSE)

#  Save plot
p <- ggplot(merged_df, aes(x = Ratio_CV_PV, y = Score)) +
  geom_point(aes(color = in_PI), alpha = 0.6, size = 3) +
  scale_color_manual(values = c("TRUE" = "grey10", "FALSE" = "grey60"), guide = "none") +
  geom_ribbon(data=newdata,
              aes(x=Ratio_CV_PV, ymin=ci_lower, ymax=ci_upper),
              fill="#f00089", alpha=0.4, inherit.aes = FALSE) +
  geom_ribbon(data=newdata,
              aes(x=Ratio_CV_PV, ymin=pi_lower, ymax=pi_upper),
              fill="#f7d417", alpha=0.2, inherit.aes = FALSE) +
  geom_line(data=newdata,
            aes(x=Ratio_CV_PV, y=fit),
            color="#dc143c", size=1, inherit.aes = FALSE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(-2, 2), expand = FALSE) +
  theme_minimal(base_size = 14)  +
  theme(
    panel.grid.major = element_line(color = "grey85", size = 0.5),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  annotate("text",
           x = min(merged_df$Ratio_CV_PV, na.rm = TRUE),
           y = 1.9,
           hjust = 0, vjust = 1,
           label = paste0("Spearman ρ = ", spearman_r, ", p = ", spearman_p),
           size = 5) +
  labs( x = "Relative distance to PV")

ggsave(filename = "../Output/Figure/Fig5D.png", plot = p, width = 6, height = 6,bg = "white",limitsize = FALSE)







