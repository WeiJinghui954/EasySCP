###########################
#### SCP Figure Code ####
###########################

#### -- Figure 5E -- ####


#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "tidyverse",
  "Seurat",
  "GSVA",
  "ggplot2",
  "circlize",
  "ComplexHeatmap",
  "ggplotify",
  "patchwork"
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

#  Read Bravo_2024 raw data
counts <- read.delim("../Data/Bravo_2024/scrna_counts.tsv", row.names = 1, check.names = FALSE)
cell_types <- sub(".*__([^_]+__.*)$", "\\1", colnames(counts))
#table(cell_types)
cells_liver_sub <- colnames(counts)[cell_types == "922e80__NucSeq_Liver_wild_steady_fresh"]
expr_matrix <- counts[, cells_liver_sub]
expr_matrix_sub <- expr_matrix[rownames(expr_matrix) %in% conserved_zonation_markers, ]

#  HSS score
expr_matrix_sub <- as.matrix(expr_matrix_sub)
params <- gsvaParam(
  expr_matrix_sub,
  HSS_gene_sets, 
  kcdf = "Gaussian" 
)
scores <- gsva(params)
cv_like_score <- scores["CV", ] - scores["PV", ]
sorted_scores <- sort(cv_like_score, decreasing = F)
sorted_scores_matrix <- cbind(Cell = names(sorted_scores),
                              Score = as.numeric(sorted_scores))

#  Save plot
cell_order <- names(sorted_scores)[order(sorted_scores, decreasing = TRUE)]
plot_df <- data.frame(
  Cell = factor(names(sorted_scores), levels = cell_order),  
  Score = as.numeric(sorted_scores)
)

p_points <- ggplot(plot_df, aes(x = Cell, y = Score, color = Score)) +
  geom_point(size = 1) +
  scale_color_gradientn(colors = c("#FDE725", "#B5DE2B", "#28AE80", "#20908C", "#440154")) +
  theme_bw() +
  labs(x = NULL, y = "Score") +
  theme(
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 18),
    axis.text.y = element_text(color = "black", size = 14),
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    plot.background = element_rect(fill = "white", color = NA),
    plot.title = element_text(hjust = 0.5, size = 16),
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    plot.margin = margin(0, 1, 0, 1)
  )
marker_genes <- c("Cyp7a1","Cyp2e1","Tenm3","Pparg","Nxpe2","Gstm2","Axin2","Lgr5","Nt5e","Tbx3os1",               #CV
                  "Slc1a2","Glul", "Ets1","Zeb2","Ptprb","Prkg1","Hgf","Ngf","Pck1","Sds","Cryl1","Hal","Hsd17b13","Cdh1",   
                  "Gls2","Etnppl","Cyp2f2","Sfxn1","Ftcd","Ass1" )   
heatmap_mat <- expr_matrix[marker_genes, cell_order, drop = FALSE]  
heatmap_mat <- t(scale(t(heatmap_mat)))

ht <- Heatmap(
  heatmap_mat,
  name = "Expr",
  col = colorRamp2(c(-2, 0, 2),  c("#007ad4", "white", "firebrick3")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_column_names = FALSE,
  show_row_names = T,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 10),
  show_heatmap_legend = F  )

ht_gg <- as.ggplot(function() draw(
  ht,
  newpage = FALSE,
  padding = unit(c(0, 0, 0, 0), "mm")  ))
final_plot <- p_points / ht_gg + plot_layout(heights = c(1, 1))
ggsave(filename = "../Output/Figure/Fig5E.png", plot = final_plot, width = 6, height = 8,
       dpi = 300,limitsize = FALSE)



