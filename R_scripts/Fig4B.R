###########################
#### SCP Figure Code ####
###########################

#### -- Figure 4B -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "ggplot2",
  "dplyr",
  "tidyr",      
  "Seurat",
  "patchwork",
  "scales",
  "tibble"      
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
load("../Output/Variables/SCP_sample_info.R")


seurat_liver_exp <- CreateSeuratObject(data =as.matrix(filtered_data), counts = as.matrix(filtered_data))
seurat_liver_exp <- ScaleData(seurat_liver_exp)
seurat_liver_exp$sample_labels = sample_info$labels
seurat_liver_exp$sample_dates = sample_info$dates

  
#  Subset protein
genes_to_keep<- c("Gba2","Cyp4a10","Plpp3",
                  "Romo1","Elovl5", "Slc22a1",                
                  "Hsd3b5","Scd1", "Hacd3",
                  
                  "Gstp1","Slc13a2","Gclc",
                  "Mtres1","Hsd17b13","Alpk1",
                  "Eef2k","Dhrs9", "Clic5"   
) 

#  Protein expression (normalized to mean)
gene_expression_long <- FetchData(seurat_liver_exp, vars = genes_to_keep) %>%
  rownames_to_column("Cell") %>%
  pivot_longer(cols = -Cell, names_to = "Gene", values_to = "Expression") %>%
  mutate(Group = seurat_liver_exp$sample_labels[match(Cell, rownames(seurat_liver_exp@meta.data))])
gene_expression_long$Gene <- factor(gene_expression_long$Gene, levels = genes_to_keep)

gene_means <- gene_expression_long %>%
  group_by(Gene) %>%
  summarise(GeneMean = mean(Expression, na.rm = TRUE))
gene_expression_long <- gene_expression_long %>%
  left_join(gene_means, by = "Gene") %>%
  mutate(Expression_norm = Expression / GeneMean)

summary_data_SCP <- gene_expression_long %>%
  group_by(Group, Gene) %>%
  summarise(
    Mean = mean(Expression_norm, na.rm = TRUE),
    SEM  = sd(Expression_norm, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    Group = factor(Group, levels = LETTERS[1:8], labels = as.character(1:8)),
    Source = "SCP"
  )

#  Plot style
source_colors <- c("SCP" = "#e50914", "ST 10X" = "#838f97", "ST Stereo" = "#ccd7dd")
source_linetypes <- c("SCP" = "solid", "ST 10X" = "solid", "ST Stereo" = "solid")
source_sizes <- c("SCP" = 1.5, "ST 10X" = 1, "ST Stereo" = 1)

#  Save plot
plots <- lapply(unique(summary_data_SCP$Gene), function(gene) {
  df_gene <- filter(summary_data_SCP, Gene == gene)
  df_gene$Group <- as.numeric(as.character(df_gene$Group))
  
  ggplot(df_gene, aes(x = Group, y = Mean, color = Source, group = Source)) +
    geom_ribbon(aes(ymin = Mean - SEM, ymax = Mean + SEM, fill = Source), alpha = 0.3, color = NA) +
    geom_line(aes(linetype = Source, size = Source)) +
    scale_color_manual(values = source_colors) +
    scale_fill_manual(values = source_colors) +
    scale_linetype_manual(values = source_linetypes) +
    scale_size_manual(values = source_sizes) +
    scale_y_continuous(
      name = "Protein expression\n(normalized to mean)",
      breaks = scales::pretty_breaks(n = 3),
      labels = scales::number_format(accuracy = 0.1)
    ) +
    scale_x_continuous(breaks = 1:8, expand = c(0, 0)) +
    labs(title = gene) +
    theme_classic() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 26),
      legend.position = "none",
      axis.text    = element_text(size = 16),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      plot.margin  = margin(t = 5, b = 5, r = 30, l = 30)
    )
})

ncol <- 3
nrow <- ceiling(length(plots) / ncol)
combined_plot <- wrap_plots(plots, ncol = ncol)

ggsave(
  filename = "../Output/Figure/Fig4B.png",
  plot = combined_plot, width = 14, height = nrow * 4, limitsize = FALSE
)

#  Save variables
save(seurat_liver_exp,file = "../Output/Variables/SCP_seurat_liver_exp.R")
