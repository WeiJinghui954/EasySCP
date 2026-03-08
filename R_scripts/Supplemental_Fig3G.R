###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 3G -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr", 
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

#  Load SCP data
load("../Output/Variables/SCP_seurat_liver_exp.R")
load("../Output/Variables/SCP_seurat_liver.R")

#  Select target protein
genes<- c("Cyp1a2","Slc1a2", "Cyp2a5","Fmo2", "Aldh3a2", "Gba2",               
          "Cyp2f2","Hsd17b6","Gls2","Tomm40l","Aldh1b1","Gldc")
seurat_genes <- rownames(seurat_liver_exp)
valid_genes <- intersect(genes, seurat_genes)

#  Extract protein expression & UMAP coordinates
seurat_liver_exp <- subset(seurat_liver_exp, sample_dates != "20241220")
gene_expression <- FetchData(seurat_liver_exp, vars = genes)
tmp_output <- data.frame(
  x = as.numeric(seurat_liver@reductions$umap.mnn@cell.embeddings[,1]),
  y = as.numeric(seurat_liver@reductions$umap.mnn@cell.embeddings[,2]),
  gene_expression
)

#  PCA
pca_coords <- Embeddings(seurat_liver, reduction = "pca")
tmp_output <- data.frame(
  x = as.numeric(pca_coords[, 1]),  
  y = as.numeric(pca_coords[, 2]),  
  gene_expression
)

# Plot
plots <- lapply(valid_genes, function(gene) {
  if (gene %in% colnames(tmp_output)) {
    ggplot(tmp_output, aes(x = x, y = y, fill = .data[[gene]])) +
      geom_point(shape = 21, color = "black", stroke = 0.5, size = 4) +
      scale_fill_gradient(
        low = "#ffffff",
        high = "#ff8000",
        name = NULL,
        breaks = c(
          min(tmp_output[[gene]], na.rm = TRUE),
          max(tmp_output[[gene]], na.rm = TRUE)
        ),
        labels = c("Low", "High"),
        guide = guide_colorbar(
          barwidth = 2,
          barheight = 65,
          ticks = FALSE,
          frame.colour = "black",
          frame.linewidth = 0.5
        )
      ) +
      labs(
        title = gene,
        x = "PC 1",
        y = "PC 2"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 30),
        legend.position = "right",
        legend.text = element_text(size = 26),
        panel.grid = element_blank(),
        axis.line = element_line(linewidth = 0.5, color = "black"),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20)
      )
  }
})

plots <- Filter(Negate(is.null), plots)

final_plot <- wrap_plots(plots, ncol = 3, guides = "collect") &
  theme(legend.position = "right")

final_plot

ggsave(filename = "../Output/Figure/Supplemental Fig3G.png", plot = final_plot, width = 18, height = 20,limitsize = FALSE)








