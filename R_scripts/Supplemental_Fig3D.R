###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Fig3D -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "ggplot2",
  "Seurat",
  "scales"
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
load("../Output/Variables/SCP_seurat_liver.R")

#  Plot PCA
plot_single_gene_on_pca <- function(seurat_obj, gene, pca_reduction = "pca") {
  
  if (!(gene %in% rownames(seurat_obj))) {
    stop(paste("Gene", gene, "not found in the Seurat object."))
  }
  
  gene_expression <- FetchData(seurat_obj, vars = gene)
  pca_coords <- Embeddings(seurat_obj, reduction = pca_reduction)
  
  tmp_output <- data.frame(
    x = pca_coords[, 1],   # PC1
    y = pca_coords[, 2],   # PC2
    Expression = gene_expression[, 1]
  )
  
  p <- ggplot(tmp_output, aes(x = x, y = y, fill = Expression)) +
    geom_point(shape = 21, color = "black", stroke = 0.5, size = 4) +
    scale_fill_gradientn(
      colours = c("#ffffff", "#ffead8", "#ff963e", "#ff8000"),
      values = scales::rescale(c(
        min(tmp_output$Expression),
        quantile(tmp_output$Expression, 0.2),
        quantile(tmp_output$Expression, 0.9),
        max(tmp_output$Expression)
      ))
    ) +
    labs(
      title = gene,
      x = "PC 1",
      y = "PC 2"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 30),
      legend.position = "none",
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 16)
    )
  
  return(p)
}

#  UMAP data
genes <- c( "Ass1",  "Ces1", "Cdh1", "Nt5e")

#  Save plot
for (gene in genes) {
  p <- plot_single_gene_on_pca(seurat_liver, gene = gene)
  
  ggsave(
    filename = paste0("../Output/Figure/Supplemental Fig3D_", gene, ".png"),
    plot = p,
    dpi = 300,
    width = 6,
    height = 5,
    bg = "white"
  )
}

