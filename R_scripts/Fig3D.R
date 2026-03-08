###########################
#### SCP Figure Code ####
###########################

#### -- Figure 3D -- ####

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

#  Plot
plot_single_gene_on_umap <- function(seurat_obj, gene, umap_reduction = "umap.mnn") {
  
  if (!(gene %in% rownames(seurat_obj))) {
    stop(paste("Gene", gene, "not found in the Seurat object."))
  }
  
  gene_expression <- FetchData(seurat_obj, vars = gene)
  umap_coords <- Embeddings(seurat_obj, reduction = umap_reduction)
  tmp_output <- data.frame(
    x = umap_coords[, 1],
    y = umap_coords[, 2],
    Expression = gene_expression[, 1]
  )
  
  p <- ggplot(tmp_output, aes(x = x, y = y, fill = Expression)) +
    geom_point(shape = 21, color = "black", stroke = 0.5, size = 4) +
    scale_fill_gradientn(
      colours = c("#ffffff","#ffead8", "#ff963e", "#ff8000"),
      values = scales::rescale(c(min(tmp_output$Expression),
                                 quantile(tmp_output$Expression, 0.2),
                                 quantile(tmp_output$Expression, 0.9),
                                 max(tmp_output$Expression)))
    )+
    labs( title = gene ,x = "UMAP 1", y = "UMAP 2")+
    theme_minimal() +
    theme(
      plot.title = element_blank(),
      legend.position = "none",
      legend.title = element_blank(),
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.text = element_text(size = 14, color = "black"),
      axis.title = element_text(size = 16)
    )
  p
  return(p)
}

#  UMAP data
genes <- c( "Ass1",  "Ces1", "Cdh1", "Nt5e")

#  Save plot
for (gene in genes) {
  p <- plot_single_gene_on_umap(seurat_liver, gene = gene)
  
  ggsave(
    filename = paste0("../Output/Figure/Fig3D_", gene, ".png"),
    plot = p,
    dpi = 300,
    width = 6,
    height = 5,
    bg = "white"
  )
}
