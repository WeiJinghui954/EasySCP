###########################
#### SCP Figure Code ####
###########################

#### -- Figure 5B -- ####


#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "tidyverse",
  "Seurat"
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
load("../Output/Variables/SCP_seurat_liver_conserved_zonation_markers.R")

#  Screening cells
cells_CV <- colnames(seurat_liver)[
  seurat_liver$sample_labels %in% "A" &
    as.character(seurat_liver$clusters.mnn) %in% "0"]
cells_PV <- colnames(seurat_liver)[
  seurat_liver$sample_labels %in%  "H" &
    as.character(seurat_liver$clusters.mnn) %in% "1"]
cells_use <- c(cells_CV, cells_PV)

seurat_liver <- JoinLayers(seurat_liver)

sub_obj <- subset(seurat_liver, cells = cells_use)
sub_obj$zone_label <- NA
sub_obj$zone_label[colnames(sub_obj) %in% cells_CV] <- "CV"
sub_obj$zone_label[colnames(sub_obj) %in% cells_PV] <- "PV"
sub_obj$zone_label <- factor(sub_obj$zone_label, levels = c("CV", "PV"))

#  Markers
Idents(sub_obj) <- "zone_label"
markers <- FindMarkers(sub_obj,
  ident.1 = "PV",ident.2 = "CV",
  logfc.threshold = 0, layer = "integrated.mnn")
markers$change <- ifelse(markers$p_val < 0.05 & markers$avg_log2FC > 0.5, "PV", 
                         ifelse(markers$p_val < 0.05 & markers$avg_log2FC < -0.5, "CV", 
                                "No_Change"))
top5_CV <- markers %>%
  filter(change == "CV") %>%            
  arrange(avg_log2FC) %>%         
  slice_head(n = 5)                  
top5_PV <- markers %>%
  filter(change == "PV") %>%             
  arrange(desc(avg_log2FC)) %>%         
  slice_head(n = 5)                    
genes_to_keep <- c(rownames(top5_CV), rownames(top5_PV))

CV_genes <- rownames(markers[markers$change == "CV", ])
PV_genes <- rownames(markers[markers$change == "PV", ])
HSS_gene_sets <- list(
  CV = CV_genes,
  PV = PV_genes)

#  Expression quantity extraction
gene_expression <- FetchData(seurat_liver, vars = genes_to_keep)
cluster_rename <- c(
  "Cluster 1" = "0",
  "Cluster 10" = "1",
  "Cluster 3" = "2",
  "Cluster 7" = "3",
  "Cluster 8" = "4",
  "Cluster 6" = "5",
  "Cluster 4" = "6",
  "Cluster 5" = "7",
  "Cluster 9" = "8",
  "Cluster 2" = "9")
seurat_liver$clusters_renamed <- plyr::mapvalues(
  x = seurat_liver$seurat_clusters,
  from = cluster_rename,
  to = names(cluster_rename))
gene_expression$Group <- seurat_liver$clusters_renamed
expr_long <- gene_expression %>%
  pivot_longer(
    cols = all_of(genes_to_keep),
    names_to = "Gene",
    values_to = "Expression")
expr_summary <- expr_long %>%
  group_by(Group, Gene) %>%
  summarise(mean_exp = mean(Expression), .groups = "drop")
expr_summary <- expr_summary %>%
  group_by(Gene) %>%
  mutate(ratio = mean_exp / sum(mean_exp)) %>%
  ungroup()
expr_summary$Group <- factor(expr_summary$Group, 
                             levels = c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5",
                                        "Cluster 6","Cluster 7","Cluster 8","Cluster 9","Cluster 10"))

#  Save plot
custom_colors <- c(
  "Cyp2e1"  = "#d4be02",
  "Cyp2c29" = "#edd502",
  "Gulo"    = "#fde40c",
  "Cyp1a2"  = "#fde725",
  "Ftl1"    = "#fdea3e",
  
  "Cyp2f2"  = "#6d0286",
  "Gls2"    = "#440154",
  "Aldh1b1" = "#30013b",
  "Sfxn1"   = "#1b0022",
  "Gldc"    = "#070008")

p <- ggplot(expr_summary, aes(x = Group, y = ratio, group = Gene, color = Gene)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  theme_classic() +
  ylab("Expression Ratio (%)") +
  scale_color_manual(
    values = custom_colors,
    breaks = names(custom_colors)   
  ) +
  theme(
    axis.line = element_line(color = "black", size = 0.6),
    axis.ticks = element_line(color = "black"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),          
    axis.ticks.x = element_line(color = "black"),  
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "right",             
    legend.title = element_blank(),
    axis.text.y = element_text(size = 12, color = "black"),  
    axis.title.y = element_text(size = 14), 
    legend.text = element_text(size = 12),
  )
ggsave(filename = "../Output/Figure/Fig5B.png", plot = p, width = 8, height = 6,limitsize = FALSE)

#  Save variables
save(HSS_gene_sets,file = "../Output/Variables/SCP_HSS_gene_sets.R")







