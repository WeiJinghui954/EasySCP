###########################
#### SCP Figure Code ####
###########################

#### -- Figure 4D -- ####


#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "ggrepel",
  "stringr",
  "SeuratWrappers"
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
load("../Output/Variables/SCP_filtered_data.R")
load("../Output/Variables/SCP_sample_info.R")

#  Create a Seurat object
seurat_liver <- CreateSeuratObject(data =as.matrix(filtered_data), counts = as.matrix(filtered_data))

seurat_liver$sample_labels = sample_info$labels
seurat_liver$dates = sample_info$dates

seurat_liver <- subset(seurat_liver, dates != "20241220")
seurat_liver[["RNA"]] <- split(seurat_liver[["RNA"]], f = seurat_liver$dates)
seurat_liver <- ScaleData(seurat_liver)
seurat_liver <- FindVariableFeatures(seurat_liver, nfeatures = 1200)
head(seurat_liver[["RNA"]]@meta.data)
seurat_liver <- RunPCA(seurat_liver)

#  Batch Integration (FastMNN)
seurat_liver <- IntegrateLayers(
  object = seurat_liver, 
  method = FastMNNIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.mnn",
  k = 85,       
  ndist = 1.14,   
  verbose = FALSE)

#  Define groups
label_mapping <- setNames(as.character(1:8), LETTERS[1:8])
seurat_liver$sample_labels <- unname(label_mapping[seurat_liver$sample_labels])
seurat_liver$Group_CVPV <- ifelse(seurat_liver$sample_labels %in% c("1", "2"), "CV_associated", 
                                  ifelse(seurat_liver$sample_labels %in% c( "7", "8"), "PV_associated", NA))
Idents(seurat_liver) <- "Group_CVPV"

#  Differential analysis
seurat_liver <- JoinLayers(seurat_liver)
markers <- FindMarkers(seurat_liver, ident.1 = "PV_associated", ident.2 = "CV_associated",
                       logfc.threshold = 0, layer = "integrated.mnn")

markers$change <- ifelse(markers$p_val_adj < 0.05 & markers$avg_log2FC > 0.5, "PV_high_expression", 
                         ifelse(markers$p_val_adj < 0.05 & markers$avg_log2FC < -0.5, "CV_high_expression", "No_Change"))

#  Calculate the average expression level
avg_expr <- AverageExpression(seurat_liver, assays = "RNA", group.by = "ident")$RNA
colnames(avg_expr) <- gsub("-", "_", colnames(avg_expr)) 
avg_expr_selected <- avg_expr[, c("PV_associated", "CV_associated")]
markers$avg_expr_PV <- avg_expr_selected[rownames(markers), "PV_associated"]
markers$avg_expr_CV <- avg_expr_selected[rownames(markers), "CV_associated"]
markers_with_rownames <- cbind(Gene = rownames(markers), markers)
markers_with_rownames$log2_avg_expr_PV <- log2(markers_with_rownames$avg_expr_PV)
markers_with_rownames$log2_avg_expr_CV <- log2(markers_with_rownames$avg_expr_CV)


#  Save plot
CV_marker <- c("Gulo", "Cyp2e1", "Cyp2a5","Slc1a2","Slc22a1","Plpp3","Ces1","Rbp4","Idh3b","Cyp1a2","Scl10a1","Oat",
               "Wnt2","Rspo3","Plpp1","Tek","Hhip","Elovl5","Cyp4a12a","Acot3","Pdk4","Plod1","Sumf1")
PV_marker <- c("Hsd17b13","Pck1","Gldc","Ass1","Asl","Cyp2f2", "Ftcd","C8a","Gls2","Ctnnbip1","Xbp1",
               "Dll4","Ngfr","Notch3","Gfpt2","Pgm2","Gstp1","Dhrs9")
custom_genes <- union(CV_marker, PV_marker)

markers_with_rownames$Protein <-rownames(markers_with_rownames)
markers_with_rownames$label <- ifelse(markers_with_rownames$Protein %in% custom_genes, markers_with_rownames$Protein, NA)
markers_with_rownames$change <- factor(markers_with_rownames$change, levels = c("CV_high_expression", "PV_high_expression", "No_Change"))

markers_with_rownames$label_color <- ifelse(
  markers_with_rownames$change == "CV_high_expression", "#D4B900",  
  ifelse(markers_with_rownames$change == "PV_high_expression", "#440154", "grey50"))
markers_with_rownames$point_size <- ifelse(is.na(markers_with_rownames$label), 1, 1.5)
markers_with_rownames$point_color <- markers_with_rownames$change

p <- ggplot(markers_with_rownames, aes(x = log2_avg_expr_PV, y = log2_avg_expr_CV)) +
  geom_point(aes(color = point_color, size = point_size), alpha = 1) +
  geom_point(
    data = subset(markers_with_rownames, !is.na(label)),
    aes(x = log2_avg_expr_PV, y = log2_avg_expr_CV,
        color = I(label_color)),size = 1.5)+
  scale_color_manual(values = c(
    "CV_high_expression" = "#D4B900",
    "PV_high_expression" = "#440154",
    "No_Change" = "grey60"), labels = c(
      "CV_high_expression" = "CV",
      "PV_high_expression" = "PV",
      "No_Change" = "No Significant Change")) +
  scale_size_identity() +  
  labs(
    x = "Log2(average PV expression)",
    y = "Log2(average CV expression)",
    color = "Gene change") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    panel.background = element_rect(fill = "white", color = "black"),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-10, 6)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), limits = c(-10, 6))+
  ggrepel::geom_text_repel(
    aes(label = label),
    color = markers_with_rownames$label_color,
    size = 5,
    fontface = "bold",
    box.padding = 0.3,
    point.padding = 0.5,
    max.overlaps = 50)
ggsave(filename =  "../Output/Figure/Fig4D.png", plot = p, width = 6, height = 6,bg="white")

#  Save variables
save(markers_with_rownames,file = "../Output/Variables/SCP_CVPV_markers_with_rownames_rev.R")










