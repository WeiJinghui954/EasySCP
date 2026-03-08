###########################
#### SCP Figure Code ####
###########################

#### -- Figure 3C -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "tidyverse",
  "dplyr",
  "Seurat",
  "SeuratWrappers",
  "ComplexHeatmap",
  "circlize",
  "ggplot2",
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

#  Read SCP liver raw data
file_path <- "../Data/SCP_Liver_Report.tsv"
dat1 <- read.delim(file_path, row.names = NULL)
dat1 <- dat1[dat1$PG.Genes != "", ] 
dat1$PG.Genes <- make.unique(dat1$PG.Genes) 
prot_cols <- grep("Quantity", colnames(dat1), value = TRUE)
prot_data <- dat1[, prot_cols] 
colnames(prot_data) <- sub(".*(\\d{8}).*(plate\\d+[_\\.][A-Za-z0-9]+).*", "\\1_\\2", colnames(prot_data))
colnames(prot_data) <- gsub("\\.", "_", colnames(prot_data))
prot_data <- as.data.frame(lapply(prot_data, function(x) as.numeric(as.character(x))))
prot_data[is.na(prot_data)] = 0
rownames(prot_data) = dat1$PG.Genes


#  Normalization
scaling_factor = 10000
df_normalized <- prot_data %>%
  mutate(across(where(is.numeric), ~ . / sum(.) * scaling_factor)) %>% 
  mutate(across(where(is.numeric), ~ log1p(.)))    
  
#  Screening proteins
sample_gates <- sapply(colnames(df_normalized), function(x) {sub("\\d+$", "", strsplit(x, split = "_")[[1]][3])}) 
names(sample_gates) <- colnames(df_normalized)
gated_proteins <- list()
gate_names <- unique(sample_gates)
for (gate in gate_names) {
  gate_samples <- which(sample_gates == gate)
  proteins_in_gate <- df_normalized[rowSums(df_normalized[, gate_samples] > 0) >= 6, ]
  proteins_in_gate <- proteins_in_gate[, gate_samples]
  gated_proteins[[gate]] <- proteins_in_gate
}
protein_names_by_gate <- lapply(gated_proteins, rownames)
common_proteins <- Reduce(union, protein_names_by_gate)
filtered_data <- df_normalized[rownames(df_normalized) %in% common_proteins, ]
colnames(filtered_data) <- sub(".*(\\d{8}).*(plate\\d+_[A-Za-z0-9]+).*", "\\1_\\2", colnames(filtered_data))

#  Extract cell information
samples <- colnames(filtered_data)
sample_info <- data.frame(samples = samples)
sample_info$labels <- gsub(".*_(.)[0-9]+.*", "\\1", samples)
sample_info$dates <- regmatches(samples, regexpr("\\d{8}", samples))
sample_info <- sample_info %>%
  mutate(plate = sub(".*_(plate\\d+_[A-Za-z0-9]+).*", "\\1", samples))

#  Create a Seurat object
disp_ratio <- apply(filtered_data, 1, var) / apply(filtered_data, 1, mean)
filtered_data_ratio <- filtered_data[disp_ratio > 0.1, ]

seurat_liver <- CreateSeuratObject(data =as.matrix(filtered_data_ratio), counts = as.matrix(filtered_data_ratio))
seurat_liver$sample_labels = sample_info$labels
seurat_liver$dates = sample_info$dates
seurat_liver <- subset(seurat_liver, dates != "20241220")

g2m.genes <- str_to_title(cc.genes$g2m.genes)
s.genes <- str_to_title(cc.genes$s.genes)

seurat_liver[["RNA"]] <- split(seurat_liver[["RNA"]], f = seurat_liver$dates)
seurat_liver <- CellCycleScoring(seurat_liver, s.features=s.genes, g2m.features=g2m.genes)
seurat_liver <- ScaleData(seurat_liver, vars.to.regress=c("S.Score","G2M.Score","nCount_RNA","nFeature_RNA"))

seurat_liver <- FindVariableFeatures(seurat_liver, nfeatures = 1200)
seurat_liver <- RunPCA(seurat_liver)

#  Batch Integration (FastMNN)
seurat_liver <- IntegrateLayers(
  object = seurat_liver, 
  method = FastMNNIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.mnn",
  k = 85,       
  ndist = 1.14,    
  verbose = FALSE
)

#  Dimensionality reduction visualization
seurat_liver <- RunUMAP(seurat_liver, reduction = "integrated.mnn", dims = 1:30, reduction.name = "umap.mnn")

#  Save plot
dim_plot <- DimPlot(
  seurat_liver,
  reduction = "umap.mnn",
  group.by = "sample_labels",
  cols = c( "#FDE725", "#B4DE2C", "#6DCD59", "#35B779",
            "#31688E", "#3E4A89", "#482878", "#440154"), 
  pt.size = 3                   # 调整点的大小（默认值为 1）  
)
umap_data <- dim_plot$data
label_mapping <- c(A = "1", B = "2", C = "3", D = "4", E = "5", F = "6", G = "7", H = "8")
umap_data$sample_labels <- as.character(label_mapping[umap_data$sample_labels])
umap_data$sample_labels <- factor(umap_data$sample_labels, levels = as.character(1:8))

umap_plot <-ggplot(umap_data, aes(x = umapmnn_1, y = umapmnn_2, fill = sample_labels)) +
  geom_point(shape = 21, color = "black", stroke = 0.5, size = 4) +  # 带黑色描边
  scale_fill_manual(values = c( "#FDE725", "#B4DE2C", "#6DCD59", "#35B779",
                                "#31688E", "#3E4A89", "#482878", "#440154")) +  
  theme_minimal() +
  labs( x = "UMAP 1", y = "UMAP 2") +
  theme(
    legend.position = "none",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 16)  
  )
umap_plot

ggsave(filename = "../Output/Figure/Fig3C.png", plot = umap_plot, dpi = 300, width = 6, height = 5,bg="white")

#  Save variables
save(filtered_data,file = "../Output/Variables/SCP_filtered_data.R")
save(sample_info,file = "../Output/Variables/SCP_sample_info.R")
save(seurat_liver,file = "../Output/Variables/SCP_seurat_liver.R")

