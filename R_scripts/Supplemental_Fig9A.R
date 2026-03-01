###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Fig9A -- ####

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

seurat_liver[["RNA"]] <- split(seurat_liver[["RNA"]], f = seurat_liver$dates)
seurat_liver <- NormalizeData(seurat_liver)
seurat_liver <- ScaleData(seurat_liver)
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

#label
dim_plot <- DimPlot(
  seurat_liver,
  reduction = "integrated.mnn",
  group.by = "sample_labels",
  cols = c( "#FDE725", "#B4DE2C", "#6DCD59", "#35B779",
            "#31688E", "#3E4A89", "#482878", "#440154"), 
  pt.size = 3                   
)
pca_data <- dim_plot$data
label_mapping <- c(A = "1", B = "2", C = "3", D = "4", E = "5", F = "6", G = "7", H = "8")
pca_data$sample_labels <- as.character(label_mapping[pca_data$sample_labels])
pca_data$sample_labels <- factor(pca_data$sample_labels, levels = as.character(1:8))

pca_plot <-ggplot(pca_data, aes(x = mnn_1, y = mnn_2, fill = sample_labels)) +
  geom_point(shape = 21, color = "black", stroke = 0.5, size = 4) +  # 带黑色描边
  scale_fill_manual(values = c( "#FDE725", "#B4DE2C", "#6DCD59", "#35B779",
                                "#31688E", "#3E4A89", "#482878", "#440154")) +  
  theme_minimal() +
  labs( x = "MNN 1", y = "MNN 2") +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),  
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    axis.text = element_text(size = 14),  
    axis.title = element_text(size = 16)  
  )
pca_plot
ggsave(filename = "../Output/Figure/Supplemental Fig9A.png", plot = pca_plot, dpi = 300, width = 6, height = 5,bg="white")
