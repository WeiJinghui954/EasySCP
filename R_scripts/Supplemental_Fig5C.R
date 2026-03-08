###########################
#### SCP Figure Code ####
###########################

#### -- Supplemental Figure 5C -- ####

#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "ggplot2",
  "dplyr",
  "tibble",
  "ggrepel",
  "viridis",
  "clusterProfiler",
  "org.Mm.eg.db"
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
load("../Output/Variables/SCP_Guilliams_gene_protein.R")
load("../Output/Variables/SCP_Xu_gene_protein.R")

#  Group A and Group B
SCP_up_ST_down_genes <- rownames(df_zonated_only[
  df_zonated_only$SCP_bias > 0 &
    df_zonated_only$ST_10X_bias < 0, 
])
SCP_down_ST_up_genes <- rownames(df_zonated_only[
  df_zonated_only$SCP_bias < 0 &
    df_zonated_only$ST_10X_bias > 0, 
])
SCP_up_ST_down_genes_stereo <- rownames(df_zonated_only_stereo[
  df_zonated_only_stereo$Group == "SCP_up_ST_down",
])
SCP_down_ST_up_genes_stereo <- rownames(df_zonated_only_stereo[
  df_zonated_only_stereo$Group == "SCP_down_ST_up",
])

#  Intersection
Protein_down_RNA_up <- intersect(
  SCP_down_ST_up_genes_stereo,
  SCP_down_ST_up_genes
)
genes2_df <- bitr(
  Protein_down_RNA_up,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

#  GO enrichment
ego_up_rna <- enrichGO(
  gene = genes2_df$ENTREZID,
  OrgDb = org.Mm.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE
)

df2 <- as.data.frame(ego_up_rna)

#  Save plot
keywords <- c(
  "cytoplasmic translation",
  "glutathione metabolic process",
  "purine ribonucleotide metabolic process",
  "ribosome biogenesis",
  "proteasome-mediated ubiquitin-dependent protein catabolic process",
  "protein polyubiquitination",
  "mitochondrial respiratory chain complex I assembly",
  "oxidative phosphorylation",
  "cellular response to oxidative stress",
  "oxidoreductase complex"
)
filtered_data <- df2 %>%
  filter(Description %in% keywords) %>%
  arrange(desc(FoldEnrichment)) %>%             
  mutate(
    Description = factor(Description, levels = rev(Description)),  
    neglog10_padj = -log10(p.adjust)
  )

p <- ggplot(filtered_data, aes(x = FoldEnrichment, y = Description, fill = neglog10_padj)) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  scale_fill_viridis(option = "E", direction = 1) +  
  labs(
    x = "FoldEnrichment",
    y = NULL,
    fill = "-log10(p.adjust)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18,color = "black"),
    axis.text.x = element_text(size = 16,color = "black"),
    axis.title.x = element_text(size = 18, color = "black"),
    axis.title.y = element_text(size = 18, color = "black"),
    plot.title = element_text(size = 10,hjust = 0.5,color = "black"),
    panel.grid = element_blank(),  
    axis.line = element_line(color = "black"),  
    legend.title = element_text(size = 16, color = "black"),
    legend.text = element_text(size = 14)
  )
p
ggsave(filename =  "../Output/Figure/Supplemental Figure 5C.png", plot = p, width = 16, height = 8,bg="white")
