###########################
#### SCP Figure Code ####
###########################

#### -- Figure 4E -- ####


#  Prepare Workspace
rm(list = ls());gc()

packages <- c(
  "tidyverse",       
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
load("../Output/Variables/SCP_CVPV_markers_with_rownames_rev.R")

#  KEGG enrichment
Types = c("CV_high_expression", "PV_high_expression")
KEGG_results <- list()
for (type in Types) {
  significant_markers <- markers_with_rownames[markers_with_rownames$change == type, ]
  significant_markers$gene = rownames(significant_markers)
  
  df <- bitr(significant_markers$gene, fromType = "SYMBOL", 
             toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
  gene_diff <- df$ENTREZID
  kk_diff <- enrichKEGG(gene = gene_diff,
                        organism = "mmu",
                        qvalueCutoff = 0.5)
  KEGG_diff <- as.data.frame(kk_diff@result)
  KEGG_diff$geneID <- as.character(sapply(KEGG_diff$geneID, function(x) {
    paste(df$SYMBOL[match(strsplit(x, "/")[[1]], as.character(df$ENTREZID))], collapse = "/")
  }))
  KEGG_diff <- KEGG_diff %>%filter(p.adjust < 0.05)
  KEGG_results[[type]] <- KEGG_diff
}

all_KEGG_results <- list()
for (group in names(KEGG_results)) {
  KEGG_diff <- KEGG_results[[group]]
  kegg_data <- KEGG_diff[, c("subcategory","ID","Description","RichFactor", "FoldEnrichment","p.adjust","GeneRatio","pvalue","qvalue","Count","geneID")]
  kegg_data$Group <- group
  all_KEGG_results[[group]] <- kegg_data
}
common_data <- do.call(rbind, all_KEGG_results)
common_data$Description <- str_replace(common_data$Description, " - Mus musculus \\(house mouse\\)", "")

#  Save plot
keywords <- c("Peroxisome",
              "Biosynthesis of unsaturated fatty acids",
              "Linoleic acid metabolism",
              "PPAR signaling pathway",
              "Fatty acid degradation",
              "Drug metabolism - cytochrome P450",
              "ABC transporters",
              "Fatty acid elongation",
              
              "Glucagon signaling pathway",
              "Biosynthesis of amino acids",
              "Glycolysis / Gluconeogenesis",
              "Arginine biosynthesis",
              "Aminoacyl-tRNA biosynthesis",
              "Histidine metabolism",
              "Oxidative phosphorylation",
              "Ferroptosis",
              "Arginine and proline metabolism")  
filtered_data <- common_data[common_data$Description %in% keywords, ]

gene_sets <- setNames(
  lapply(seq_along(filtered_data$geneID), function(i) {
    list(
      genes = unlist(str_split(filtered_data$geneID[i], "/")),
      description = filtered_data$Description[i] 
    )
  }),filtered_data$ID )
filtered_data <- filtered_data %>%
  filter(!(Group == "CV_high_expression" & Description == "Ferroptosis")) %>%
  filter(!(Group == "PV_high_expression" & Description == "Drug metabolism - cytochrome P450"))
filtered_data <- filtered_data %>%
  mutate(FoldEnrichment = ifelse(Group == "CV_high_expression",
                                 -FoldEnrichment, FoldEnrichment)) %>%arrange(FoldEnrichment)
filtered_data$Group_Description <- paste(filtered_data$Description, filtered_data$Group, sep = " - ")
filtered_data$Group_Description <- factor(filtered_data$Group_Description, levels = unique(filtered_data$Group_Description))
filtered_data <- filtered_data %>%mutate(Description = sub("^.*/", "", Description))

p <- ggplot(filtered_data, aes(x = Group_Description, y = FoldEnrichment, fill = Group)) +
  geom_bar(stat = "identity", width = 0.8, color = "black") +
  scale_fill_manual(values = c("PV_high_expression" = "#440154", "CV_high_expression" = "#D4B900"),
                    labels = c("PV_high_expression" = "PV", "CV_high_expression" = "CV")) +
  coord_flip() +
  scale_x_discrete(labels = filtered_data$Description) + 
  labs(
    x = NULL,
    y = "FoldEnrichment",
    fill = "Group"
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
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )

ggsave(filename =  "../Output/Figure/Fig4E.png", plot = p, width = 12, height = 8,bg="white")















